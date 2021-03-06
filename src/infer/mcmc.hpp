//!
//! Contains the implementation of Monte Carlo Markov Chain simulations.
//!
//! \file infer/mcmc.cpp
//! \author Lachlan McCalman
//! \author Darren Shen
//! \date 2014
//! \license Affero General Public License version 3 or later
//! \copyright (c) 2014, NICTA
//!

#pragma once

#include <limits>
#include <random>
#include <functional>
#include <iomanip>
#include <chrono>
#include <Eigen/Dense>
#include <boost/circular_buffer.hpp>
#include <iostream>
#include <glog/logging.h>

#include "db/db.hpp"
#include "app/settings.hpp"
#include "infer/chainarray.hpp"
#include "infer/diagnostics.hpp"
#include "distrib/multigaussian.hpp"
#include "comms/transport.hpp"

namespace ph = std::placeholders;

namespace stateline
{
  namespace mcmc
  {
    //! Represents a Markov Chain Monte Carlo sampler that returns samples from
    //! a distribution.
    //!
    class Sampler
    {
    public:
      //! Create a new sampler.
      //!
      //! \param s Settings to tune various parameters of the MCMC.
      //! \param d Settings for configuring the database for storing samples.
      //! \param stateDim The number of dimensions in each sample.
      //! \param interrupted A flag used to monitor whether the sampler has been interrupted.
      //!
      Sampler(const MCMCSettings& s, const DBSettings& d, uint stateDim, volatile bool& interrupted)
          : db_(d),
            chains_(s.stacks, s.chains, s.initialTempFactor, s.proposalInitialSigma, s.initialSigmaFactor, db_, s.cacheLength, d.recover),
            lengths_(s.stacks * s.chains, 0),
            propStates_(s.stacks * s.chains, stateDim),
            locked_(s.stacks * s.chains, false),
            nextChainBeta_(s.stacks * s.chains),
            nAcceptsGlobal_(s.stacks * s.chains, 0),
            nSwapsGlobal_(s.stacks * (s.chains), 0),
            nSwapAttemptsGlobal_(s.stacks * s.chains, 0),
            sigmas_(s.stacks * s.chains),
            betas_(s.stacks * s.chains),
            acceptRates_(s.stacks * s.chains),
            swapRates_(s.stacks * s.chains),
            lowestEnergies_(s.stacks * s.chains),
            SX_(s.stacks * s.chains),
            SX2_(s.stacks * s.chains),
            qcovL_(s.stacks * s.chains),
            s_(s),
            recover_(d.recover),
            numOutstandingJobs_(0),
            interrupted_(interrupted),
            context_(1)
      {
        // Initialise for logging purposes
        for (uint i = 0; i < s.chains * s.stacks; i++)
        {
          lengths_[i] = chains_.length(i);
          sigmas_[i] = chains_.sigma(i);
          betas_[i] = chains_.beta(i);
          nextChainBeta_[i] = chains_.beta(i);
          acceptRates_[i] = 1;
          swapRates_[i] = 0;
          lowestEnergies_[i] = std::numeric_limits<double>::infinity();
          acceptBuffers_.push_back(boost::circular_buffer<bool>(s.adaptionLength));
          swapBuffers_.push_back(boost::circular_buffer<bool>(s.adaptionLength / s.swapInterval + 1));
          nAcceptsGlobal_[i] = 1;
          acceptBuffers_[i].push_back(true); // first state always accepts
          nSwapsGlobal_[i] = 0;
          nSwapAttemptsGlobal_[i] = 1;
          swapBuffers_[i].push_back(false); // gets rid of a nan, not really needed
        }
      }

      //! Run the sampler for a period of time.
      //!
      //! \param policy Async policy to evaluate states.
      //! \param initialStates Initial chain states. Ignored if recovering.
      //! \param propFn The proposal function.
      //! \param propPDFFn The proposal PDF function.
      //! \param numSeconds The number of seconds to run the MCMC for.
      //!
      template<class AsyncPolicy, class PropFn, class PropPdfFn>
      void run(AsyncPolicy &policy, const std::vector<Eigen::VectorXd>& initialStates, PropFn &propFn, PropPdfFn &propPdfFn, uint numSeconds)
      {
        using namespace std::chrono;

        // Used for publishing statistics to visualisation server.
        // RS 2018/09/12:  Since we're not using such a server, comment out
        // to enable us to run multiple simulations at once.
        // zmq::socket_t publisher(context_, ZMQ_PUB);
        // publisher.bind("tcp://*:5556");

        // Record the starting time of the MCMC
        steady_clock::time_point startTime = steady_clock::now();

        // Initialise the chains if we're not recovering
        if (!recover_)
        {
          initialise(policy, initialStates);
        }

        // Start all the chains from hottest to coldest
        for (uint i = 0; i < chains_.numTotalChains(); i++)
        {
          uint c = chains_.numTotalChains() - i - 1;
          propose(policy, c, propFn);
        }

        // Initialise the convergence criteria
        uint stateDim = initialStates[0].size();
        EpsrConvergenceCriteria cc(chains_.numStacks(), stateDim);

        // RS 2018/03/09:  Initialize the chain means and covariances.
        for (uint i = 0; i < chains_.numTotalChains(); i++)
        {
          SX_[i] = Eigen::VectorXd::Zero(stateDim);
          SX2_[i] = Eigen::MatrixXd::Zero(stateDim, stateDim);
          qcovL_[i] = Eigen::MatrixXd::Identity(stateDim, stateDim);
        }

        // Listen for replies. As soon as a new state comes back,
        // add it to the corresponding chain, and submit a new proposed state
        auto lastLogTime = steady_clock::now();
        auto lastPrintTime = steady_clock::now();
        while (duration_cast<seconds>(steady_clock::now() - startTime).count() < numSeconds && !interrupted_)
        {
          std::pair<uint, double> result;

          // Wait a for reply
          try
          {
            result = policy.retrieve();
          }
          catch (...)
          {
            VLOG(3) << "Comms error -- probably shutting down";
          }

          // If we weree interrupted this state will be garbage
          if (interrupted_)
            break;

          numOutstandingJobs_--;
          uint id = result.first;
          double energy = result.second;

          // Check if this chain is either the coldest or the hottest
          bool isHottestChainInStack = id % chains_.numChains() == chains_.numChains() - 1;
          bool isColdestChainInStack = id % chains_.numChains() == 0;

          // Handle the new proposal and add a new state to the chain
	      double proposalDensity = propPdfFn(propStates_.row(id), sigmas_[id]);
	      double lastsampleDensity = propPdfFn(chains_.lastState(id).sample, sigmas_[id]);
	      double logDensityRatio = proposalDensity - lastsampleDensity;
          VLOG(3) << "logDensityRatio = " << logDensityRatio;
          State propState { propStates_.row(id), energy, logDensityRatio, chains_.beta(id), false, SwapType::NoAttempt };
          bool propAccepted = chains_.append(id, propState);
          lengths_[id] += 1;
          updateAccepts(id, propAccepted);

          // Update the convergence test if this is the coldest chain in a stack
          if (isColdestChainInStack && chains_.numChains() > 1 && chains_.numStacks() > 1)
          {
            cc.update(id / chains_.numChains(), chains_.lastState(id).sample);
          }

          // Check if this chain was locked. If it was locked, it means that
          // the chain above (hotter) locked it so that it can swap with it
          if (locked_[id])
          {
            // Try swapping this chain with the one above it
            bool swapAccepted = chains_.swap(id, id + 1);
            updateSwaps(id + 1, swapAccepted);
            // Unlock this chain, propgating the lock downwards
            unlock(policy, id, propFn);
          }
          else if (isHottestChainInStack && lengths_[id] % s_.swapInterval == 0 && chains_.numChains() > 1)
          {
            // The hottest chain is ready to swap. Lock the next chain
            // to prevent it from proposing any more
            locked_[id - 1] = true;
          }
          else
          {
            // This chain is not locked, so we can propose
            try
            {
              propose(policy, id, propFn);
            }
            catch (...)
            {
              VLOG(3) << "Comms error -- probably shutting down";
            }
          }

          // Check again after a new interaction with comms
          if (interrupted_)
            break;

          // Log the best energy state so far
          lowestEnergies_[id] = std::min(lowestEnergies_[id], chains_.lastState(id).energy);

          // RS 2018/03/09:  Update the chain covariance.
          updateChainsums(id);

          // Check if we need to adapt the step size for this chain
          if (lengths_[id] % s_.proposalAdaptInterval == 0)
          {
            adaptSigma(id);
          }

          // Update the temperature which might have changed while waiting
          chains_.setBeta(id, nextChainBeta_[id]);
          betas_[id] = nextChainBeta_[id];
          // Check for adapting the temperatures of all the chains but the 1st
          if (lengths_[id] % s_.betaAdaptInterval == 0 && !isColdestChainInStack)
          {
            adaptBeta(id);
          }

          // Update the accept and swap rates
          // RS 2018/09/18:  now writes out every 2 sec (instead of 50 msec),
          // to keep log files manageable
          if (duration_cast<milliseconds>(steady_clock::now() - lastLogTime).count() > 2000)
          {
            lastLogTime = steady_clock::now();

            std::stringstream s;
            s << "\n\nChainID  Length  MinEngy  CurrEngy    Sigma      AcptRt    GlbAcptRt    Beta     SwapRt   GlbSwapRt\n";
            s << "-----------------------------------------------------------------------------------------------------\n";
            for (uint i = 0; i < chains_.numTotalChains(); i++)
            {
              if (i % chains_.numChains() == 0 && i != 0)
                s << '\n';
              s << std::setprecision(6) << std::showpoint << i << " " << std::setw(9) << lengths_[i] << " " << std::setw(10)
                  << lowestEnergies_[i] << " " << std::setw(10) << chains_.lastState(i).energy << " " << std::setw(10) << sigmas_[i] << " "
                  << std::setw(10) << acceptRates_[i] << " " << std::setw(10) << nAcceptsGlobal_[i] / (double) lengths_[i] << " "
                  << std::setw(10) << betas_[i] << " " << std::setw(10) << swapRates_[i] << " " << std::setw(10)
                  << nSwapsGlobal_[i] / (double) nSwapAttemptsGlobal_[i] << " \n";
            }

            // Quick and dirty way to get the data to the visualisation server
            // RS 2018/09/12:  commented out to enable multiple concurrent runs.
            // comms::sendString(publisher, s.str());

            if (duration_cast<milliseconds>(steady_clock::now() - lastPrintTime).count() > 2000)
            {
              lastPrintTime = steady_clock::now();

              LOG(INFO)<< s.str() << "\n";

              if (chains_.numStacks() > 1 && chains_.numChains() > 1)
              {
                if (cc.hasConverged())
                {
                  LOG(INFO)<< "converged: true ("<< cc.rHat().transpose() << " < 1.1)";
                }
                else
                {
                  LOG(INFO) << "converged: false ("<< cc.rHat().transpose() << " > 1.1)";
                }
              }
            }
          }
        }

        // Time limit reached. We need to now retrieve all outstanding job results.
        if (!interrupted_)
        {
          while (numOutstandingJobs_--)
          {
            auto result = policy.retrieve();
            uint id = result.first;
            double energy = result.second;
	    double proposalDensity = propPdfFn(propStates_.row(id), sigmas_[id]);
	    double lastsampleDensity = propPdfFn(chains_.lastState(id).sample, sigmas_[id]);
	    double logDensityRatio = proposalDensity - lastsampleDensity;
            State propState { propStates_.row(id), energy, logDensityRatio, chains_.beta(id), false, SwapType::NoAttempt };
            bool propAccepted = chains_.append(id, propState);
            lengths_[id] += 1;
            updateAccepts(id, propAccepted);
          }
        }
        else
        {
          LOG(INFO)<< "Inference interrupted by user: Exiting.";
        }

        // Manually flush any chain states that are in memory to disk
        for (uint i = 0; i < chains_.numTotalChains(); i++)
        {
          chains_.flushCache(i);
        }

        if (cc.hasConverged())
        {
          LOG(INFO)<< "MCMC has converged";
        }
        else
        {
          LOG(INFO) << "WARNING: MCMC has not converged";
        }

        LOG(INFO)<<"Length of chain 0: " << lengths_[0];
    }

    //! Get the MCMC chain array.
    //!
    //! \return A copy of the chain array.
    //!
    ChainArray chains()
    {
      return chains_;
    }

    //! Get the MCMC chain array.
    //!
    //! \return A const reference to the chain array.
    //!
    const ChainArray &chains() const
    {
      return chains_;
    }

    //! RS 2018/08/10:  Get the parameter mean of MCMC chain id.
    //! Relevant mostly for multivariate adaptive Metropolis proposals.
    //!
    //! \return The chain mean.
    //!
    const Eigen::VectorXd chainmean(uint id)
    {
      uint n = lengths_[id];
      return ((n < 1) ? 0.0 : 1.0/n) * SX_[id];
    }

    //! RS 2018/08/10:  Get the parameter covariance of MCMC chain id.
    //! Relevant mostly for multivariate adaptive Metropolis proposals.
    //!
    //! \return The chain covariance.
    //!
    const Eigen::MatrixXd chaincov(uint id) const
    {
      int stateDim = SX2_[id].rows();
      uint n = lengths_[id];
      Eigen::MatrixXd eid = 1.0e-10*Eigen::MatrixXd::Identity(stateDim, stateDim);
      return (n < 2) ? eid : (SX2_[id] - SX_[id]*SX_[id].transpose()/n)/(n-1) + eid;
    }

  private:

    //! Initialise the sampler.
    //!
    //! \param policy Async policy to evaluate states.
    //! \param initialStates A list containing the initial states of the chains.
    //!
    template <class AsyncPolicy>
    void initialise(AsyncPolicy &policy, const std::vector<Eigen::VectorXd>& initialStates)
    {
      // Evaluate the initial states of the chains
      for (uint i = 0; i < chains_.numTotalChains(); i++)
      {
        propStates_.row(i) = initialStates[i];
        policy.submit(i, propStates_.row(i));
      }

      // Retrieve the energies and temperatures for the initial states
      for (uint i = 0; i < chains_.numTotalChains(); i++)
      {
        auto result = policy.retrieve();
        uint id = result.first;
        double energy = result.second;
        double logDensityRatio = 1.0;
        State s
        { initialStates[id], energy, logDensityRatio, chains_.beta(id), true, SwapType::NoAttempt};
        chains_.initialise(id, s);
      }
    }

    //! Propose a new state through an async policy.
    //!
    //! \param policy Async policy to evaluate proposals.
    //! \param id The id of the chain that is proposing.
    //! \param propFn The proposal function.
    //!
    template <class AsyncPolicy, class PropFn>
    void propose(AsyncPolicy &policy, uint id, PropFn &propFn)
    {
      //! RS 2018/03/12:  If we're doing adaptive Metropolis, change over
      //! to the AM proposal function if we've got at least some target
      //! number of samples in this chain.
      int stateDim = chains_.lastState(id).sample.size();
      if (qcovL_[id].rows() == 0)
      {
        VLOG(2) << "qcov = NULL; reverting to I(" << stateDim << "," << stateDim << ")";
        qcovL_[id] = Eigen::MatrixXd::Identity(stateDim, stateDim);
      }
      propStates_.row(id) = propFn(chains_.lastState(id).sample, chains_.sigma(id), qcovL_[id]);
      policy.submit(id, propStates_.row(id));
      numOutstandingJobs_++;
    }

    //! Unlock a chain and reactivate any chains that were waiting for it.
    //!
    //! \param policy Async policy to evulate proposals.
    //! \param id The id of the chain to unlock.
    //! \param propFn The proposal function.
    //!
    template <class AsyncPolicy, class PropFn>
    void unlock(AsyncPolicy &policy, uint id, PropFn &propFn)
    {
      // Unlock this chain
      locked_[id] = false;

      // The hotter chain no longer has to wait for this chain, so
      // it can propose new state
      propose(policy, id + 1, propFn);

      // Check if this was the coldest chain
      if (id % chains_.numChains() != 0)
      {
        // Lock the chain that is below (colder) than this.
        locked_[id - 1] = true;
      }
      else
      {
        // This is the coldest chain and there is no one to swap with
        propose(policy, id, propFn);
      }
    }

    //! Return a new proposal step size for a chain.
    //!
    //! \param id The id of the chain to be adapted.
    //! \return The new proposal step size.
    //!
    double adaptSigma(uint id)
    {
      double acceptRate = acceptRates_[id];
      double oldSigma= chains_.sigma(id);
      double factor = std::pow(acceptRate / s_.proposalOptimalAccept, s_.proposalAdaptRate);
      double boundFactor = std::min(std::max(factor, s_.proposalMinFactor), s_.proposalMaxFactor);
      double gamma = s_.adaptionLength/(double)(s_.adaptionLength+lengths_[id]);

      double newSigma = oldSigma * std::pow(boundFactor, gamma);
      VLOG(2) << "Adapting Sigma" << id <<":" << oldSigma << "->" << newSigma << " @acceptrate:" << acceptRate;
      chains_.setSigma(id, newSigma);
      sigmas_[id] = newSigma;

      // Ensure higher temperature chains have larger sigmas than chains colder than it
      if (id % chains_.numChains() != 0 && newSigma < chains_.sigma(id - 1))
      {
        newSigma = chains_.sigma(id - 1);
      }

      return newSigma;
    }

    //! Adapt a new temperature for a chain.
    //!
    //! \param id The id of the chain to be adapted.
    //!
    void adaptBeta(uint id)
    {
      // Adapt the temperature
      double swapRate = swapRates_[id];
      double rawFactor = std::pow(swapRate/s_.betaOptimalSwapRate, s_.betaAdaptRate);
      double boundedFactor = std::min( std::max(rawFactor, s_.betaMinFactor), s_.betaMaxFactor);
      double beta = chains_.beta(id);
      double lowerBeta = chains_.beta(id-1);// temperature changes propogate UP
      double factor = 1.0/std::max(boundedFactor, 2*beta/(beta + lowerBeta));

      // Set the temperature for this chain (because it hasn't proposed yet)
      double gamma = s_.adaptionLength/(double)(s_.adaptionLength+lengths_[id]);
      double newbeta = chains_.beta(id) * std::pow(factor, gamma);
      chains_.setBeta(id,newbeta);
      nextChainBeta_[id] = newbeta;
      betas_[id] = newbeta;

      // Just for logging
      VLOG(2) << "Adapting Beta" << id << ":" << beta << "->" << newbeta << " @swaprate:" << swapRate;

      // Loop through the other temperatures
      uint coldestChainId = (uint)(id / (double)chains_.numChains()) * chains_.numChains();
      uint hottestChainId = coldestChainId + chains_.numChains()-1;

      // Set the next temperatures for the other chains (as they're busy)
      for (uint i = id+1; i <= hottestChainId; i++)
      {
        nextChainBeta_[i] = nextChainBeta_[i] * std::pow(factor, gamma);
      }
    }

    //! RS 2018/03/09:  Update the chain covariances.
    //!
    //! \param id The id of the chain to be updated.
    //!
    void updateChainsums(uint id)
    {
      uint n = lengths_[id];
      if (n == 0) return;

      // Update sums used to make the covariance
      VLOG(3) << "updating chainsums";
      Eigen::VectorXd Xk = chains_.lastState(id).sample;
      SX_[id] += Xk;
      SX2_[id] += Xk*Xk.transpose();

      // Build *correlation* (whitened covariance) matrix of chain
      uint amL = s_.adaptAMLength;
      if (n > 0 && n % amL == 0)
      {
        int stateDim = Xk.size();
        qcovL_[id] = (n*chaincov(id) + amL*Eigen::MatrixXd::Identity(stateDim, stateDim)) / (1.0*(n + amL));
        qcovL_[id] = qcovL_[id].llt().matrixL();
      }

      // Regularize with round Gaussian so that the proposal transitions
      // smoothly from isotropic MHRW as it gains samples
      VLOG(3) << "done updating chainsums";
    }

    void updateAccepts(uint id, bool acc)
    {
      uint oldSize = acceptBuffers_[id].size();
      double oldRate = acceptRates_[id];
      bool isFull = acceptBuffers_[id].full();
      bool lastAcc = acceptBuffers_[id].front();

      // Now push on the new state
      acceptBuffers_[id].push_back(acc);

      // Compute the new rate
      uint newSize = acceptBuffers_[id].size();
      double delta = ((int)acc - (int)(lastAcc&&isFull))/(double)newSize;
      double scale = oldSize/(double)newSize;
      acceptRates_[id] = std::max(oldRate*scale + delta, 0.0);
      nAcceptsGlobal_[id] += acc;
      //double comparison = 1.0 + (0.5/newSize);
      double comparison = 1.0;
      if (acceptRates_[id] > comparison)
      {
        std::cout << "oldSize: " << oldSize << "\n"
        << "isFull:" << isFull << "\n"
        << "newSize:" << newSize << "\n"
        << "lastAcc:" << lastAcc << "\n"
        << "delta:" << delta << "\n"
        << "scale:" << scale << "\n"
        << "oldRate:" << oldRate << "\n"
        << "rate:" << acceptRates_[id] << "\n"
        << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    void updateSwaps(uint id, bool sw)
    {
      uint oldSize = swapBuffers_[id].size();
      double oldRate = swapRates_[id];
      bool isFull = swapBuffers_[id].full();
      bool lastSw = swapBuffers_[id].front();

      // Now push back the new state
      swapBuffers_[id].push_back(sw);

      // Compute the new rate
      uint newSize = swapBuffers_[id].size();
      double delta = ((int)sw - (int)(lastSw&&isFull))/(double)newSize;
      double scale = oldSize/(double)newSize;
      swapRates_[id] = std::max(oldRate*scale + delta, 0.0);
      nSwapsGlobal_[id] += sw;// for global rate
      nSwapAttemptsGlobal_[id] += 1;
    }

    // The MCMC chain database
    db::Database db_;

    // The MCMC chain wrapper
    ChainArray chains_;

    // lengths of the chains
    std::vector<uint> lengths_;

    // Matrix of proposed states
    Eigen::MatrixXd propStates_;

    // Whether a chain is locked. A locked chain will wait for any outstanding
    // job results and propagate the lock.
    std::vector<bool> locked_;

    // Cache the next temperature as computed by lower chains in the stack
    std::vector<double> nextChainBeta_;

    // Keep track of the swaps and accepts for adaption
    std::vector<unsigned long long> nAcceptsGlobal_;
    std::vector<unsigned long long> nSwapsGlobal_;
    std::vector<unsigned long long> nSwapAttemptsGlobal_;

    std::vector<boost::circular_buffer<bool>> acceptBuffers_;
    std::vector<boost::circular_buffer<bool>> swapBuffers_;

    // For logging purposes
    std::vector<double> sigmas_;
    std::vector<double> betas_;
    std::vector<double> acceptRates_;
    std::vector<double> swapRates_;
    std::vector<double> lowestEnergies_;

    // RS 2018/08/10:  Accumulated sums for chain mean and covariance.
    std::vector<Eigen::VectorXd> SX_;
    std::vector<Eigen::MatrixXd> SX2_;

    // RS 2018/08/10:  Adaptive multivariate MHRW proposal covariance.
    // RS 2018/08/23:  Now represented by its Cholesky factor.
    std::vector<Eigen::MatrixXd> qcovL_;

    // The MCMC settings
    MCMCSettings s_;

    // Recovering?
    bool recover_;

    // How many jobs haven't been retrieved?
    uint numOutstandingJobs_;

    // Whether an interrupt signal has been sent to the sampler.
    volatile bool& interrupted_;

    zmq::context_t context_;
  };
}
// namespace mcmc
}// namespace obsidian
