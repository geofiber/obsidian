//!
//! Contains the implementation for the magnetic forward model.
//!
//! \file fwdmodel/magnetic.cpp
//! \author Darren Shen
//! \author Alistair Reid
//! \date 2014
//! \license Affero General Public License version 3 or later
//! \copyright (c) 2014, NICTA
//!

#include <glog/logging.h>
#include "magnetic.hpp"
#include "world/voxelise.hpp"

namespace obsidian
{
  namespace fwd
  {
    //! Generate a cache object for a magnetic forward model.
    //!
    //! \param boundaryInterpolation The world model interpolation parameters.
    //! \param worldSpec The world model specification.
    //! \param magSpec The forward model specification.
    //! \returns Forward model cache object.
    //!
    template<>
    MagCache generateCache<ForwardModel::MAGNETICS>(const std::vector<world::InterpolatorSpec>& boundaryInterpolation,
                                                    const WorldSpec& worldSpec, const MagSpec& magSpec)
    {
      LOG(INFO)<< "Caching mag sensitivity...";

      const VoxelSpec& magVox = magSpec.voxelisation;
      const Eigen::VectorXd& magB = magSpec.backgroundField;

      // Cache the query
      // Query MUST use internalGrid2DX :: X, then y, smaller to larger value
      const world::Query magQuery(boundaryInterpolation, worldSpec, magVox.xResolution, magVox.yResolution,
          magVox.zResolution, world::SamplingStrategy::noAA);

      GravmagInterpolatorParams interpParams = makeInterpParams(magVox, magSpec.locations, worldSpec);

      return// MagCache
      {
        boundaryInterpolation,
        magQuery,
        fwd::magSens(magQuery.edgeX, magQuery.edgeY,
            magQuery.edgeZ, interpParams.gridLocations, magB(0), magB(1), magB(2)),
        interpParams.sensorIndices,
        interpParams.sensorWeights
      };
    }

    //! Run a magnetic forward model.
    //! 
    //! \param spec The forward model specification.
    //! \param cache The forward model cache generated by generateCache().
    //! \param world The world model parameters.
    //! \returns Forward model results.
    //!
    template<>
    MagResults forwardModel<ForwardModel::MAGNETICS>(const MagSpec & spec, const MagCache& cache, const WorldParams& world)
    {
      Eigen::MatrixXd suscepts = world::getVoxels(cache.boundaryInterpolation, world, cache.query, RockProperty::Susceptibility);
      MagResults results;
      results.readings = fwd::computeField(cache.sensitivityMatrix, cache.sensorIndices, cache.sensorWeights, world::flatten(suscepts));
      return results;
    }

    namespace detail
    {
      //! Calculate the magnetic sensitivity at a particular position
      //! relative to the origin.
      //! 
      //! \param x, y, z The coordinates of the position.
      //! \param bx, by, bz The magnetic field at this position
      //!
      double magSensFunc(double x, double y, double z, double bx, double by, double bz)
      {
        // Compute the displacement
        double r = std::sqrt(x * x + y * y + z * z) + 1e-13;

        // Compute the normalisation factor for the magnetic field
        double normB = std::sqrt(bx * bx + by * by + bz * bz);

        return ((2 * by * bz * std::log(x + r)) + (2 * bz * bx * std::log(y + r)) + (2 * by * bx * std::log(z + r))
            + (bz * bz - by * by) * std::atan((x * z) / (y * r)) + (bz * bz - bx * bx) * std::atan((y * z) / (x * r))) / normB;
      }
    }

    Eigen::MatrixXd magSens(const Eigen::VectorXd &xEdges, const Eigen::VectorXd &yEdges, const Eigen::VectorXd &zEdges,
                            const Eigen::MatrixXd &locations, const double &bX, const double &bY, const double &bZ)
    {
      return magSens(xEdges, yEdges, zEdges, locations, Eigen::VectorXd::Constant(xEdges.rows(), bX),
                     Eigen::VectorXd::Constant(yEdges.rows(), bY), Eigen::VectorXd::Constant(zEdges.rows(), bZ));
    }

    Eigen::MatrixXd magSens(const Eigen::VectorXd &xEdges, const Eigen::VectorXd &yEdges, const Eigen::VectorXd &zEdges,
                            const Eigen::MatrixXd &locations, const Eigen::VectorXd &bX, const Eigen::VectorXd &bY,
                            const Eigen::VectorXd &bZ)
    {
      return detail::computeSensitivity(xEdges, yEdges, zEdges, bX, bY, -bZ, locations, detail::magSensFunc);
    }

  } // namespace fwd
} // namespace obsidian
