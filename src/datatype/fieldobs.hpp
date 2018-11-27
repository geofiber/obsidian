//!
//! Contains the basic GDF datatypes for geological field observation sensors.
//!
//! \file datatype/fieldobs.hpp
//! \author Richard Scalzo
//! \date May, 2014
//! \license Affero General Public License version 3 or later
//! \copyright (c) 2018 NICTA + USyd
//!
/**
 *
 *
 * @file fieldobs.hpp
 * @author Richard Scalzo
 * @date 2018
 */

#pragma once

#include "datatype/forwardmodels.hpp"

namespace obsidian
{
  /**
   * Initial Parameters for FieldObs forward model. These parameters are required
   * only at the beginning of a set of simulations, to define the observation
   * locations and the prism boundaries the discretise the world
   */
  struct FieldObsSpec
  {
    // shape: nLocations x 2 (x, y)
    Eigen::MatrixXd locations;

    // Noise
    NoiseSpec noise;
  };

  struct FieldObsCache
  {
    std::vector<world::InterpolatorSpec> boundaryInterpolation;
    world::Query query;
  };

  /**
   * Update Parameters for initialised field observation forward model.
   */
  struct FieldObsParams
  {
    bool returnSensorData;
  };

  /**
   * Contains the results of a forward Seismic 1d simulation.
   */
  struct FieldObsResults
  {
    double likelihood;
    /**
     * Vector containing index of layer intersected at the surface
     */
    Eigen::VectorXd readings;
  };

  template<>
  struct Types<ForwardModel::FIELDOBS>
  {
    using Spec = FieldObsSpec;
    using Cache = FieldObsCache;
    using Params = FieldObsParams;
    using GlobalParams = WorldParams;
    using Results = FieldObsResults;
  };

  template<>
  inline std::string forwardModelLabel<ForwardModel::FIELDOBS>()
  {
    return "Field Observation";
  }

} // namespace obsidian

