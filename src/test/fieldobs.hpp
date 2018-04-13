/**
 * Contains common testing utility functions for field observations
 *
 * @file fieldobs.hpp
 * @author Richard Scalzo
 * @date 2018-03-22
 * @license General Public License version 3 or later
 * @copyright (c) 2013, 2018 NICTA + USyd
 */

#pragma once

#include "datatype/datatypes.hpp"
#include "test/common.hpp"

namespace obsidian
{

  bool operator==(const FieldObsSpec& g, const FieldObsSpec& p)
  {
    return (g.locations == p.locations) && (g.noiseProb == p.noiseProb);
  }

  bool operator==(const FieldObsParams& g, const FieldObsParams& p)
  {
    return (g.returnSensorData == p.returnSensorData);
  }

  bool operator==(const FieldObsResults& g, const FieldObsResults& p)
  {
    return (g.likelihood == p.likelihood) && (g.readings == p.readings);
  }

  inline std::ostream& operator<<(std::ostream& os, const FieldObsSpec& spec)
  {
    os << "Field Observation SPEC" << "  LOCATIONS : " << spec.locations << " NOISE "
        << spec.noiseProb << std::endl;
    return os;
  }

  inline std::ostream& operator<<(std::ostream& os, const FieldObsResults& results)
  {
    os << "Field Observation RESULTS "<< "  READINGS " << results.readings << " LIKELIHOOD " << results.likelihood << std::endl;
    return os;
  }

  template<>
  inline void generateVariations<FieldObsSpec>(std::function<void(FieldObsSpec)> test)
  {
    for (uint l :
    { 0, 1, 5, 20 })
    {
      FieldObsSpec spec;
      spec.locations = testing::randomMatrix(l, 2);
      spec.noiseProb = testing::randomDouble();
      test(spec);
    }
  }

  template<>
  inline void generateVariations<FieldObsParams>(std::function<void(FieldObsParams)> test)
  {
    for (bool u :
    { true, false })
    {
      FieldObsParams param;
      param.returnSensorData = u;
      test(param);
    }
  }

  template<>
  inline void generateVariations<FieldObsResults>(std::function<void(FieldObsResults)> test)
  {
    for (uint u :
    { 0, 5, 20 })
    {
      FieldObsResults result;
      result.likelihood = testing::randomDouble();
      result.readings = testing::randomMatrix(u, 1);
      test(result);
    }
  }
}

