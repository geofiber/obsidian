//!
//! Test field observation 2d sensor forward model.
//!
//! \file fwdmodel/testfieldobs.cpp
//! \author Richard Scalzo
//! \date 2018
//! \license Affero General Public License version 3 or later
//! \copyright (c) 2014, 2018 NICTA + USyd
//!

#include <gtest/gtest.h>

#include "datatype/datatypes.hpp"
#include "world/interpolate.hpp"
#include "fwdmodel/fwd.hpp"
#include "test/world.hpp"

namespace obsidian
{
  class SensorTest: public ::testing::Test
  {
  public:
    std::vector<WorldSpec> worldSpecs;
    std::vector<WorldParams> worldParams;
    std::vector<FieldObsSpec> sensorSpec;
    std::vector<FieldObsParams> sensorParams;
    std::vector<FieldObsResults> sensorResults;

    virtual void SetUp()
    {
      WorldSpec worldSpec_;
      WorldParams worldParams_;

      testing::initWorld(worldSpec_, worldParams_, -10, 10, 20, -10, 10, 20, 0, 20, 5, [](double x, double y, uint boundary)
      {
        return boundary;
      },
                         [](double x, double y, double boundary)
                         {
                           return 0;
                         },
                         [](uint layer, uint property)
                         {
                           if (property == static_cast<uint>(RockProperty::PWaveVelocity))
                           {
                             return 1.0;
                           }
                           return 0.0;
                         });

      worldSpecs.push_back(worldSpec_);
      worldParams.push_back(worldParams_);

      FieldObsSpec spec;
      spec.locations.resize(2, 3);
      for (uint sensor = 0; sensor < 2; sensor++)
      {
        spec.locations(sensor, 0) = sensor * sin(sensor / 5.0 * 2.0 * M_PI);
        spec.locations(sensor, 1) = sensor * cos(sensor / 5.0 * 2.0 * M_PI);
        spec.locations(sensor, 2) = 0;
        Eigen::VectorXi interface(4);
      }
      sensorSpec.push_back(spec);
      FieldObsParams params;
      params.returnSensorData = true;
      sensorParams.push_back(params);
      FieldObsResults results;
      Eigen::VectorXd readings(4);
      readings << 0, 1, 2, 3;
      results.readings.push_back(readings);
      results.readings.push_back(readings);
      sensorResults.push_back(results);
    }
  };

TEST_F(SensorTest, test)
{
  for (size_t world = 0; world < worldSpecs.size(); world++)
  {
    WorldSpec worldSpec = worldSpecs[world];
    std::vector<world::InterpolatorSpec> layerInterp = world::worldspec2Interp(worldSpec);

    WorldParams worldParam = worldParams[world];

    for (uint sensor = 0; sensor < sensorSpec.size(); sensor++)
    {
      FieldObsSpec spec = sensorSpec[sensor];
      FieldObsCache cache = fwd::generateCache<ForwardModel::FIELDOBS>(layerInterp, worldSpec, spec);
      FieldObsResults results = fwd::forwardModel<ForwardModel::FIELDOBS>(spec, cache, worldParam);
      FieldObsResults expected = sensorResults[sensor];
      double error = 0;
      for (uint r = 0; r < results.readings.size(); r++)
      {
        error += (results.readings[r] - expected.readings[r]).norm();
      }
      EXPECT_LT(error, 0.1);
    }
  }
}
}

int main(int argc, char **argv)
{
::testing::InitGoogleTest(&argc, argv);
return RUN_ALL_TESTS();
}
