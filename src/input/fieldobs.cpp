//!
//! Input implementations related to field observations forward model.
//!
//! \file fieldobs.cpp
//! \author Richard Scalzo
//! \date 2018
//! \license General Public License version 3 or later
//! \copyright (c) 2014 + 2018, NICTA + USyd
//!

#include "common.hpp"

namespace obsidian
{

  template<>
  std::string configHeading<ForwardModel::FIELDOBS>()
  {
    return "fieldobs";
  }

  template<>
  void initSensorInputFileOptions<ForwardModel::FIELDOBS>(po::options_description & options)
  {
    options.add_options() //
    ("fieldobs.enabled", po::value<bool>(), "enable sensor") //
    ("fieldobs.sensorLocations", po::value<std::string>(), "sensor locations") //
    ("fieldobs.sensorReadings", po::value<std::string>(), "sensor readings") //
    ("fieldobs.noiseProb", po::value<double>(), "noise probability that field obs is wrong") //
  }

  template<>
  FieldObsSpec parseSpec(const po::variables_map& vm, const std::set<ForwardModel> & sensorsEnabled)
  {
    FieldObsSpec spec;
    if (sensorsEnabled.count(ForwardModel::FIELDOBS))
    {
      std::vector<std::vector<std::string>> data = io::csv::readRaw(vm["fieldobs.sensorLocations"].as<std::string>());
      spec.locations = readFixed<double, 2, 0>(data);
      spec.noiseProb = vm["fieldobs.noiseProb"].as<double>();
    }
    return spec;
  }

  template<>
  po::variables_map write<>(const std::string & prefix, FieldObsSpec spec, const po::options_description & od)
  {
    std::vector<std::vector<std::string>> data;
    writeFixed<double>(data, spec.locations);
    io::csv::writeRaw(prefix + "sensorLocations.csv", data);

    return build_vm(po::variables_map(), od, "fieldobs",
                    { { "sensorLocations", prefix + "sensorLocations.csv" },
                      { "noiseProb", io::to_string(spec.noiseProb) } });
  }

  //! @note the sensor params don't actually have anything in them at the moment so we don't need to do any parsing
  //!
  template<>
  FieldObsParams parseSimulationParams(const po::variables_map& vm, const std::set<ForwardModel> & sensorsEnabled)
  {
    return FieldObsParams();
  }

  template<>
  FieldObsResults parseSensorReadings(const po::variables_map& vm, const std::set<ForwardModel> & sensorsEnabled)
  {
    FieldObsResults s;
    if (sensorsEnabled.count(ForwardModel::FIELDOBS))
    {
      std::vector<std::vector<std::string>> data = io::csv::readRaw(vm["fieldobs.sensorReadings"].as<std::string>());
      s.readings = readRagged<uint, 0>(data);
      s.likelihood = 0.0;
    }
    return s;
  }
  template<>
  po::variables_map write<>(const std::string & prefix, FieldObsResults g, const po::options_description & od)
  {
    std::vector<std::vector<std::string>> data;
    writeRagged<double>(data, g.readings);
    io::csv::writeRaw(prefix + "sensorReadings.csv", data);
    return build_vm(po::variables_map(), od, "fieldobs", { { "sensorReadings", prefix + "sensorReadings.csv" } });
  }

  template<>
  void enableProperties<ForwardModel::FIELDOBS>(Eigen::VectorXi & propMasksMask)
  {
  }

  template<>
  prior::FieldObsParamsPrior parsePrior(const po::variables_map& vm, const std::set<ForwardModel> & sensorsEnabled)
  {
    return prior::FieldObsParamsPrior();
  }

  template<>
  bool validateSensor(const WorldSpec & world, const FieldObsSpec & spec, const FieldObsResults & result)
  {
    bool valid = true;
    if (spec.locations.rows() == 0)
    {
      LOG(ERROR)<< "input: no field observation locations specified. Disable forward model if it is not used.";
      valid = false;
    }
    if (static_cast<uint>(spec.locations.cols()) != 2)
    {
      LOG(ERROR)<< "input: locations in contact point must have two cols (x, y).";
      valid = false;
    }
    for (uint l = 0; l < spec.locations.rows(); l++)
    {
      if (spec.locations(l, 0) < world.xBounds.first || spec.locations(l, 0) > world.xBounds.second
          || spec.locations(l, 1) < world.yBounds.first || spec.locations(l, 1) > world.yBounds.second)
      {
        LOG(ERROR)<< "input: contact point location " << (l + 1) << " is out of world bounds";
        valid = false;
      }
    }
    if (spec.noiseProb <= 0 || spec.noiseProb >= 1)
    {
      LOG(ERROR)<< "input: field observation noise probability must be between 0 and 1";
      valid = false;
    }
    if (static_cast<uint>(spec.locations.rows()) != result.readings.size())
    {
      LOG(ERROR)<< "input: different number of readings for field observation results (" << result.readings.size() << ") to location specified. Remove or add more locations ("<< spec.locations.rows() <<").";
      valid = false;
    }
    for (uint l = 0; l < result.readings.size(); l++)
    {
      if (result.readings[l] < 0 or result.readings[l] >= world.boundaries.size())
      {
        LOG(ERROR)<< "input: field observation must be a valid world boundary index between 0 and" << world.boundaries.size();
        valid = false;
      }
    }
    return valid;
  }
}
