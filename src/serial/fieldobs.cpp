//!
//! \file serial/fieldobs.cpp
//! \author Richard Scalzo
//! \date March, 2018
//! \license Affero General Public License version 3 or later
//! \copyright (c) 2014, 2018 NICTA + USyd
//!

#include "serial/seismic.hpp"
#include "serial/serialtypes.pb.h"
#include "serial/utility.hpp"

namespace obsidian
{
  namespace comms
  {

    std::string serialise(const FieldObsSpec& g)
    {
      FieldObsSpecProtobuf pb;
      pb.set_numlocations(g.locations.rows());
      pb.set_locations(matrixString(g.locations));
      pb.set_noiseprob(g.noiseProb);
      return protobufToString(pb);
    }

    void unserialise(const std::string& s, FieldObsSpec& g)
    {
      FieldObsSpecProtobuf pb;
      pb.ParseFromString(s);

      g.locations = stringMatrix(pb.locations(), pb.numlocations());
      g.noiseProb = pb.noiseprob();
    }
    std::string serialise(const FieldObsParams& g)
    {
      FieldObsParamsProtobuf pb;
      pb.set_returnsensordata(g.returnSensorData);
      return protobufToString(pb);
    }
    void unserialise(const std::string& s, FieldObsParams& g)
    {
      FieldObsParamsProtobuf pb;
      pb.ParseFromString(s);
      g.returnSensorData = pb.returnsensordata();
    }
    std::string serialise(const FieldObsResults& g)
    {
      FieldObsResultsProtobuf pb;
      pb.set_likelihood(g.likelihood);
      if (g.readings.size() > 0)
      {
        pb.set_numreadings(g.readings.size());
        pb.set_readings(matrixString(g.readings));
      }
      return protobufToString(pb);
    }
    void unserialise(const std::string& s, FieldObsResults& g)
    {
      FieldObsResultsProtobuf pb;
      pb.ParseFromString(s);
      g.likelihood = pb.likelihood();
      if (pb.has_readings())
      {
        g.readings = stringMatrix(pb.readings(), pb.numreadings());
      }
    }
  }
}
