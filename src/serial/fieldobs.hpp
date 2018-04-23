//!
//! Field Observation Forward Model serialisation.
//!
//! \file serial/fieldobs.hpp
//! \author Richard Scalzo
//! \date March, 2018
//! \license Affero General Public License version 3 or later
//! \copyright (c) 2014, 2018 NICTA + USyd
//!

#pragma once

#include "datatype/datatypes.hpp"

#include <string>

namespace obsidian
{
  namespace comms
  {
    std::string serialise(const FieldObsSpec& g);
    std::string serialise(const FieldObsParams& g);
    std::string serialise(const FieldObsResults& g);

    void unserialise(const std::string& s, FieldObsSpec& g);
    void unserialise(const std::string& s, FieldObsParams& g);
    void unserialise(const std::string& s, FieldObsResults& g);
  }
}
