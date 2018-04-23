/**
 * Contains the basic GDF datatypes for communication between forward models.
 *
 * @file testfieldobs.hpp
 * @author Richard Scalzo
 * @date 2018
 * @license General Public License version 3 or later.
 * @copyright (c) 2014, 2018 NICTA + USyd
 */

#include "test/fieldobs.hpp"
#include "test/input.hpp"
#include "input.hpp"

namespace obsidian
{

TEST_F(InputTest, testSpec)
{
  generateVariations<FieldObsSpec>(testSpecCSV<ForwardModel::FIELDOBS>);
}

TEST_F(InputTest, testResult)
{
  generateVariations<FieldObsResults>(testResultsCSV<ForwardModel::FIELDOBS>);
}

}
