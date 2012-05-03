/******************************************************************************
 OSRA: Optical Structure Recognition

 This is a U.S. Government work (2007-2012) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2012 (igorf@helix.nih.gov)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

#include <sstream> // std:ostringstream
#include <string> // std:string
#include <vector> // std::vector


#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/reaction.h>

using namespace OpenBabel;
using namespace std;


// Function: convert_page_to_reaction(
//
// Create a reaction representation for input vector of structures
//
// Parameters:
//      page_of_structures - input vector of reactants, intermediates and products
//      output_format - format of the returned result, i.e. rsmi or cmlr
//
// Returns:
//      resulting reaction in the format set up by output_format parameter
//
string convert_page_to_reaction(const vector<string> &page_of_structures, const string &output_format)
{
  string reaction;
  OBConversion conv;
  conv.SetInAndOutFormats("mol",output_format.c_str());
  ostringstream strstr;

  if (page_of_structures.size() > 1)
    {
      OBReaction *react =  new OBReaction; // OBReaction doesn't seem to have a desctructor!!! 
      OBMol reactant;
      conv.ReadString(&reactant, page_of_structures[0]);
      react->AddReactant(shared_ptr<OBMol>(&reactant));
      OBMol product;
      conv.ReadString(&product, page_of_structures[1]);
      react->AddProduct(shared_ptr<OBMol>(&product));
      strstr << conv.WriteString(react, true);
      reaction = strstr.str();
      reactant.Clear();
      product.Clear();
      //delete react;
    }

  return(reaction);
}
