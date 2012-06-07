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

// Header: osra_reaction.h
//
// Defines functions dealing with generating a reaction type output
//

#define SUBSTITUTE_REACTION_FORMAT "mol"


//
// Section: Functions
//

// Function: arrange_reactions
//
// Create a reaction representation for input vector of structures
//
// Parameters:
//      arrows - a vector of arrow_t objects representing arrows found during segmentation
//      page_of_boxes - a vector of box_t objects representing bounding boxes of molecules
//      pluses - a vector of plus sing centers
//      results - a vector of strings to represent output results
//      page_of_structures - input vector of reactants, intermediates and products
//      output_format - format of the returned result, i.e. rsmi or cmlr
//


void arrange_reactions(vector<arrow_t> &arrows, const vector<box_t> &page_of_boxes, const vector<point_t> &pluses, vector<string> &results,
		       const vector<string> &page_of_structures, const string &output_format);
