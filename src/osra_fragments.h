/******************************************************************************
 OSRA: Optical Structure Recognition Application

 This is a U.S. Government work (2007-2011) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2011 (igorf@helix.nih.gov)

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

// Header: osra_fragments.h
//
// Declares operations on molecular fragments
//

#include "osra_openbabel.h"

//struct: fragment_s
// used by <osra.cpp::populate_fragments()> to split chemical structure into unconnected molecules.
struct fragment_s
{
  //int: x1,y1,x2,y2
  //top left and bottom right coordinates of the fragment
  int x1, y1, x2, y2;
  //array: atom
  //vector of atom indices for atoms in a molecule of this fragment
  vector<int> atom;
};
//typedef: fragment_t
//defines fragment_t type based on fragment_s struct
typedef struct fragment_s fragment_t;

vector<vector<int> > find_fragments(const vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom);
int reconnect_fragments(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double avg);
vector<fragment_t> populate_fragments(const vector<vector<int> > &frags, const vector<atom_t> &atom);
bool comp_fragments(const fragment_t &aa, const fragment_t &bb);
