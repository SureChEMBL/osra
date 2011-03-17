/******************************************************************************
 OSRA: Optical Structure Recognition

 This is a U.S. Government work (2007-2010) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2010 (igorf@helix.nih.gov)

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

// Header: osra.h
//
// Defines types and functions of OSRA library.
//

#include <string> // std::string
#include <ostream> // std:ostream

using namespace std;

//
// Section: Functions
//

// Function: osra_process_image()
//
// Parameters:
//      image_data - the binary image
//
// Returns:
//      0, if processing was completed successfully
int osra_process_image(
#ifdef OSRA_LIB
  const char *image_data,
  int image_length,
  ostream &structure_output_stream,
#else
  const string &input_file,
  const string &output_file,
#endif
  int rotate = 0,
  bool invert = false,
  int input_resolution = 0,
  double threshold = 0,
  int do_unpaper = 0,
  bool jaggy = false,
  bool adaptive = false,
  const string &output_format = "smi",
  const string &embedded_format = "",
  bool show_confidence = false,
  bool show_resolution_guess = false,
  bool show_page = false,
  bool show_coordinates = false,
  bool show_avg_bond_length = false,
  const string &osra_dir = "",
  const string &spelling_file = "",
  const string &superatom_file = "",
  bool debug = false,
  bool verbose = false,
  const string &output_image_file_prefix = "",
  const string &resize = ""
);
