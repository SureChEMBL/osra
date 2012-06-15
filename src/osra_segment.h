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
// Header: osra_segment.h
//
// Declares page segmentation functions
//

#ifndef OSRA_SEGMENT_H
#define OSRA_SEGMENT_H

#include <list> // sdt::list
#include <vector> // std::vector
#include <Magick++.h>

using namespace std;
using namespace Magick;

// struct: point_s
//      a point of the image, used by image segmentation routines
struct point_s
{
  // int: x,y
  //    coordinates of the image point
  int x, y;
};
// typedef: point_t
//      defines point_t type based on point_s struct
typedef struct point_s point_t;

// struct: box_s
//      encompassing box structure for image segmentation
struct box_s
{
  // int: x1, y1, x2, y2
  //    coordinates of top-left and bottom-right corners
  int x1, y1, x2, y2;
  // array: c
  //    vector of points in the box
  vector<point_t> c;
};
// typedef: box_t
//      defines box_t type based on box_s struct
typedef struct box_s box_t;

// struct: arrow_s
// coordinates of tail and head of an arrow
struct arrow_s
{
  // point_t: tail, head
  // tail and head of an arrow as points
  point_t tail,head;
  string agent;
  bool linebreak;
  bool reversible;
  bool remove;
};
// typedef: arrow_t
// defines arrow_t type based on arrow_s struct
typedef struct arrow_s arrow_t;


//
// Section: Functions
//

// Function: find_segments()
//
// Performs page segmentation to different regions (text/graphics/linear etc.)
//
// Parameters:
// image - page image
// threshold - black-white binarization threshold
// bgColor - background color
// adaptive - flag set if adaptive thresholding has been used in grayscale conversion
// is_reaction - flag set if we're looking for reaction-specific symbols (arrows, plus signs etc.)
// arrows - a vector of arrows found during segmentation
// pluses - a vector of plus centers found during segmentation
// verbose - flag set for verbose reporting
//
// Returns:
// A list of clusters, each of which is a list of  connected segments each of which is a list of points
list<list<list<point_t> > > find_segments(const Image &image, double threshold, const ColorGray &bgColor, bool adaptive, bool is_reaction, vector<arrow_t> &arrows, vector<point_t> &pluses, bool verbose);

// Function: prune_clusters()
//
// Prunes the list of clusters and retains only molecular structure images
//
// Parameters:
// clusters - a list of clusters detected by <find_segments()>
// boxes - a vector of <box_t> objects for molecular structure images
//
// Returns:
// Number of molecular structure images
int prune_clusters(list<list<list<point_t> > > &clusters, vector<box_t> &boxes);
#endif
