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
// File: osra_segment.cpp
//
// Defines page segmentation functions
//

#include <iostream> // std::ostream, std::cout
#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX
#include <limits.h> // INT_MAX
#include <algorithm> // std::min(double, double), std::max(double, double)

#include "osra.h"
#include "osra_common.h"
#include "osra_segment.h"


unsigned int distance_between_points(const point_t &p1, const point_t &p2)
{
  return max(abs(p1.x - p2.x), abs(p1.y - p2.y));
}

unsigned int distance_between_segments(const vector<point_t> &s1, const vector<point_t> &s2)
{
  int r = INT_MAX;
  int d;

  for (vector<point_t>::const_iterator i = s1.begin(); i != s1.end(); i++)
    for (vector<point_t>::const_iterator j = s2.begin(); j != s2.end(); j++)
      {
        d = distance_between_points(*i, *j);
        if (d < r)
          r = d;
      }

  /*
  unsigned int ii, jj;
  #pragma omp parallel
  {
  	int priv_min = INT_MAX;
  #pragma omp for
  	for (ii = 0; ii < s1.size(); ii++) {
  		for (jj = 0; jj < s2.size(); jj++) {
  			int d = distance_between_points(s1[ii], s2[jj]);
  			if (d < priv_min)
  				priv_min = d;
  		}
  	}

  	//#pragma omp flush (r)
  	if (priv_min < r) {
  #pragma omp critical
  		{
  			if (priv_min < r)
  				r = priv_min;
  		}
  	}
  }
  */

  return r;
}

void find_connected_components(const Image &image, double threshold, const ColorGray &bgColor,
                               vector<list<point_t> > &segments, vector<vector<point_t> > &margins, bool adaptive)
{
  point_t p;
  list<point_t> points;
  int speckle_area = 2;
  if (adaptive)
    {
      int speckle_side = min(image.columns(), image.rows()) / 200;
      speckle_area = speckle_side * speckle_side;
      if (speckle_area < 2) speckle_area = 2;
    }

  vector<vector<int> > tmp(image.columns(), vector<int> (image.rows(), 0));

  for (unsigned int i = 0; i < image.columns(); i++)
    for (unsigned int j = 0; j < image.rows(); j++)
      if (get_pixel(image, bgColor, i, j, threshold) == 1) // populate with low threshold for future anisotropic smoothing
        tmp[i][j] = 1;


  for (unsigned int i = 0; i < image.columns(); i++)
    for (unsigned int j = 0; j < image.rows(); j++)
      if (tmp[i][j] == 1)
        {
          tmp[i][j] = 2;
          p.x = i;
          p.y = j;
          points.push_back(p);
          list<point_t> new_segment;
          vector<point_t> new_margin;
          int counter = 0;
          point_t p1;
          while (!points.empty())
            {
              p = points.back();
              points.pop_back();
              new_segment.push_back(p);
              tmp[p.x][p.y] = -1;
              bool on_the_margin = false;

              // "k" should be in range "[0 .. image.columns) intercepted with [p.x - 1 .. p.x + 2)" ==> "p.x + 2" should be positive ==> "p.x >= -1"
              // "l" should be in range "[0 .. image.rows)    intercepted with [p.y - 1 .. p.y + 2)" ==> "p.y + 2" should be positive ==> "p.y >= -1"
              if (p.x >= -1 && p.y >= -1)
                {
                  unsigned int x_lower = p.x > 1 ? p.x - 1 : 0; // "k" cannot be less then zero
                  unsigned int y_lower = p.y > 1 ? p.y - 1 : 0; // "l" cannot be less then zero
                  unsigned int x_upper = p.x + 2;
                  if (x_upper > image.columns())
                    x_upper = image.columns();
                  unsigned int y_upper = p.y + 2;
                  if (y_upper > image.rows())
                    y_upper = image.rows();

                  for (int k = x_lower; k < x_upper; k++)
                    for (int l = y_lower; l < y_upper; l++)
                      {
                        if (tmp[k][l] == 1)
                          {
                            p1.x = k;
                            p1.y = l;
                            points.push_back(p1);
                            tmp[k][l] = 2;
                          }
                        else if ((int) k != p.x && (int) l != p.y && tmp[k][l] == 0)
                          {
                            on_the_margin = true;
                          }
                      }
                }

              if (on_the_margin && (new_margin.size() < PARTS_IN_MARGIN || (counter % PARTS_IN_MARGIN) == 0))
                new_margin.push_back(p);
              if (on_the_margin)
                counter++;
            }
          if (segments.size() > MAX_SEGMENTS)
            return;
          if (new_segment.size() > speckle_area)
            {
              segments.push_back(new_segment);
              margins.push_back(new_margin);
            }
        }
}

unsigned int area_ratio(unsigned int a, unsigned int b)
{
  double r = max(a, b) / min(a, b);
  return (unsigned int) r;
}

void build_distance_matrix(const vector<vector<point_t> > &margins, unsigned int max_dist,
                           vector<vector<int> > &distance_matrix, vector<vector<int> > &features, const vector<list<point_t> > &segments,
                           unsigned int max_area_ratio, vector<vector<int> > &area_matrix)
{
  unsigned int d;
  unsigned int ar;

  for (unsigned int s1 = 0; s1 < margins.size(); s1++)
    for (unsigned int s2 = s1 + 1; s2 < margins.size(); s2++)
      if (distance_between_points(margins[s1].front(), margins[s2].front()) < (PARTS_IN_MARGIN
          * margins[s1].size() + PARTS_IN_MARGIN * margins[s2].size()) / 2 + max_dist)
        {
          d = distance_between_segments(margins[s1], margins[s2]);
          if (d < max_dist)
            {
              distance_matrix[s1][s2] = d;
              distance_matrix[s2][s1] = d;
              ar = area_ratio(segments[s1].size(), segments[s2].size());
              //cout << ar << endl;
              area_matrix[s1][s2] = ar;
              area_matrix[s2][s1] = ar;
              if (ar < max_area_ratio && d < max_dist)
                features[ar][d]++;
            }
        }
}

list<list<list<point_t> > > build_explicit_clusters(const list<list<int> > &clusters,
    const vector<list<point_t> > &segments)
{
  list<list<list<point_t> > > explicit_clusters;
  for (list<list<int> >::const_iterator c = clusters.begin(); c != clusters.end(); c++)
    {
      list<list<point_t> > set_of_segments;
      for (list<int>::const_iterator s = c->begin(); s != c->end(); s++)
        if (!segments[*s].empty())
          set_of_segments.push_back(segments[*s]);
      if (!set_of_segments.empty())
        explicit_clusters.push_back(set_of_segments);
    }
  return explicit_clusters;
}

void remove_separators(vector<list<point_t> > &segments, vector<vector<point_t> > &margins, double max_aspect,
                       unsigned int size)
{
  vector<list<point_t> >::iterator s;
  vector<vector<point_t> >::iterator m;
  s = segments.begin();
  m = margins.begin();

  while (s != segments.end() && m != margins.end())
    {
      if (s->size() <= size)
        {
          s++;
          m++;
          continue;
        }

      int stop = INT_MAX, sleft = INT_MAX, sbottom = 0, sright = 0;
      for (list<point_t>::iterator p = s->begin(); p != s->end(); p++)
        {
          if (p->x < sleft)
            sleft = p->x;
          if (p->x > sright)
            sright = p->x;
          if (p->y < stop)
            stop = p->y;
          if (p->y > sbottom)
            sbottom = p->y;
        }
      double aspect = 0;

      if (sright != sleft)
        aspect = 1. * (sbottom - stop+1) / (sright - sleft+1); // where did right and left come from?
      if (aspect > max_aspect || aspect < 1. / max_aspect)
        {
          s = segments.erase(s);
          m = margins.erase(m);
        }
      else
        {
          s++;
          m++;
        }
    }
}

void remove_tables(vector<list<point_t> > &segments, vector<vector<point_t> > &margins, unsigned int size)
{
  vector<list<point_t> >::iterator s;
  vector<vector<point_t> >::iterator m;
  s = segments.begin();
  m = margins.begin();

  while (s != segments.end() && m != margins.end())
    {
      if (m->size() <= size)
        {
          s++;
          m++;
          continue;
        }

      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      int border_count = 0;
      for (vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
        {
          if (p->x < left)
            left = p->x;
          if (p->x > right)
            right = p->x;
          if (p->y < top)
            top = p->y;
          if (p->y > bottom)
            bottom = p->y;
        }

      double aspect = FLT_MAX;
      if (right != left)
        aspect = 1. * (bottom - top) / (right - left);
      if (aspect >= MAX_ASPECT || aspect <= MIN_ASPECT)
        {
          s++;
          m++;
          continue;
        }

      for (vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
        if (p->x - left < 2 || right - p->x < 2 || p->y - top < 2 || bottom - p->y < 2)
          {
            border_count++;
          }

      //cout << border_count << " " << 2 * (right - left) + 2 * (bottom - top) << endl;

      if (border_count > BORDER_COUNT)
        {
          s = segments.erase(s);
          m = margins.erase(m);
        }
      else
        {
          s++;
          m++;
        }
    }
}

list<list<int> > assemble_clusters(const vector<vector<point_t> > &margins, int dist,
                                   const vector<vector<int> > &distance_matrix, vector<int> &avail, bool text,
                                   const vector<vector<int> > &area_matrix)
{
  list<list<int> > clusters;
  list<int> bag;

  for (unsigned int s = 0; s < margins.size(); s++)
    if (avail[s] == 1)
      {
        bag.push_back(s);
        avail[s] = 2;
        list<int> new_cluster;
        while (!bag.empty())
          {
            int c = bag.back();
            bag.pop_back();
            new_cluster.push_back(c);
            avail[c] = 0;
            for (unsigned int i = 0; i < margins.size(); i++)
              if (avail[i] == 1 && distance_matrix[c][i] < dist)
                // && (!text || area_matrix[i][c] <= 10))
                {
                  bag.push_back(i);
                  avail[i] = 2;
                }
          }
        clusters.push_back(new_cluster);
      }

  return (clusters);
}

void remove_text_blocks(const list<list<int> > &clusters, const vector<list<point_t> > &segments, vector<int> &avail)
{
  for (list<list<int> >::const_iterator c = clusters.begin(); c != clusters.end(); c++)
    {
      unsigned int area = 0, square_area = 0;
      double ratio = 0, aspect = 0;
      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      bool fill_below_max = false;

      for (list<int>::const_iterator i = c->begin(); i != c->end(); i++)
        if (!segments[*i].empty())
          {
            int stop = INT_MAX, sleft = INT_MAX, sbottom = 0, sright = 0;
            for (list<point_t>::const_iterator p = segments[*i].begin(); p != segments[*i].end(); p++)
              {
                if (p->x < sleft)
                  sleft = p->x;
                if (p->x > sright)
                  sright = p->x;
                if (p->y < stop)
                  stop = p->y;
                if (p->y > sbottom)
                  sbottom = p->y;
              }

            area = segments[*i].size();
            square_area = (sbottom - stop+1) * (sright - sleft+1);

            if (square_area != 0)
              ratio = 1. * area / square_area;

            if (ratio < MAX_RATIO && ratio > 0)
              fill_below_max = true;

            if (sleft < left)
              left = sleft;
            if (sright > right)
              right = sright;
            if (stop < top)
              top = stop;
            if (sbottom > bottom)
              bottom = sbottom;
          }

      if (c->size() > TEXT_LINE_SIZE)
        {
          if (right != left)
            aspect = 1. * (bottom - top) / (right - left);
          if (aspect < MIN_ASPECT || aspect > MAX_ASPECT || !fill_below_max)
            for (list<int>::const_iterator i = c->begin(); i != c->end(); i++)
              avail[*i] = -1;
        }
    }
}

int locate_first_min(const vector<int> &stats)
{
  int peak = 1;

  for (unsigned int j = 3; j < stats.size(); j++)
    if (stats[j] > stats[j - 1] && stats[j] > stats[j + 1])
      {
        peak = j;
        break;
      }
  int dist = peak;
  for (unsigned int j = peak; j < stats.size(); j++)
    if (stats[j] < stats[j - 1] && stats[j] < stats[j + 1])
      {
        dist = j;
        break;
      }
  return (dist);
}

int locate_max_entropy(const vector<vector<int> > &features, unsigned int max_area_ratio, unsigned int max_dist,
                       vector<int> &stats)
{
  vector<double> entropy(max_area_ratio, 0);

  for (unsigned int i = 1; i < max_area_ratio; i++)
    {
      int count = 0;
      for (unsigned int j = 2; j < max_dist; j++)
        if (features[i][j] == 0)
          count++;
        else
          stats[j]++;

      if (count > 0)
        {
          double probability = 1. * count / (max_dist - 2);
          entropy[i] -= probability * log(probability);
        }
    }
  int start_b = 1;
  for (unsigned int i = 2; i < max_area_ratio; i++)
    {
      if (entropy[i] > entropy[start_b])
        start_b = i;
    }
  return (start_b);
}


void find_arrows_pluses(const vector<vector<point_t> > &margins)
{
  const int len=50;
  for (int i=0; i<margins.size(); i++)
    {
      vector<int> hist(len,0);
      point_t center;
      center.x=0; center.y=0;
      int l=margins[i].size();
      for (int j=0; j<l; j++)
	{
	  center.x += margins[i][j].x;
	  center.y += margins[i][j].y;
	}
      center.x /=l;  // Find the center of mass for the segment margin
      center.y /=l;
      for (int j=0; j<l; j++)
	{
	  int dx = margins[i][j].x-center.x;
	  int dy = margins[i][j].y-center.y;
	  double r=(double)sqrt(dx*dx+dy*dy);
	  double theta=0.;
	  if (dx!=0 || dy!=0)
	    theta = atan2(dy,dx);

	  int bin = (theta+M_PI)*len/(2*M_PI);
	  if (bin>=len) bin -= len;
	  hist[bin]++;          // build a histogram of occurencies in polar coordinates
	}
      int top_pos=0;
      int top_value=0;
      for (int k=0; k<len;k++)
	if (hist[k]>=top_value)
	  {
	    top_pos=k;                       // find the position of the highest peak
	    top_value=hist[k];
	  }
      if (top_value>5)
	{
	  vector<int> peaks(1,top_pos);
	  vector<int> values(1,top_value);
	  for (int k=1; k<len;k++)
	    {
	      int pos=k+top_pos;
	      if (pos>=len) pos -= len;
	      int after=pos+1;
	      int before=pos-1;
	      if (after>=len) after -=len;
	      if (before<0) before +=len;
	      if (hist[before]<hist[pos] && hist[after]<hist[pos] && hist[pos]>=top_value/2)  // find all peaks at least half as high as the top-most
		{
		  peaks.push_back(pos);
		  values.push_back(hist[pos]);
		}
	    }
	  
	  for (int j=0; j<peaks.size(); j++)          // check outside of the peaks is essentially zero
	    for (int k=peaks[j]-2; k<=peaks[j]+2; k++)
	    {
	      int kk=k;
	      if (kk<0) kk += len;
	      if (kk>=len) kk -=len;
	      hist[kk]=0;
	    }
	  bool low=true;
	  for(int k=0; k<len; k++)
	    if (hist[k]>3) low=false;
	  if (low)
	    {
	      if (peaks.size() == 2 && double(values[0])/values[1]>1.3 && fabs(len/2 - fabs(peaks[1]-peaks[0]))<=1)  // only two peaks are present at 180 degrees and at least 1.3 times height difference
		{
		  // we found an arrow!
		}
	    }
	}
    }
  //  exit(0);
}

list<list<list<point_t> > > find_segments(const Image &image, double threshold, const ColorGray &bgColor, bool adaptive, bool is_reaction, bool verbose)
{
  vector<list<point_t> > segments;
  vector<vector<point_t> > margins;
  list<list<list<point_t> > > explicit_clusters;

  // 1m34s

  find_connected_components(image, threshold, bgColor, segments, margins, adaptive);

  if (verbose)
    cout << "Number of segments: " << segments.size() << '.' << endl;

  if (segments.size() > MAX_SEGMENTS)
    {
      segments.clear();
      margins.clear();
    }
  if (is_reaction)
    find_arrows_pluses(margins);

  remove_separators(segments, margins, SEPARATOR_ASPECT, SEPARATOR_AREA);

  remove_tables(segments, margins, SEPARATOR_AREA);

  // 2m22s

  unsigned int max_dist = MAX_DIST;
  unsigned int max_area_ratio = MAX_AREA_RATIO;
  vector<vector<int> > distance_matrix(segments.size(), vector<int> (segments.size(), INT_MAX));
  vector<vector<int> > area_matrix(segments.size(), vector<int> (segments.size(), INT_MAX));
  vector<vector<int> > features(max_area_ratio, vector<int> (max_dist, 0));

  build_distance_matrix(margins, max_dist, distance_matrix, features, segments, max_area_ratio, area_matrix);

  // 2m53s

  vector<int> avail(margins.size(), 1);

  /*
  unsigned int ar;

  for (unsigned int i = 0; i < margins.size(); i++)
  	for (unsigned int j = i + 1; j < margins.size(); j++) {
  		ar = area_ratio(segments[i].size(), segments[j].size());
  		if (ar < max_area_ratio && distance_matrix[i][j] < max_dist)
  			features[ar][distance_matrix[i][j]]++;
  	}
  */

  // 5m53s -> new 4m15s

  vector<int> stats(max_dist, 0);
  int entropy_max = locate_max_entropy(features, max_area_ratio, max_dist, stats);

  int dist = SINGLE_IMAGE_DIST;
  if (entropy_max > THRESHOLD_LEVEL && !adaptive)
    {
      vector<int> text_stats(max_dist, 0);
      for (unsigned int j = 2; j < max_dist; j++)
        {
          text_stats[j] = features[1][j];
          //cout << j << " " << text_stats[j] << endl;
        }

      int dist_text = locate_first_min(text_stats);

      const list<list<int> > &text_blocks = assemble_clusters(margins, dist_text, distance_matrix, avail, true,
                                            area_matrix);
      remove_text_blocks(text_blocks, segments, avail);

      dist = 2 * dist_text;
    }

  for (unsigned int i = 0; i < margins.size(); i++)
    if (avail[i] != -1)
      avail[i] = 1;

  const list<list<int> > &clusters = assemble_clusters(margins, dist, distance_matrix, avail, false, area_matrix);

  explicit_clusters = build_explicit_clusters(clusters, segments);
  return explicit_clusters;
}

/*
void remove_brackets(int left, int right, int top, int bottom, list<list<list<point_t> > >::iterator c) {
	vector < vector<bool> > tmp(right - left + 1, vector<bool> (bottom - top + 1, false));
	vector < vector<bool> > global_pic(right - left + 1, vector<bool> (bottom - top + 1, false));

	for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
		for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
			global_pic[p->x - left][p->y - top] = true;

	bool found = false;
	//Image t(Geometry(right - left + 1, bottom - top + 1), "white");

	for (int i = left + FRAME; i < right - FRAME; i++)
		for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++) {
			vector<point_t> set;
			int x1 = INT_MAX, y1 = INT_MAX, x2 = 0, y2 = 0;
			for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
				if (p->x < i && i + (i - p->x) < right && global_pic[i + (i - p->x) - left][p->y - top]) {
					set.push_back(*p);
					if (p->x < x1)
						x1 = p->x;
					if (p->x > x2)
						x2 = p->x;
					if (p->y < y1)
						y1 = p->y;
					if (p->y > y2)
						y2 = p->y;
				}

			if (set.size() > 100 && (i - x2) > 5 && (x2 - x1) > 3 && (y2 - y1) > 5 && (x2 - x1) < (y2 - y1)) {
				int x = x2 - x1 + 1;
				int y = y2 - y1 + 1;
				int f = 1;
				if (y > 40)
					f = y / 40;
				x /= f;
				y /= f;

				unsigned char *pic = (unsigned char *) malloc(x * y);
				for (int j = 0; j < x * y; j++)
					pic[j] = 255;
				for (unsigned int p = 0; p < set.size(); p++)
					if ((set[p].y - y1) / f < y && (set[p].x - x1) / f < x)
						pic[((set[p].y - y1) / f) * x + (set[p].x - x1) / f] = 0;
				bool res = detect_bracket(x, y, pic);
				if (res) {
					//cout << set.size() << endl;
					for (unsigned int p = 0; p < set.size(); p++) {
						tmp[set[p].x - left][set[p].y - top] = true;
						tmp[i + (i - set[p].x) - left][set[p].y - top] = true;
						//t.pixelColor(set[p].x - left, set[p].y - top, "black");
						//t.pixelColor(i + (i - set[p].x) - left, set[p].y - top, "black");
					}
					found = true;
				}
			}
		}

	if (found) {
		//t.write("t.png");
		//exit(0);
		list<list<point_t> >::iterator s1 = c->begin();
		while (s1 != c->end()) {
			list<point_t>::iterator p1 = s1->begin();
			while (p1 != s1->end())
				if (tmp[p1->x - left][p1->y - top])
					p1 = s1->erase(p1);
				else
					p1++;
			if (s1->size() > 0)
				s1++;
			else
				s1 = c->erase(s1);
		}
	}
}
*/

int prune_clusters(list<list<list<point_t> > > &clusters, vector<box_t> &boxes)
{
  int n_boxes = 0;
  list<list<list<point_t> > >::iterator c = clusters.begin();

  while (c != clusters.end())
    {
      unsigned int area = 0, square_area = 0;
      double ratio = 0, aspect = 0;
      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      bool fill_below_max = false;
      for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
        {
          int stop = INT_MAX, sleft = INT_MAX, sbottom = 0, sright = 0;
          for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
            {
              if (p->x < sleft)
                sleft = p->x;
              if (p->x > sright)
                sright = p->x;
              if (p->y < stop)
                stop = p->y;
              if (p->y > sbottom)
                sbottom = p->y;
            }

          area = s->size();
          square_area = (sbottom - stop+1) * (sright - sleft+1);
          ratio = 0;
          if (square_area != 0)
            ratio = 1. * area / square_area;

          if (ratio < MAX_RATIO && ratio > 0)
            fill_below_max = true;

          if (sleft < left)
            left = sleft;
          if (sright > right)
            right = sright;
          if (stop < top)
            top = stop;
          if (sbottom > bottom)
            bottom = sbottom;
        }

      if (right != left)
        aspect = 1. * (bottom - top) / (right - left);

      if (fill_below_max && aspect > MIN_ASPECT && aspect < MAX_ASPECT)
        {
          box_t b1;
          boxes.push_back(b1);
          boxes[n_boxes].x1 = left;
          boxes[n_boxes].y1 = top;
          boxes[n_boxes].x2 = right;
          boxes[n_boxes].y2 = bottom;

          //remove_brackets(left, right, top, bottom, c);

          for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
            for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
              boxes[n_boxes].c.push_back(*p);
          c++;
          n_boxes++;
        }
      else
        {
          c = clusters.erase(c);
        }
    }
  return (n_boxes);
}

