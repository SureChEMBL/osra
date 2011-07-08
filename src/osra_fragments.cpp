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

// File osra_fragments.cpp
//
// Defines operations on molecular fragments
//

#include <float.h> // FLT_MAX
#include <limits.h> // INT_MAX
#include "osra.h"
#include "osra_fragments.h"

double atom_distance(const vector<atom_t> &atom, int a, int b)
{
  return (distance(atom[a].x, atom[a].y, atom[b].x, atom[b].y));
}

/**
 * TODO: Returning the vector from the stack causes copy constructor to trigger, which is inefficient.
 * Consider passing the vector as a reference.
 */
vector<vector<int> > find_fragments(const vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom)
{
  vector<vector<int> > frags;
  vector<int> pool;
  int n = 0;

  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && atom[bond[i].a].exists && atom[bond[i].b].exists)
      pool.push_back(i);

  while (!pool.empty())
    {
      frags.resize(n + 1);
      frags[n].push_back(bond[pool.back()].a);
      frags[n].push_back(bond[pool.back()].b);
      pool.pop_back();
      bool found = true;

      while (found)
        {
          found = false;
          unsigned int i = 0;
          while (i < pool.size())
            {
              bool found_a = false;
              bool found_b = false;
              bool newfound = false;
              for (unsigned int k = 0; k < frags[n].size(); k++)
                {
                  if (frags[n][k] == bond[pool[i]].a)
                    found_a = true;
                  else if (frags[n][k] == bond[pool[i]].b)
                    found_b = true;
                }
              if (found_a && !found_b)
                {
                  frags[n].push_back(bond[pool[i]].b);
                  pool.erase(pool.begin() + i);
                  found = true;
                  newfound = true;
                }
              if (!found_a && found_b)
                {
                  frags[n].push_back(bond[pool[i]].a);
                  pool.erase(pool.begin() + i);
                  found = true;
                  newfound = true;
                }
              if (found_a && found_b)
                {
                  pool.erase(pool.begin() + i);
                  newfound = true;
                }
              if (!newfound)
                i++;
            }
        }
      n++;
    }
  return (frags);
}

int reconnect_fragments(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double avg)
{
  vector<vector<int> > frags = find_fragments(bond, n_bond, atom);

  if (frags.size() <= 3)
    {
      for (unsigned int i = 0; i < frags.size(); i++)
        if (frags[i].size() > 2)
          for (unsigned int j = i + 1; j < frags.size(); j++)
            if (frags[j].size() > 2)
              {
                double l = FLT_MAX;
                int atom1 = 0, atom2 = 0;
                for (unsigned int ii = 0; ii < frags[i].size(); ii++)
                  for (unsigned int jj = 0; jj < frags[j].size(); jj++)
                    {
                      double d = atom_distance(atom, frags[i][ii], frags[j][jj]);
                      if (d < l)
                        {
                          l = d;
                          atom1 = frags[i][ii];
                          atom2 = frags[j][jj];
                        }
                    }
                if (l < 1.1 * avg && l > avg / 3)
                  {
                    bond_t b1;
                    bond.push_back(b1);
                    bond[n_bond].a = atom1;
                    bond[n_bond].exists = true;
                    bond[n_bond].type = 1;
                    bond[n_bond].b = atom2;
                    bond[n_bond].curve = atom[atom1].curve;
                    bond[n_bond].hash = false;
                    bond[n_bond].wedge = false;
                    bond[n_bond].up = false;
                    bond[n_bond].down = false;
                    bond[n_bond].Small = false;
                    bond[n_bond].arom = false;
                    bond[n_bond].conjoined = false;
                    n_bond++;
                  }
                if (l < avg / 3)
                  {
                    atom[atom2].x = atom[atom1].x;
                    atom[atom2].y = atom[atom1].y;
                  }
              }
    }

  return (n_bond);
}


/**
 * TODO: Returning the vector from the stack causes copy constructor to trigger, which is inefficient.
 * Consider passing the vector as a reference.
 */
vector<fragment_t> populate_fragments(const vector<vector<int> > &frags, const vector<atom_t> &atom)
{
  vector<fragment_t> r;

  for (unsigned int i = 0; i < frags.size(); i++)
    {
      fragment_t f;
      f.x1 = INT_MAX;
      f.x2 = 0;
      f.y1 = INT_MAX;
      f.y2 = 0;

      for (unsigned j = 0; j < frags[i].size(); j++)
        {
          f.atom.push_back(frags[i][j]);
          if ((int) atom[frags[i][j]].x < f.x1)
            f.x1 = (int) atom[frags[i][j]].x;
          if ((int) atom[frags[i][j]].x > f.x2)
            f.x2 = (int) atom[frags[i][j]].x;
          if ((int) atom[frags[i][j]].y < f.y1)
            f.y1 = (int) atom[frags[i][j]].y;
          if ((int) atom[frags[i][j]].y > f.y2)
            f.y2 = (int) atom[frags[i][j]].y;
        }
      r.push_back(f);
    }
  return (r);
}

bool comp_fragments(const fragment_t &aa, const fragment_t &bb)
{
  if (aa.y2 < bb.y1)
    return (true);
  if (aa.y1 > bb.y2)
    return (false);
  if (aa.x1 > bb.x1)
    return (false);
  if (aa.x1 < bb.x1)
    return (true);

  return (false);
}
