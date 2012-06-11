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

#include "osra_segment.h"
#include "osra_reaction.h"
#include "osra_common.h"

using namespace OpenBabel;
using namespace std;

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
string convert_page_to_reaction(const vector<string> &page_of_structures, const string &output_format, const vector <int> &reactants, const vector <int> &products, string value)
{
  string reaction;
  OBConversion conv;
  conv.SetInAndOutFormats(SUBSTITUTE_REACTION_FORMAT,output_format.c_str());
  ostringstream strstr;
  
  OBReaction react;
  for (int j=0; j<reactants.size(); j++)
    {
      shared_ptr<OBMol> reactant(new OBMol);
      conv.ReadString(reactant.get(), page_of_structures[reactants[j]]);
      react.AddReactant(reactant);
    }
  for (int j=0; j<products.size(); j++)
    {
      shared_ptr<OBMol> product(new OBMol);
      conv.ReadString(product.get(), page_of_structures[products[j]]);
      react.AddProduct(product);
    }
  //	  react.AddAgent(transition);

  trim(value);
  if (!value.empty())
    {
      OBPairData *label = new OBPairData;
      label->SetAttribute("OSRA_REACTION_AGENT");
      label->SetValue(value.c_str());
      react.SetData(label);
    }
  strstr << conv.WriteString(&react, true);
  if (output_format == "rsmi" && !strstr.str().empty() && !value.empty())
    strstr << " " << value;
  reaction = strstr.str();

  return(reaction);
}

string convert_to_smiles_agent_structure(string structure)
{
  OBConversion conv;
  conv.SetInAndOutFormats(SUBSTITUTE_REACTION_FORMAT,"smi");
  string result;
  OBMol mol;
  if (conv.ReadString(&mol,structure))
    result = conv.WriteString(&mol,true);
  return result;
}

void linear_arrow_sort(vector<arrow_t> &arrows)
{
  point_t start;
  start.x=0;
  start.y=0;
  vector<arrow_t> new_arrows;
  while (!arrows.empty())
    {
      double d=FLT_MAX;
      int closest=0;
      for (int i=0; i<arrows.size(); i++)
	if (distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y)<d)
	  {
	    d = distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y);
	    closest = i;
	  }
      new_arrows.push_back(arrows[closest]);
      start=arrows[closest].head;
      arrows.erase(arrows.begin()+closest);
    }
  arrows = new_arrows;
}

double distance_from_box(const point_t &p, const box_t &b)
{
  return(distance(p.x,p.y,(b.x1+b.x2)/2,(b.y1+b.y2)/2));
}

bool comp_structures(const pair<int,box_t> &a, const pair<int,box_t> &b)
{
  if (a.second.y2 < b.second.y1)
    return (true);
  if (a.second.x2 < b.second.x1)
    return (true);
  return (false);
}

void arrange_reactions(vector<arrow_t> &arrows, const vector<box_t> &page_of_boxes, const vector<point_t> &pluses, vector<string> &results,
		       const vector<string> &page_of_structures,  const string &output_format)
{
  vector < vector<pair<int,box_t> > > before(arrows.size()+1);
  // arrange arrows in head to tail fashion
  linear_arrow_sort(arrows);
  //for (int i=0; i<arrows.size(); i++)
  //  cout<<arrows[i].tail.x<<","<<arrows[i].tail.y<<" "<<arrows[i].head.x<<","<<arrows[i].head.y<<endl;

  // arrange structures to best fit between arrows
  for (int i=0; i<page_of_boxes.size(); i++)
    {
      double rt = FLT_MAX;
      int j_tail=0;
      double rh = FLT_MAX;
      int j_head=0;
      bool agent_structure = false;
      for (int j=0; j<arrows.size(); j++)
	{
	  double r = distance_from_box(arrows[j].tail, page_of_boxes[i]);
	  double ry = distance_from_bond_y(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,(page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2);
	  double l = distance(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y);
	  double rx = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,(page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2);
	  double cr = distance((page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2,(arrows[j].tail.x+arrows[j].head.x)/2,(arrows[j].tail.y+arrows[j].head.y)/2);
	  if (rx>0 && rx<l && cr<max(page_of_boxes[i].x2-page_of_boxes[i].x1, page_of_boxes[i].y2-page_of_boxes[i].y1))
	    {
	      agent_structure = true;
	      arrows[j].agent += convert_to_smiles_agent_structure(page_of_structures[i]);
	      break;
	    }
	  if (fabs(ry)<min(page_of_boxes[i].x2-page_of_boxes[i].x1, page_of_boxes[i].y2-page_of_boxes[i].y1))
	    {
	      if (r<rt)
		{
		  rt = r;
		  j_tail = j;
		}
	      r = distance_from_box(arrows[j].head, page_of_boxes[i]);
	      if (r<rh)
		{
		  rh = r;
		  j_head = j;
		}
	    }
	}
      if (rh<FLT_MAX || rt<FLT_MAX && !agent_structure)
	{
	  if (rt<rh)
	    before[j_tail].push_back(make_pair(i,page_of_boxes[i]));
	  else
	    before[j_head+1].push_back(make_pair(i,page_of_boxes[i]));
	}

    }


  // products can be on the next line
  for (int i=1; i<before.size(); i++)
    {
      if (before[i].empty())
	for (int j=0; j<before[i-1].size(); j++)
	  if (before[i-1][j].second.y1 > arrows[i-1].head.y)
	  {
	    before[i].push_back(before[i-1][j]);
	    before[i-1][j].first=-1;
	  }
    }


  for (int i=0; i<before.size(); i++)
    sort(before[i].begin(),before[i].end(),comp_structures);


  vector < vector <bool> > is_plus(page_of_boxes.size(), vector <bool> (page_of_boxes.size(), false));
  // arrange plus signs between boxes
  for (int i=0; i<before.size(); i++)
    for (int j=1; j<before[i].size(); j++)
      {
	int k = before[i][j].first;
	int l = before[i][j-1].first;
	if (k>=0 && l>=0)
	  {
	    box_t b=before[i][j].second;
	    box_t a=before[i][j-1].second;
	    for (int m=0; m<pluses.size(); m++)
	      {
		double d=distance_from_bond_y((a.x1+a.x2)/2,(a.y1+a.y2)/2,(b.x1+b.x2)/2,(b.y1+b.y2)/2,pluses[m].x, pluses[m].y);
		if (fabs(d)<min(a.y2-a.y1,b.y2-b.y1)/2 && pluses[m].x>a.x2 && pluses[m].x<b.x1)
		  {
		    is_plus[k][l] = true;
		    is_plus[l][k] = true;
		  }
		// after plus things can be on the next line
		d = pluses[m].y - (a.y2 + a.y1)/2;
		if (pluses[m].x>a.x2 && fabs(d)<(a.y2-a.y1)/2 && b.y1>a.y2)
		  {
		    is_plus[k][l] = true;
		    is_plus[l][k] = true;
		  }
		
	      }
	  }
      }

  // extract reactions, if any
  for (int i=0; i<arrows.size(); i++)
    {
      vector <int> r,p;
      int ii = before[i].size()-1;
      while (ii>=0 && before[i][ii].first<0) ii--;

      if (ii>=0 && before[i][ii].first>=0)
	r.push_back(before[i][ii].first);

      for (int j=ii-1; j>=0; j--)
	{
	  int k = before[i][j].first;
	  int l = before[i][j+1].first;
	  if (k>=0 && l>=0)
	    {
	      if (is_plus[k][l])
		r.push_back(k);
	      else
		break;

	    }
	}

      if (!before[i+1].empty())
	{
	  ii = 0;
	  while (ii<before[i+1].size() && before[i+1][0].first<0) ii++;

	  if (ii<before[i+1].size() && before[i+1][ii].first>=0)
	    p.push_back(before[i+1][ii].first);
	  for (int j=ii+1; j<before[i+1].size(); j++)
	    {
	      int k = before[i+1][j].first;
	      int l = before[i+1][j-1].first;
	      if (k>=0 && l>=0)
		{
		  if (is_plus[k][l])
		    p.push_back(k);
		  else
		    break;
		}
	    }
	}

      if (!r.empty() && !p.empty())
	{
	  string result=convert_page_to_reaction(page_of_structures,output_format, r, p, arrows[i].agent);
	  trim(result);
	  if (!result.empty())
	    results.push_back(result);
	}
    }
}
