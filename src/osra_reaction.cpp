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
#include <set>

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
string convert_page_to_reaction(const vector<string> &page_of_structures, const string &output_format, const vector <int> &reactants, const vector <int> &products, string value, bool reversible)
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
  if (reversible) react.SetReversible(true);

  trim(value);
  if (!value.empty())
    {
      //      OBPairData *label = new OBPairData;
      //      label->SetAttribute("OSRA_REACTION_AGENT");
      //      label->SetValue(value.c_str());
      //      react.SetData(label);
      react.SetComment(value);
    }
  strstr << conv.WriteString(&react, true);
  if (output_format == "rsmi" && !strstr.str().empty() && !value.empty())
    strstr << " " << value;
  reaction = strstr.str();

  return(reaction);
}

string convert_to_smiles_agent_structure(const string &structure)
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
      bool found = false;
      for (int i=0; i<arrows.size(); i++)
	{
	  bool linebreak = false;
	  if (!new_arrows.empty())
	    linebreak = arrows[i].head.x>arrows[i].tail.x && new_arrows.back().head.x>new_arrows.back().tail.x 
	      && abs(arrows[i].head.y-arrows[i].tail.y)<5  && abs(new_arrows.back().head.y-new_arrows.back().tail.y)<5
	      && min(arrows[i].head.y,arrows[i].tail.y)-max(new_arrows.back().head.y,new_arrows.back().tail.y)>MAX_FONT_HEIGHT;
	      //&& arrows[i].tail.x<new_arrows.back().head.x;
	  
	  if (distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y)<d && (!linebreak || new_arrows.back().linebreak))
	    {
	      d = distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y);
	      closest = i;
	      found = true;
	    }
	}
     
      if (found)
	{
	  new_arrows.push_back(arrows[closest]);
	  start=arrows[closest].head;
	  arrows.erase(arrows.begin()+closest);
	}
      else
	{
	  start.x = 0;
	  start.y = 0;
	  if (!new_arrows.empty())
	    {
	      new_arrows.back().linebreak = true;
	      start.y = max(new_arrows.back().tail.y,new_arrows.back().head.y);
	    }
	}
    }
  arrows = new_arrows;
}

void check_the_last_arrow_linebreak(vector<arrow_t> &arrows,const vector<box_t> &page_of_boxes)
{
  if (arrows.back().head.x>arrows.back().tail.x && abs(arrows.back().head.y-arrows.back().tail.y)<5)
    {
      bool linebreak = true;
      for (int i=0; i<page_of_boxes.size(); i++)
	if (page_of_boxes[i].x1 > arrows.back().head.x && page_of_boxes[i].x1 - arrows.back().head.x < MAX_DISTANCE_BETWEEN_ARROWS
	    && page_of_boxes[i].y1 < arrows.back().head.y && page_of_boxes[i].y2 > arrows.back().head.y)
	  linebreak = false;
      arrows.back().linebreak = linebreak;
    }
}

void mark_reversible(vector<arrow_t> &arrows)
{
  if (arrows.size()<2) return;
  for (int i=0; i<arrows.size(); i++)
    for (int j=i+1; j<arrows.size(); j++)
      if (angle4(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,
		 arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y)<-D_T_TOLERANCE)
      {
	double d1=distance_from_bond_y(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].tail.x,arrows[j].tail.y);
	double d2=distance_from_bond_y(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].head.x,arrows[j].head.y);
	if (max(fabs(d1),fabs(d2))<2*MAX_FONT_HEIGHT)
	  {
	    
	    double l = distance(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y);
	    double l1 = distance_from_bond_x_a(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].tail.x,arrows[j].tail.y);
	    double l2 = distance_from_bond_x_a(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].head.x,arrows[j].head.y);
	    if (fabs(l-l1)<5 && fabs(l2)<5)
	      {
		double x1=0,y1=0, x2=1,y2=1;
		if (i>0)
		  {
		    x1=arrows[i-1].tail.x;
		    y1=arrows[i-1].tail.y;
		    int k = i+1;
		    while (k<arrows.size() && k==j) k++;
		    if (k==arrows.size())
		      {
			x2=arrows[i-1].head.x;
			y2=arrows[i-1].head.y;
		      }
		    else 
		      {
			x1=arrows[i-1].head.x;
			y1=arrows[i-1].head.y;
			x2=arrows[k].tail.x;
			y2=arrows[k].tail.y;
		      }
		  }
		double a1 = angle4(x1,y1,x2,y2,arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y);
		double a2 = angle4(x1,y1,x2,y2,arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y);
		if (a1>a2)
		  {
		    arrows[i].reversible = true;
		    arrows[j].remove=true;
		  }
		else
		  {
		    arrows[j].reversible = true;
		    arrows[i].remove=true;
		  }
	      }
	  }
      }
      
  vector<arrow_t>::iterator i=arrows.begin();
  while (i!=arrows.end())
    {
      if (i->remove) i=arrows.erase(i);
      else i++;
    }
}

double distance_from_box(const point_t &p, const box_t &b)
{
  if (p.x < b.x1 && p.y < b.y1)  return(distance(p.x,p.y,b.x1,b.y1));
  if (p.x > b.x1 && p.x < b.x2 && p.y < b.y1)  return(b.y1-p.y);
  if (p.x > b.x2 && p.y < b.y1)  return(distance(p.x,p.y,b.x2,b.y1));
  if (p.x < b.x1 && p.y > b.y1 && p.y < b.y2)  return(b.x1-p.x);
  if (p.x > b.x2 && p.y > b.y1 && p.y < b.y2)  return(p.x-b.x2);
if (p.x < b.x1 && p.y > b.y2)  return(distance(p.x,p.y,b.x1,b.y2));
if (p.x > b.x1 && p.x < b.x2 && p.y > b.y2)  return(p.y-b.y2);
if (p.x > b.x2 && p.y > b.y2)  return(distance(p.x,p.y,b.x2,b.y2));
return 0;
}

double distance_between_boxes(const box_t &a, const box_t &b)
{
  point_t p;
  p.x = a.x1;
  p.y = a.y1;
  double dab = distance_from_box(p,b);
  p.x = a.x1;
  p.y = a.y2;
 dab = min(dab,distance_from_box(p,b));
  p.x = a.x2;
  p.y = a.y1;
  dab = min(dab,distance_from_box(p,b));
  p.x = a.x2;
  p.y = a.y2;
  dab = min(dab,distance_from_box(p,b));

  p.x = b.x1;
  p.y = b.y1;
  dab = min(dab,distance_from_box(p,a));
  p.x = b.x2;
  p.y = b.y1;
  dab = min(dab,distance_from_box(p,a));
  p.x = b.x1;
  p.y = b.y2;
  dab = min(dab,distance_from_box(p,a));
  p.x = b.x2;
  p.y = b.y2;
  dab = min(dab,distance_from_box(p,a));
  return(dab);
}

vector<int>  sort_boxes_one_by_one(point_t p,  vector<int> b, const vector<box_t> &page_of_boxes, vector<int> *a=NULL)
{
  vector<int> t;
  while (!b.empty())
    {
      double d = FLT_MAX;
      int min_j = 0;
      if (t.empty() && a!=NULL && !a->empty())
	{
	  int ii = a->back();
	  for (int j=0; j<b.size(); j++)
	    {
	      int jj = b[j];
	      if (d > distance_between_boxes(page_of_boxes[ii],page_of_boxes[jj]))
		{
		  d = distance_between_boxes(page_of_boxes[ii],page_of_boxes[jj]);
		  min_j = j;
		}
	    }
	  a=NULL;
	}
      else
	for (int j=0; j<b.size(); j++)
	  if (d>distance_from_box(p,page_of_boxes[b[j]]))
	    {
	      d=distance_from_box(p,page_of_boxes[b[j]]);
	      min_j = j;
	    }
      int jj = b[min_j];
      if (!t.empty())
	{
	  int ii = t.back();
	  d = distance_between_boxes(page_of_boxes[ii],page_of_boxes[jj]);
	}
      if (d<MAX_DISTANCE_BETWEEN_ARROWS)
	{
	  t.push_back(jj);
	  p.x = (page_of_boxes[jj].x1+page_of_boxes[jj].x2)/2;
	  p.y = (page_of_boxes[jj].y1+page_of_boxes[jj].y2)/2;
	}
      b.erase(b.begin()+min_j);
    }
  return t;
}

void sort_boxes_from_arrows(const vector<arrow_t> &arrows,  vector < vector<int> > &before, const vector<box_t> &page_of_boxes)
{
  if (before.empty() || arrows.empty() || page_of_boxes.empty()) return;
  vector<int> t=sort_boxes_one_by_one(arrows[0].tail,before[0],page_of_boxes);
  reverse(t.begin(),t.end());
  before[0]=t;

  for (int i=1; i<before.size(); i++)
    {
      vector<int> t1;
      if (!arrows[i-1].linebreak)
	t1 = sort_boxes_one_by_one(arrows[i-1].head,before[i],page_of_boxes);
      else
	{
	  t1 = sort_boxes_one_by_one(arrows[i-1].tail,before[i],page_of_boxes,&before[i-1]);
	  vector<int> t3;
	  while (!t1.empty())
	    {
	      double d=FLT_MAX;
	      int min_k=0;
	      for (int k=0; k<t1.size(); k++)
		{
		  int kk=t1[k];
		  if (d>page_of_boxes[kk].x1)
		    {
		      d = page_of_boxes[kk].x1;
		      min_k=k;
		    }
		}
	      t3.push_back(t1[min_k]);
	      t1.erase(t1.begin()+min_k);
	    }
	  t1 = t3;
	}

      vector<int> t2=sort_boxes_one_by_one(arrows[i].tail,before[i],page_of_boxes);
      reverse(t2.begin(),t2.end());
      set<int> d(t1.begin(),t1.end());
      before[i] = t1;
      for (int j=0; j<t2.size(); j++)
	if (d.find(t2[j])==d.end())
	  {
	    before[i].push_back(t2[j]);
	    d.insert(t2[j]);
	  }
    }
}

void arrange_structures_between_arrows(vector<arrow_t> &arrows,  vector < vector<int> > &before, const vector<box_t> &page_of_boxes, const vector<string> &page_of_structures)
{
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
	  double rx1 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x1,page_of_boxes[i].y1);
	  double rx2 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x1,page_of_boxes[i].y2);
	  double rx3 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x2,page_of_boxes[i].y1);
	  double rx4 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x2,page_of_boxes[i].y2);
	  double cr = distance((page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2,(arrows[j].tail.x+arrows[j].head.x)/2,(arrows[j].tail.y+arrows[j].head.y)/2);
	  if (rx1>=0 && rx1<=l && rx2>=0 && rx2<=l && rx3>=0 && rx3<=l && rx4>=0 && rx4<=l && 
	      cr<max(page_of_boxes[i].x2-page_of_boxes[i].x1, page_of_boxes[i].y2-page_of_boxes[i].y1))
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

      if ((rh<FLT_MAX|| rt<FLT_MAX) && !agent_structure)
	{
	  if (rt<rh)
	    before[j_tail].push_back(i);
	  else
	    before[j_head+1].push_back(i);
	}
      else if (!agent_structure)
	{
	  double rh = FLT_MAX;
	  int j_head=0;
	  for (int j=0; j<arrows.size(); j++)
	    if (arrows[j].linebreak)
	      {
		double r = page_of_boxes[i].y1-arrows[j].head.y;
		if (r>0 && r<rh)
		  {
		    rh = r;
		    j_head = j;
		  }
	      }
	  double rt = FLT_MAX;
	  int j_tail=0;
	  for (int j=0; j<arrows.size(); j++)
	    if (j==0 || arrows[j-1].linebreak)
	      {
		double r = arrows[j].tail.y - page_of_boxes[i].y2;
		if (r>0 && r<rt)
		  {
		    rt = r;
		    j_tail = j;
		  }
	      }		
	  if (rh<FLT_MAX|| rt<FLT_MAX)
	    {
	      if (rt<rh)
		before[j_tail].push_back(i);
	      else
		before[j_head+1].push_back(i);
	    }
	}
    }
 
}

void arrange_reactions(vector<arrow_t> &arrows, const vector<box_t> &page_of_boxes, const vector<point_t> &pluses, vector<string> &results,
		       const vector<string> &page_of_structures,  const string &output_format)
{
  vector < vector<int> > before;
  if (arrows.empty() || page_of_boxes.empty()) return;

  // Find average distance between nearest arrows and standard deviation
  vector<arrow_t> arrows_by_closest;
  vector<double> dist_arrows;
  double avg_dist_arrow=0, avg_dist_arrow_squared=0, std_dev_arrow=0;
  int n_arrows = arrows.size();
  point_t start;
  start.x=0;
  start.y=0;
  
  while (!arrows.empty())
    {
      double d = FLT_MAX;
      int min_i=0;
      for (int i=0; i<arrows.size(); i++)
	if (d>min(distance(start.x,start.y,arrows[i].head.x,arrows[i].head.y),distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y)))
	  {
	    d = min(distance(start.x,start.y,arrows[i].head.x,arrows[i].head.y),distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y));
	    min_i = i;
	  }
      arrows_by_closest.push_back(arrows[min_i]);
      dist_arrows.push_back(d);
      start = arrows[min_i].head;
      avg_dist_arrow += d;
      avg_dist_arrow_squared +=d*d;
      arrows.erase(arrows.begin()+min_i);
    }
  avg_dist_arrow /= n_arrows;
  avg_dist_arrow_squared /= n_arrows;
  std_dev_arrow = sqrt(avg_dist_arrow_squared - avg_dist_arrow*avg_dist_arrow);

  mark_reversible(arrows_by_closest);

  // Break arrows into close-knit groups
  vector < vector <arrow_t> > arrow_groups;
  int i=0;
  while (i<n_arrows)
    {
      vector <arrow_t> group;
      group.push_back(arrows_by_closest[i]);
      i++;
      while (i<n_arrows && dist_arrows[i]<avg_dist_arrow+2*std_dev_arrow)
	{
	  group.push_back(arrows_by_closest[i]);
	  i++;
	}
      arrow_groups.push_back(group);
    }

  // arrange arrows in head to tail fashion
  for (int i=0; i<arrow_groups.size(); i++)
    {
      linear_arrow_sort(arrow_groups[i]);
      check_the_last_arrow_linebreak(arrow_groups[i],page_of_boxes);
    }
  // combine groups into a flat list of arrows
  arrows.clear();
  for (int i=0; i<arrow_groups.size(); i++)
    {
      // arrange structures to best fit between arrows
      vector < vector<int> > before_group(arrow_groups[i].size()+1);
      arrange_structures_between_arrows(arrow_groups[i],before_group,page_of_boxes,page_of_structures);
      sort_boxes_from_arrows(arrow_groups[i],before_group,page_of_boxes);

      for (int j=0; j<arrow_groups[i].size(); j++)
	arrows.push_back(arrow_groups[i][j]);

      if (i!=0)
	{
	  for (int k=0; k<before_group[0].size(); k++)
	    before.back().push_back(before_group[0][k]);
	  for (int j=1; j<before_group.size(); j++)
	    before.push_back(before_group[j]);
	}
      else
	for (int j=0; j<before_group.size(); j++)
	  before.push_back(before_group[j]);

    }

  //   for (int i=0; i<arrows.size(); i++)
  // cout<<arrows[i].tail.x<<","<<arrows[i].tail.y<<" "<<arrows[i].head.x<<","<<arrows[i].head.y<<" "<<arrows[i].linebreak<<endl;

 
  vector < vector <bool> > is_plus(page_of_boxes.size(), vector <bool> (page_of_boxes.size(), false));
  // arrange plus signs between boxes
  for (int i=0; i<before.size(); i++)
    for (int j=1; j<before[i].size(); j++)
      {
	int l = before[i][j-1];
	int k = before[i][j];
	if (k>=0 && l>=0)
	  {
	    box_t a = page_of_boxes[l];
	    box_t b = page_of_boxes[k];
	    for (int m=0; m<pluses.size(); m++)
	      {
		double d=distance_from_bond_y((a.x1+a.x2)/2,(a.y1+a.y2)/2,(b.x1+b.x2)/2,(b.y1+b.y2)/2,pluses[m].x, pluses[m].y);
		if (fabs(d)<min(a.y2-a.y1,b.y2-b.y1)/2 && ((pluses[m].x>a.x2 && pluses[m].x<b.x1) ||  (pluses[m].x>b.x2 && pluses[m].x<a.x1)))
		  {
		    is_plus[k][l] = true;
		    is_plus[l][k] = true;
		  }
		// after plus things can be on the next line
		d = pluses[m].y - (a.y2 + a.y1)/2;
		if (pluses[m].x>a.x2 && fabs(d)<(a.y2-a.y1)/2 && b.y1>a.y2 && pluses[m].x-a.x2<MAX_DISTANCE_BETWEEN_ARROWS)
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
      while (ii>=0 && before[i][ii]<0) ii--;

      if (ii>=0 && before[i][ii]>=0)
	r.push_back(before[i][ii]);

      for (int j=ii-1; j>=0; j--)
	{
	  int k = before[i][j];
	  int l = before[i][j+1];
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
	  while (ii<before[i+1].size() && before[i+1][0]<0) ii++;

	  if (ii<before[i+1].size() && before[i+1][ii]>=0)
	    p.push_back(before[i+1][ii]);
	  for (int j=ii+1; j<before[i+1].size(); j++)
	    {
	      int k = before[i+1][j];
	      int l = before[i+1][j-1];
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
	  string result=convert_page_to_reaction(page_of_structures,output_format, r, p, arrows[i].agent,arrows[i].reversible);
	  trim(result);
	  if (!result.empty())
	    results.push_back(result);
	}
    }
}
