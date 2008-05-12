/*********************************************************************
  OSRA: Optical Structure Recognition
  
  This is a U.S. Government work (2007-2008) and is therefore not subject to 
  copyright.  
  However, portions of this work were obtained from a GPL or GPL-compatiple source.   
  Created by Igor Filippov, 2007-2008 (igorf@helix.nih.gov)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
  USA

*********************************************************************/
#include <Magick++.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <list>

//#include <omp.h>

using namespace std;
extern "C" {
#include "potracelib.h"
#include "pgm2asc.h"
}

#include "openbabel/mol.h"
#include "openbabel/obconversion.h" 
using namespace OpenBabel;

#include "CmdLine.h"

#include <float.h>

#include <algorithm>
#include <cstdio>
#include <vector>

#include "osra.h"


using namespace Magick;
#define PI 3.14159265358979323846


#define BM_WORDSIZE ((int)sizeof(potrace_word))
#define BM_WORDBITS (8*BM_WORDSIZE)
#define BM_HIBIT (((potrace_word)1)<<(BM_WORDBITS-1))
#define bm_scanline(bm, y) ((bm)->map + (y)*(bm)->dy)
#define bm_index(bm, x, y) (&bm_scanline(bm, y)[(x)/BM_WORDBITS])
#define bm_mask(x) (BM_HIBIT >> ((x) & (BM_WORDBITS-1)))
#define bm_range(x, a) ((int)(x) >= 0 && (int)(x) < (a))
#define bm_safe(bm, x, y) (bm_range(x, (bm)->w) && bm_range(y, (bm)->h))
#define BM_USET(bm, x, y) (*bm_index(bm, x, y) |= bm_mask(x))
#define BM_UCLR(bm, x, y) (*bm_index(bm, x, y) &= ~bm_mask(x))
#define BM_UPUT(bm, x, y, b) ((b) ? BM_USET(bm, x, y) : BM_UCLR(bm, x, y))
#define BM_PUT(bm, x, y, b) (bm_safe(bm, x, y) ? BM_UPUT(bm, x, y, b) : 0)

/* return new un-initialized bitmap. NULL with errno on error */
static potrace_bitmap_t *bm_new(int w, int h) {
  potrace_bitmap_t *bm;
  int dy = (w + BM_WORDBITS - 1) / BM_WORDBITS;

  bm = (potrace_bitmap_t *) malloc(sizeof(potrace_bitmap_t));
  if (!bm) {
    return NULL;
  }
  bm->w = w;
  bm->h = h;
  bm->dy = dy;
  bm->map = (potrace_word *) malloc(dy * h * BM_WORDSIZE);
  if (!bm->map) {
    free(bm);
    return NULL;
  }
  return bm;
}



double distance(double x1, double y1, double x2, double y2)
{
  return(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)));
}

double angle4(double x1,double y1, double x2, double y2, 
	      double x3, double y3, double x4, double y4)
{
  double p,l1,l2,cos;

  p=(x1-x2)*(x3-x4)+(y1-y2)*(y3-y4);
  l1=distance(x1,y1,x2,y2);
  l2=distance(x4,y4,x3,y3);
  cos=p/(l1*l2);
  return(cos);
}


void remove_disconnected_atoms(atom_t *atom, bond_t *bond, int n_atom, int n_bond)
{
  for (int i=0;i<n_atom;i++)
      {
        if (atom[i].exists)
          {
            atom[i].exists=false;
            for (int j=0;j<n_bond;j++)
              {
                if ((bond[j].exists) && (i==bond[j].a || i==bond[j].b))
                  {
                    atom[i].exists=true;
                  }
              }
          }
      }
}

void remove_zero_bonds(bond_t *bond, int n_bond,atom_t *atom)
{
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists)
      {
	for (int j=0;j<n_bond;j++)
	  if ((bond[j].exists) && (j!=i) && 
	      ((bond[i].a==bond[j].a && bond[i].b==bond[j].b) ||
	       (bond[i].a==bond[j].b && bond[i].b==bond[j].a)))
	    bond[j].exists=false;
	if (bond[i].a==bond[i].b) bond[i].exists=false;
	if (!atom[bond[i].a].exists || !atom[bond[i].b].exists) bond[i].exists=false;
      }
}



int getPixel(Image image, ColorGray bg,unsigned int x, unsigned int y, double THRESHOLD)
{
  if ((x<image.columns()) && (y<image.rows()))
    {
            ColorGray c=image.pixelColor(x,y);
            if (fabs(c.shade()-bg.shade())>THRESHOLD) return(1);
    }
  return(0);
}



void delete_curve(atom_t *atom,bond_t *bond,int n_atom,int n_bond, potrace_path_t *curve)
{
  for(int i=0;i<n_atom;i++)
    {
      if (atom[i].curve==curve) {atom[i].exists=false;}
    }
  for(int i=0;i<n_bond;i++)
    {
      if (bond[i].curve==curve) {bond[i].exists=false;}
    }
}

void delete_curve_with_children(atom_t *atom,bond_t *bond,int n_atom,int n_bond, 
				potrace_path_t *p)
{
  delete_curve(atom,bond,n_atom,n_bond,p);
  potrace_path_t *child=p->childlist;
  while (child !=NULL)
    {
      delete_curve(atom,bond,n_atom,n_bond,child);
      child=child->sibling;
    }
}


void  delete_bonds_in_char(bond_t *bond,int n_bond,atom_t *atom,
			   double left,double top,double right,double bottom)
{
  for (int j=0;j<n_bond;j++)
    if (bond[j].exists && atom[bond[j].a].x>=left &&
	atom[bond[j].a].x<=right && atom[bond[j].a].y>=top &&
	atom[bond[j].a].y<=bottom && atom[bond[j].b].x>=left &&
	atom[bond[j].b].x<=right && atom[bond[j].b].y>=top &&
	atom[bond[j].b].y<=bottom)
      bond[j].exists=false;
}

double angle_between_bonds(bond_t *bond,int i,int j,atom_t *atom)
{
  return(angle4(atom[bond[i].a].x,atom[bond[i].a].y,atom[bond[i].b].x,atom[bond[i].b].y,
	       atom[bond[j].a].x,atom[bond[j].a].y,atom[bond[j].b].x,atom[bond[j].b].y));
}

double angle_between_connected_bonds(bond_t *bond,int i,int j,atom_t *atom)
{
  double x1,y1,x2,y2,x3,y3,a;
  if (bond[i].a==bond[j].a)
    {
      x1=atom[bond[i].a].x;
      y1=atom[bond[i].a].y;
      x2=atom[bond[i].b].x;
      y2=atom[bond[i].b].y;
      x3=atom[bond[j].b].x;
      y3=atom[bond[j].b].y;
    }
  else if (bond[i].a==bond[j].b)
    {
      x1=atom[bond[i].a].x;
      y1=atom[bond[i].a].y;
      x2=atom[bond[i].b].x;
      y2=atom[bond[i].b].y;
      x3=atom[bond[j].a].x;
      y3=atom[bond[j].a].y;
    }
  else if (bond[i].b==bond[j].a)
    {
      x1=atom[bond[i].b].x;
      y1=atom[bond[i].b].y;
      x2=atom[bond[i].a].x;
      y2=atom[bond[i].a].y;
      x3=atom[bond[j].b].x;
      y3=atom[bond[j].b].y;
    }
  else if (bond[i].b==bond[j].b)
    {
      x1=atom[bond[i].b].x;
      y1=atom[bond[i].b].y;
      x2=atom[bond[i].a].x;
      y2=atom[bond[i].a].y;
      x3=atom[bond[j].a].x;
      y3=atom[bond[j].a].y;
    }
  else {return(0);}
  a=angle4(x2,y2,x1,y1,x3,y3,x1,y1);
  a=acos(a)* 180.0 / PI;
  return(a);
}


double bond_length(bond_t *bond, int i,atom_t *atom)
{
  return(distance(atom[bond[i].a].x,atom[bond[i].a].y,atom[bond[i].b].x,atom[bond[i].b].y));
}



double distance_from_bond_y(double x0,double y0,double x1,double y1,double x,double y)
{
  double d1=distance(x0,y0,x1,y1);
  double cos=(x1-x0)/d1;
  double sin=-(y1-y0)/d1;
  double h=-(x-x0)*sin-(y-y0)*cos;
  return(h);
}

double distance_between_bonds(bond_t *bond,int i,int j,atom_t *atom,double thickness)
{
  double y3=distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,atom[bond[i].b].x,
				 atom[bond[i].b].y,atom[bond[j].a].x,atom[bond[j].a].y);
  double y4=distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,atom[bond[i].b].x,
				 atom[bond[i].b].y,atom[bond[j].b].x,atom[bond[j].b].y);
  if (fabs(y3-y4)>thickness) return(FLT_MAX);
  return(max(fabs(y3),fabs(y4)));
}

double distance_from_bond_x_a(double x0,double y0,double x1,double y1,double x,double y)
{
  double d1=distance(x0,y0,x1,y1);
  double cos=(x1-x0)/d1;
  double sin=-(y1-y0)/d1;
  double l=(x-x0)*cos-(y-y0)*sin;
  return(l);
}

double distance_from_bond_x_b(double x0,double y0,double x1,double y1,double x,double y)
{
  double d1=distance(x0,y0,x1,y1);
  double cos=(x1-x0)/d1;
  double sin=-(y1-y0)/d1;
  double l=(x-x0)*cos-(y-y0)*sin;
  return(l-d1);
}

void bond_end_swap(bond_t *bond, int i)
{
  int t=bond[i].a;
  bond[i].a=bond[i].b;
  bond[i].b=t;
}

bool bonds_within_each_other(bond_t *bond,int ii,int jj,atom_t *atom)
{
  int i,j;
  bool res=false;
  if (bond_length(bond,ii,atom)>bond_length(bond,jj,atom))
    {
      i=ii;
      j=jj;
    }
  else
    {
      i=jj;
      j=ii;
    }
  double x1=atom[bond[i].a].x;
  double x2=atom[bond[i].b].x;
  double y1=atom[bond[i].a].y;
  double y2=atom[bond[i].b].y;
  double d1=bond_length(bond,i,atom);
  double x3=distance_from_bond_x_a(x1,y1,x2,y2,atom[bond[j].a].x,atom[bond[j].a].y);
  double x4=distance_from_bond_x_a(x1,y1,x2,y2,atom[bond[j].b].x,atom[bond[j].b].y);
  if ((x3+x4)/2>0 && (x3+x4)/2<d1) res=true;
  return(res);
}

int num_comp(const void *a,const void *b)
{
  double aa=double(*(double*)a);
  double bb=double(*(double*)b);
  if (aa<bb) return(-1);
  if (aa==bb) return(0);
  if (aa>bb) return(1);
  return(0);
}


double percentile75(bond_t *bond, int n_bond,atom_t *atom)
{
  double a[MAX_ATOMS];
  int n=0;
  for(int i=0;i<n_bond;i++)
    if (bond[i].exists)
      {
	a[n++]=bond_length(bond,i,atom);
      }
  qsort(a,n,sizeof(double),num_comp);
  double pos=0.75*n;
  return(a[int(pos)]);
}

bool alone(bond_t* bond,int i,double avg)
{
  bool alone=false;
  potrace_path_t *p=bond[i].curve;
  if ((p->sign == int('+')) && (p->area<2*avg)) alone=true;
  return(alone);
}


bool no_white_space(int ai,int bi,int aj, int bj, atom_t *atom,Image image,
		    double threshold,ColorGray bgColor)
{
  double xx[4];
  double dx1=atom[bi].x-atom[ai].x;
  double dy1=atom[bi].y-atom[ai].y;
  double dx2=atom[bj].x-atom[aj].x;
  double dy2=atom[bj].y-atom[aj].y;
  double k1,k2;
  int s=0,w=0;

  if (fabs(dx1)>fabs(dy1))
    {
      xx[0]=atom[ai].x;
      xx[1]=atom[bi].x;
      xx[2]=atom[aj].x;
      xx[3]=atom[bj].x;
      qsort(xx,4,sizeof(double),num_comp);
      k1=dy1/dx1;
      k2=dy2/dx2;
      int d=(dx1>0?1:-1);

      for(int x=int(atom[ai].x);x!=int(atom[bi].x);x+=d)
	if(x>xx[1] && x<xx[2])
	  {
	   double p1=(x-atom[ai].x)*k1+atom[ai].y;
	   double p2=(x-atom[aj].x)*k2+atom[aj].y;
	   if (fabs(p2-p1)<1) continue;
	   int dp=(p2>p1?1:-1);

	   for(int y=int(p1)+dp;y!=int(p2);y+=dp)
	     {
	       s++;
	       if (getPixel(image,bgColor,x,y,threshold)==0) w++;
	     }
	  }
    }
  else  
    {
      xx[0]=atom[ai].y;
      xx[1]=atom[bi].y;
      xx[2]=atom[aj].y;
      xx[3]=atom[bj].y;
      qsort(xx,4,sizeof(double),num_comp);
      k1=dx1/dy1;
      k2=dx2/dy2;
      int d=(dy1>0?1:-1);

      for(int y=int(atom[ai].y);y!=int(atom[bi].y);y+=d)
	if(y>xx[1] && y<xx[2])
	  {
	   double p1=(y-atom[ai].y)*k1+atom[ai].x;
	   double p2=(y-atom[aj].y)*k2+atom[aj].x;
	   if (fabs(p2-p1)<1) continue;
	   int dp=(p2>p1?1:-1);
	   for(int x=int(p1)+dp;x!=int(p2);x+=dp)
	     {
	       s++;
	       if (getPixel(image,bgColor,x,y,threshold)==0) w++;
	     }
	  }
    }
  if (s==0) return(true);
  if ((1.*w)/s>WHITE_SPACE_FRACTION) return(false);
  else return(true);
}


double skeletize(atom_t *atom,bond_t *bond,int n_bond,Image image,
		 double threshold,ColorGray bgColor)
{
  double thickness=0;
  double a[MAX_ATOMS];
  int n=0;
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists)
      {
	double l1=bond_length(bond,i,atom);
	for (int j=0;j<n_bond;j++)
	  if (i!=j && bond[j].exists && bonds_within_each_other(bond,i,j,atom))
	    {
	      double tt=distance_between_bonds(bond,i,j,atom,FLT_MAX); 
	      if ((fabs(angle_between_bonds(bond,i,j,atom))>D_T_TOLERANCE 
		   && no_white_space(bond[i].a,bond[i].b,bond[j].a,bond[j].b,atom,
				     image,threshold,bgColor) && tt<MAX_BOND_THICKNESS)
		  || tt<4)
		{
		  double l2=bond_length(bond,j,atom);
		  a[n++]=tt;
		  if (l1<l2)
		    {
		      bond[i].exists=false;
		      bond[j].type=1;
		      if (bond[i].arom) bond[j].arom=true;
		      break;
		    }
		  else
		    {
		      bond[j].exists=false;
		      bond[i].type=1;
		      if (bond[j].arom) bond[i].arom=true;
		    }
		}
	    }
      }
  qsort(a,n,sizeof(double),num_comp);
  thickness=a[n/2];
  return(thickness);
}


int double_triple_bonds(atom_t *atom,bond_t *bond,int n_bond,double avg,int &n_atom,
			double thickness, double &max_dist_double_bond)
{
  double a[MAX_ATOMS];
  int n=0;
  
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists)
      {
	double l1=bond_length(bond,i,atom);
	bond[i].conjoined=false;
	for (int j=i+1;j<n_bond;j++)
	  if ((bond[j].exists) 
	      && (fabs(angle_between_bonds(bond,i,j,atom))>D_T_TOLERANCE))
	    {
	      double l2=bond_length(bond,j,atom);
	      double dbb=distance_between_bonds(bond,i,j,atom,thickness);
	      if (dbb<2*avg/3 && l1>avg/3 && l2>avg/3 && 
		  bonds_within_each_other(bond,i,j,atom))
		a[n++]=dbb;
	    }
      }
  qsort(a,n,sizeof(double),num_comp);
  max_dist_double_bond=a[3*n/4];
  if (max_dist_double_bond<1) max_dist_double_bond=avg/3;
  else max_dist_double_bond+=2;

  for (int i=0;i<n_bond;i++)
    if (bond[i].exists)
      {
	double l1=bond_length(bond,i,atom);
	for (int j=i+1;j<n_bond;j++)
	  if ((bond[j].exists) && (fabs(angle_between_bonds(bond,i,j,atom))>D_T_TOLERANCE))
	    {
	      double l2=bond_length(bond,j,atom);
	      if ((distance_between_bonds(bond,i,j,atom,thickness)<=max_dist_double_bond)
		  && (bonds_within_each_other(bond,i,j,atom))
		  )
		{
		  if (l1>avg && l1>1.5*l2 && l2>0.5*avg)
		    {
		      double aa=distance(atom[bond[i].a].x,atom[bond[i].a].y,
					 atom[bond[j].a].x,atom[bond[j].a].y);
		      double ab=distance(atom[bond[i].a].x,atom[bond[i].a].y,
					 atom[bond[j].b].x,atom[bond[j].b].y);
		      double ba=distance(atom[bond[i].b].x,atom[bond[i].b].y,
					 atom[bond[j].a].x,atom[bond[j].a].y);
		      double bb=distance(atom[bond[i].b].x,atom[bond[i].b].y,
					 atom[bond[j].b].x,atom[bond[j].b].y);
		      double da=min(aa,ab);
		      double db=min(ba,bb);
		      if (da>0.5*l2)
			{
			  double x=atom[bond[i].a].x+
			    (atom[bond[i].b].x-atom[bond[i].a].x)*da/l1;
			  double y=atom[bond[i].a].y+
			    (atom[bond[i].b].y-atom[bond[i].a].y)*da/l1;
			  bond[n_bond].a=bond[i].a;
			  bond[n_bond].exists=true;
			  bond[n_bond].type=1;
			  bond[n_bond].curve=bond[i].curve;
			  bond[n_bond].hash=false;
			  bond[n_bond].wedge=false;
			  bond[n_bond].up=false;
			  bond[n_bond].down=false;
			  bond[n_bond].Small=false;
			  bond[n_bond].arom=false;
			  atom[n_atom].x=x;
			  atom[n_atom].y=y;
			  atom[n_atom].label=" ";
			  atom[n_atom].exists=true;
			  atom[n_atom].curve=bond[i].curve;
			  atom[n_atom].n=0;
			  atom[n_atom].corner=false;
			  bond[i].a=n_atom;
			  n_atom++;
			  if (n_atom>=MAX_ATOMS) n_atom--;
			  bond[n_bond].b=bond[i].a;
			  n_bond++;
			  if (n_bond>=MAX_ATOMS) n_bond--;
			}
		      if (db>0.5*l2)
			{
			  double x=atom[bond[i].b].x+
			    (atom[bond[i].a].x-atom[bond[i].b].x)*db/l1;
			  double y=atom[bond[i].b].y+
			    (atom[bond[i].a].y-atom[bond[i].b].y)*db/l1;
			  bond[n_bond].a=bond[i].b;
			  bond[n_bond].exists=true;
			  bond[n_bond].type=1;
			  bond[n_bond].curve=bond[i].curve;
			  bond[n_bond].hash=false;
			  bond[n_bond].wedge=false;
			  bond[n_bond].up=false;
			  bond[n_bond].down=false;
			  bond[n_bond].Small=false;
			  bond[n_bond].arom=false;
			  atom[n_atom].x=x;
			  atom[n_atom].y=y;
			  atom[n_atom].label=" ";
			  atom[n_atom].exists=true;
			  atom[n_atom].curve=bond[i].curve;
			  atom[n_atom].n=0;
			  atom[n_atom].corner=false;
			  bond[i].b=n_atom;
			  n_atom++;
			  if (n_atom>=MAX_ATOMS) n_atom--;
			  bond[n_bond].b=bond[i].b;
			  n_bond++;
			  if (n_bond>=MAX_ATOMS) n_bond--;
			}
		      bond[j].exists=false;
		      bond[i].type+=bond[j].type;
		      if (bond[j].arom) bond[i].arom=true;
		      if (bond[j].curve==bond[i].curve) bond[i].conjoined=true;
		    }
		  else if (l2>avg && l2>1.5*l1 && l1>0.5*avg)
		    {
		      double aa=distance(atom[bond[j].a].x,atom[bond[j].a].y,
					 atom[bond[i].a].x,atom[bond[i].a].y);
		      double ab=distance(atom[bond[j].a].x,atom[bond[j].a].y,
					 atom[bond[i].b].x,atom[bond[i].b].y);
		      double ba=distance(atom[bond[j].b].x,atom[bond[j].b].y,
					 atom[bond[i].a].x,atom[bond[i].a].y);
		      double bb=distance(atom[bond[j].b].x,atom[bond[j].b].y,
					 atom[bond[i].b].x,atom[bond[i].b].y);
		      double da=min(aa,ab);
		      double db=min(ba,bb);
		      if (da>0.5*l1)
			{
			  double x=atom[bond[j].a].x+
			    (atom[bond[j].b].x-atom[bond[j].a].x)*da/l2;
			  double y=atom[bond[j].a].y+
			    (atom[bond[j].b].y-atom[bond[j].a].y)*da/l2;
			  bond[n_bond].a=bond[j].a;
			  bond[n_bond].exists=true;
			  bond[n_bond].type=1;
			  bond[n_bond].curve=bond[j].curve;
			  bond[n_bond].hash=false;
			  bond[n_bond].wedge=false;
			  bond[n_bond].up=false;
			  bond[n_bond].down=false;
			  bond[n_bond].Small=false;
			  bond[n_bond].arom=false;
			  atom[n_atom].x=x;
			  atom[n_atom].y=y;
			  atom[n_atom].label=" ";
			  atom[n_atom].exists=true;
			  atom[n_atom].curve=bond[j].curve;
			  atom[n_atom].n=0;
			  atom[n_atom].corner=false;
			  bond[j].a=n_atom;
			  n_atom++;
			  if (n_atom>=MAX_ATOMS) n_atom--;
			  bond[n_bond].b=bond[j].a;
			  n_bond++;
			  if (n_bond>=MAX_ATOMS) n_bond--;
			}
		      if (db>0.5*l1)
			{
			  double x=atom[bond[j].b].x+
			    (atom[bond[j].a].x-atom[bond[j].b].x)*db/l2;
			  double y=atom[bond[j].b].y+
			    (atom[bond[j].a].y-atom[bond[j].b].y)*db/l2;
			  bond[n_bond].a=bond[j].b;
			  bond[n_bond].exists=true;
			  bond[n_bond].type=1;
			  bond[n_bond].curve=bond[j].curve;
			  bond[n_bond].hash=false;
			  bond[n_bond].wedge=false;
			  bond[n_bond].up=false;
			  bond[n_bond].down=false;
			  bond[n_bond].Small=false;
			  bond[n_bond].arom=false;
			  atom[n_atom].x=x;
			  atom[n_atom].y=y;
			  atom[n_atom].label=" ";
			  atom[n_atom].exists=true;
			  atom[n_atom].curve=bond[j].curve;
			  atom[n_atom].n=0;
			  atom[n_atom].corner=false;
			  bond[j].b=n_atom;
			  n_atom++;
			  if (n_atom>=MAX_ATOMS) n_atom--;
			  bond[n_bond].b=bond[j].b;
			  n_bond++;
			  if (n_bond>=MAX_ATOMS) n_bond--;
			}
		      bond[i].exists=false;
		      bond[j].type+=bond[i].type;
		      if (bond[i].arom) bond[j].arom=true;
		      if (bond[j].curve==bond[i].curve) bond[j].conjoined=true;
		      break;
		    }
		  else
		    {
		      if (l1>l2)
			{
			  bond[j].exists=false;
			  if (l2>l1/2)
			    {  
			      bond[i].type+=bond[j].type;
			      if (bond[j].curve==bond[i].curve) bond[i].conjoined=true;
			    }
			  if (bond[j].arom) bond[i].arom=true;
			}
		      else
			{
			  bond[i].exists=false;
			  if (l1>l2/2)  
			    {
			      bond[j].type+=bond[i].type;
			      if (bond[j].curve==bond[i].curve) bond[j].conjoined=true;
			    }
			  if (bond[i].arom) bond[j].arom=true;
			  break;
			}
		    }
		}
	    }
      }
  return(n_bond);
}

bool chlorine(bond_t *bond, atom_t *atom,int i, letters_t *letters,int n_letters, 
	     int max_font_height,int min_font_height)
{
  bool res=false;
  double x=(atom[bond[i].a].x+atom[bond[i].b].x)/2;
  double y=(atom[bond[i].a].y+atom[bond[i].b].y)/2;
  double r=bond_length(bond,i,atom)/2;
  if ((bond_length(bond,i,atom)<max_font_height) &&
      (bond_length(bond,i,atom)>min_font_height) &&
      (fabs(atom[bond[i].a].x-atom[bond[i].b].x)<fabs(atom[bond[i].a].y-atom[bond[i].b].y)))
    {
      for (int j=0;j<n_letters;j++)
	{
	  if ((distance(x,y,letters[j].x,letters[j].y)<r+letters[j].r) &&
	      (fabs(y-letters[j].y)<min(r,letters[j].r)))
	    { 
	      res=true;
	    }
	}
    }
	 
  return(res);
}

int remove_small_bonds(bond_t *bond, int n_bond,atom_t *atom, 
			letters_t *letters, int n_letters, int max_font_height,
		       int min_font_height,double avg)
{
  for (int i=0;i<n_bond;i++)
    if ((bond[i].exists) && (bond[i].type==1))
      {
	bool al=alone(bond,i,avg);
	if (bond_length(bond,i,atom)<V_DISPLACEMENT)
	  {
	    bond[i].exists=false;
	  }
	else if ((al) && (chlorine(bond,atom,i,letters,n_letters,max_font_height,min_font_height)))
	  {
	    letters[n_letters].a='l';
	    letters[n_letters].x=(atom[bond[i].a].x+atom[bond[i].b].x)/2;
	    letters[n_letters].y=(atom[bond[i].a].y+atom[bond[i].b].y)/2;
	    letters[n_letters].r=bond_length(bond,i,atom)/2;
	    letters[n_letters].free=true;
	    n_letters++;
	    if (n_letters>=MAX_ATOMS) n_letters--;
	    bond[i].exists=false;
	  }
      }
  return(n_letters);
} 


int comp_lbonds(const void *a,const void *b)
{
  lbond_t *aa=(lbond_t *) a;
  lbond_t *bb=(lbond_t *) b;
  if (aa->x<bb->x) return(-1);
  if (aa->x==bb->x) return(0);
  if (aa->x>bb->x) return(1);
  return(0);
}

int comp_letters(const void *a,const void *b)
{
  letters_t *aa=(letters_t *) a;
  letters_t *bb=(letters_t *) b;
  if (aa->x<bb->x) return(-1);
  if (aa->x==bb->x) return(0);
  if (aa->x>bb->x) return(1);
  return(0);
}

bool terminal_bond(int a,int b,bond_t *bond,int n_bond)
{
  bool terminal=true;
  for (int l=0;l<n_bond;l++)
    if (l!=b && bond[l].exists && 
	(bond[l].a==a || bond[l].b==a))
      terminal=false;
  return(terminal);
}

int assemble_labels(letters_t *letters,int n_letters,label_t *label)
{
  lbond_t lbond[MAX_ATOMS];
  int n_lbond=0;
  int n_label=0;
  qsort(letters,n_letters,sizeof(letters_t),comp_letters);
  for (int i=0;i<n_letters;i++)
    {
     for (int j=i+1;j<n_letters;j++)
       if (
	   (distance(letters[i].x,letters[i].y,letters[j].x,letters[j].y)<2*max(letters[i].r,letters[j].r) && 
	   (((fabs(letters[i].y-letters[j].y)<min(letters[i].r,letters[j].r))) ||
	    ((fabs(letters[i].y-letters[j].y)<(letters[i].r+letters[j].r)) &&
	     (((letters[i].y<letters[j].y) && (isdigit(letters[j].a))) ||
	      ((letters[j].y<letters[i].y) && (isdigit(letters[i].a))))))) ||
	   (distance(letters[i].x,letters[i].y,letters[j].x,letters[j].y)<1.5*(letters[i].r+letters[j].r) && 
	    (letters[i].a=='-' || letters[i].a=='+' || letters[j].a=='-' || letters[j].a=='+')))
	{
	  lbond[n_lbond].a=i;
	  lbond[n_lbond].b=j;
	  lbond[n_lbond].x=letters[i].x;

	  letters[i].free=false;
	  letters[j].free=false;
	  lbond[n_lbond].exists=true;
	  n_lbond++;
	  if (n_lbond>=MAX_ATOMS) n_lbond--;
	  break;
	}
    }
  qsort(lbond,n_lbond,sizeof(lbond_t),comp_lbonds);
  
  for (int i=0;i<n_lbond;i++)
    if (lbond[i].exists)
      {
	bool found_left=false;
	label[n_label].x1=FLT_MAX;
	label[n_label].y1=FLT_MAX;
	label[n_label].r1=0;
	label[n_label].x2=FLT_MAX;
	label[n_label].y2=FLT_MAX;
	label[n_label].r2=0;
	label[n_label].a=letters[lbond[i].a].a;
	label[n_label].a+=letters[lbond[i].b].a;
	if (!isdigit(letters[lbond[i].a].a) && letters[lbond[i].a].a!='-'
	    && letters[lbond[i].a].a!='+' && !found_left)
	  {
	    label[n_label].x1=letters[lbond[i].a].x;
	    label[n_label].y1=letters[lbond[i].a].y;
	    label[n_label].r1=letters[lbond[i].a].r;
	    found_left=true;
	  }
	if (!isdigit(letters[lbond[i].b].a) && letters[lbond[i].b].a!='-'
	    && letters[lbond[i].b].a!='+' && !found_left)
	  {
	    label[n_label].x1=letters[lbond[i].b].x;
	    label[n_label].y1=letters[lbond[i].b].y;
	    label[n_label].r1=letters[lbond[i].b].r;
	    found_left=true;
	  }
	if (!isdigit(letters[lbond[i].a].a) && letters[lbond[i].a].a!='-'
	    && letters[lbond[i].a].a!='+')
	  {
	    label[n_label].x2=letters[lbond[i].a].x;
	    label[n_label].y2=letters[lbond[i].a].y;
	    label[n_label].r2=letters[lbond[i].a].r;
	  }
	if (!isdigit(letters[lbond[i].b].a) && letters[lbond[i].b].a!='-'
	    && letters[lbond[i].b].a!='+')
	  {
	    label[n_label].x2=letters[lbond[i].b].x;
	    label[n_label].y2=letters[lbond[i].b].y;
	    label[n_label].r2=letters[lbond[i].b].r;
	  }
	lbond[i].exists=false;
	int last=lbond[i].b;
	 for (int j=i+1;j<n_lbond;j++)
	   if ((lbond[j].exists) && (lbond[j].a==last))
	     {
		label[n_label].a+=letters[lbond[j].b].a;
		if (!isdigit(letters[lbond[j].a].a) && letters[lbond[j].a].a!='-'
		    && letters[lbond[j].a].a!='+' && !found_left)
		  {
		    label[n_label].x1=letters[lbond[j].a].x;
		    label[n_label].y1=letters[lbond[j].a].y;
		    label[n_label].r1=letters[lbond[j].a].r;
		    found_left=true;
		  }
		if (!isdigit(letters[lbond[j].b].a) && letters[lbond[j].b].a!='-'
		    && letters[lbond[j].b].a!='+' && !found_left)
		  {
		    label[n_label].x1=letters[lbond[j].b].x;
		    label[n_label].y1=letters[lbond[j].b].y;
		    label[n_label].r1=letters[lbond[j].b].r;
		    found_left=true;
		  }
		if (!isdigit(letters[lbond[j].a].a) && letters[lbond[j].a].a!='-'
		    && letters[lbond[j].a].a!='+')
		  {
		    label[n_label].x2=letters[lbond[j].a].x;
		    label[n_label].y2=letters[lbond[j].a].y;
		    label[n_label].r2=letters[lbond[j].a].r;
		  }
		if (!isdigit(letters[lbond[j].b].a) && letters[lbond[j].b].a!='-'
		    && letters[lbond[j].b].a!='+')
		  {
		    label[n_label].x2=letters[lbond[j].b].x;
		    label[n_label].y2=letters[lbond[j].b].y;
		    label[n_label].r2=letters[lbond[j].b].r;
		  }
		last=lbond[j].b;
		lbond[j].exists=false;
	     }

	 bool cont=true;
	 string charges="";
	 while (cont)
	  {
	    cont=false;
	    string::size_type pos=label[n_label].a.find_first_of('-');
	    if (pos!=string::npos)
	      {
		label[n_label].a.erase(pos,1);
		charges+="-";
		cont=true;
	      }
	    pos=label[n_label].a.find_first_of('+');
	    if (pos!=string::npos)
	      {
		label[n_label].a.erase(pos,1);
		charges+="+";
		cont=true;
	      }
	  }
	 label[n_label].a+=charges;
	 n_label++;
	 if (n_label>=MAX_ATOMS) n_label--;
      }
  return(n_label);
}

void extend_terminal_bond_to_label(atom_t *atom,letters_t *letters,int n_letters,
			  bond_t *bond, int n_bond,label_t *label,int n_label,
				   double avg, double maxh, double max_dist_double_bond)
{
 for (int j=0;j<n_bond;j++)
    if (bond[j].exists)
      {
	bool not_corner_a=terminal_bond(bond[j].a,j,bond,n_bond);
	bool not_corner_b=terminal_bond(bond[j].b,j,bond,n_bond);
	if (atom[bond[j].a].label!=" ") not_corner_a=false;
	if (atom[bond[j].b].label!=" ") not_corner_b=false;
	double xa=atom[bond[j].a].x;
	double ya=atom[bond[j].a].y;
	double xb=atom[bond[j].b].x;
	double yb=atom[bond[j].b].y;
	double bl=bond_length(bond,j,atom);
	double minb=FLT_MAX;
	bool found1=false,found2=false;
	int l1=-1,l2=-1;
	if (not_corner_a)
	  {
	    for (int i=0;i<n_label;i++)
	      if ((label[i].a)[0]!='+' && (label[i].a)[0]!='-')
		{
		  double d1=distance_from_bond_x_a(xa,ya,xb,yb,
						   label[i].x1,label[i].y1);
		  double d2=distance_from_bond_x_a(xa,ya,xb,yb,
						   label[i].x2,label[i].y2);
		  double h1=fabs(distance_from_bond_y(xa,ya,xb,yb,
						      label[i].x1,label[i].y1));
		  double h2=fabs(distance_from_bond_y(xa,ya,xb,yb,
						      label[i].x2,label[i].y2));
		  double y_dist=maxh;
		  if (bond[j].type==2) y_dist+=max_dist_double_bond;
		  if (fabs(d1)<fabs(d2))   // 1st side of label closer
		    {
		      double nb=fabs(d1)-label[i].r1;
		      if (nb<=avg && h1<=y_dist+label[i].r1 && nb<minb && d1<bl/2)
			{
			  found1=true;
			  l1=i;
			  minb=nb;
			}
		    }
		  else         // 2nd side closer
		    {
		      double nb=fabs(d2)-label[i].r2;
		      if (nb<=avg && h2<=y_dist+label[i].r2 && nb<minb && d2<bl/2)
			{
			  found1=true;
			  l1=i;
			  minb=nb;
			}
		    }
		}
	    for (int i=0;i<n_letters;i++)
	      if (letters[i].free && letters[i].a!='+' && letters[i].a!='-')
		{
		  double d=distance_from_bond_x_a(xa,ya,xb,yb,letters[i].x,letters[i].y);
		  double y_dist=maxh;
		  if (bond[j].type==2) y_dist+=max_dist_double_bond;
		  double h=fabs(distance_from_bond_y(xa,ya,xb,yb,
						     letters[i].x,letters[i].y));
		  double nb=fabs(d)-letters[i].r;
		  if (nb<=avg && h<=y_dist+letters[i].r && nb<minb && d<bl/2)
		    {
		      found2=true;
		      l2=i;
		      minb=nb;
		    }
		}
	    if (found2)
	      {
		atom[bond[j].a].label=toupper(letters[l2].a);;
		atom[bond[j].a].x=letters[l2].x;
		atom[bond[j].a].y=letters[l2].y;
	      }
	    else if (found1)
	      {
		atom[bond[j].a].label=label[l1].a;
		atom[bond[j].a].x=(label[l1].x1+label[l1].x2)/2;
		atom[bond[j].a].y=(label[l1].y1+label[l1].y2)/2;
	      }
	  }
	if (not_corner_b)
	  {
	    found1=false,found2=false;
	    minb=FLT_MAX;
	    for (int i=0;i<n_label;i++)
	      if ((label[i].a)[0]!='+' && (label[i].a)[0]!='-' && i!=l1)
		{
		  double d1=distance_from_bond_x_b(xa,ya,xb,yb,label[i].x1,label[i].y1);
		  double d2=distance_from_bond_x_b(xa,ya,xb,yb,label[i].x2,label[i].y2);
		  double h1=fabs(distance_from_bond_y(xa,ya,xb,yb,
						      label[i].x1,label[i].y1));
		  double h2=fabs(distance_from_bond_y(xa,ya,xb,yb,
						      label[i].x2,label[i].y2));
		  double y_dist=maxh;

		  if (bond[j].type==2) y_dist+=max_dist_double_bond;
		  if (fabs(d1)<fabs(d2))   // 1st side of label closer
		    {
		      double nb=fabs(d1)-label[i].r1;  // end "b" and 1st side
		      y_dist+=label[i].r1;
		      if (nb<=avg && h1<=y_dist && nb<minb && d1>-bl/2)
			{
			  found1=true;
			  l1=i;
			  minb=nb;
			}
		    }
		  else         // 2nd side closer
		    {
		      double nb=fabs(d2)-label[i].r2;   // end "b" and 2nd side
		      y_dist+=label[i].r2;
		      if (nb<=avg && h2<=y_dist && nb<minb && d2>-bl/2)
			{
			  found1=true;
			  l1=i;
			  minb=nb;
			}
		    }
		}
	    for (int i=0;i<n_letters;i++)
	      if (letters[i].free && letters[i].a!='+' && letters[i].a!='-' && i!=l2)
		{
		  double d=distance_from_bond_x_b(xa,ya,xb,yb,letters[i].x,letters[i].y);
		  double nb=fabs(d)-letters[i].r; // distance between end "b" and letter
		  double y_dist=maxh+letters[i].r;
		  if (bond[j].type==2) y_dist+=max_dist_double_bond;
		  double h=fabs(distance_from_bond_y(xa,ya,xb,yb,
						     letters[i].x,letters[i].y));
		  if (nb<=avg && h<=y_dist && nb<minb && d>-bl/2)
		    {
		      found2=true;
		      l2=i;
		      minb=nb;
		    }
		}

	    if (found2)
	      {
		atom[bond[j].b].label=toupper(letters[l2].a);;
		atom[bond[j].b].x=letters[l2].x;
		atom[bond[j].b].y=letters[l2].y;
	      }
	    else if (found1)
	      {
		atom[bond[j].b].label=label[l1].a;
		atom[bond[j].b].x=(label[l1].x1+label[l1].x2)/2;
		atom[bond[j].b].y=(label[l1].y1+label[l1].y2)/2;
	      }
	  }
      }
}

void extend_terminal_bond_to_bonds(atom_t *atom,int n_atom,bond_t *bond, int n_bond,
				   double avg, double maxh,double max_dist_double_bond)
{
 for (int j=0;j<n_bond;j++)
    if (bond[j].exists)
      {
	bool not_corner_a=terminal_bond(bond[j].a,j,bond,n_bond);
	bool not_corner_b=terminal_bond(bond[j].b,j,bond,n_bond);
	double xa=atom[bond[j].a].x;
	double ya=atom[bond[j].a].y;
	double xb=atom[bond[j].b].x;
	double yb=atom[bond[j].b].y;
	double bl=bond_length(bond,j,atom);
	double minb=FLT_MAX;
	bool found=false;
	int l=-1;
	if (not_corner_a)
	  for (int i=0;i<n_atom;i++)
	    if (atom[i].exists && i!=bond[j].a && i!=bond[j].b)
	      {
		double d=distance_from_bond_x_a(xa,ya,xb,yb,atom[i].x,atom[i].y);
		double h=fabs(distance_from_bond_y(xa,ya,xb,yb,atom[i].x,atom[i].y));
		double y_dist=maxh;
		if (bond[j].type==2 && !bond[j].conjoined) y_dist+=max_dist_double_bond;
		if (fabs(d)<=avg/2 && h<=y_dist && fabs(d)<minb && d<bl/2)
		  {
		    found=true;
		    l=i;
		    minb=fabs(d);
		  }
	      }	      
	if (found)
	      {
		atom[l].x=(atom[bond[j].a].x+atom[l].x)/2;
		atom[l].y=(atom[bond[j].a].y+atom[l].y)/2;
		bond[j].a=l;
	      }
	found=false;
	if (not_corner_b)
	  for (int i=0;i<n_atom;i++)
	    if (atom[i].exists && i!=l && i!=bond[j].a && i!=bond[j].b)
	      {
		double d=distance_from_bond_x_b(xa,ya,xb,yb,atom[i].x,atom[i].y);
		double h=fabs(distance_from_bond_y(xa,ya,xb,yb,atom[i].x,atom[i].y));
		double y_dist=maxh;
		if (bond[j].type==2 && !bond[j].conjoined) y_dist+=max_dist_double_bond;
		if (fabs(d)<=avg/2 && h<=y_dist && fabs(d)<minb && d>-bl/2)
		  {
		    found=true;
		    l=i;
		    minb=fabs(d);
		  }
	      }	      
	if (found)
	      {
		atom[l].x=(atom[bond[j].b].x+atom[l].x)/2;
		atom[l].y=(atom[bond[j].b].y+atom[l].y)/2;
		bond[j].b=l;
	      }

      }
}



void assign_charge(atom_t *atom, bond_t *bond, int n_atom,int n_bond)
{
  for (int j=0;j<n_bond;j++)
    if (bond[j].exists && (!atom[bond[j].a].exists || !atom[bond[j].b].exists))
      bond[j].exists=false;

  for (int i=0;i<n_atom;i++)
    if (atom[i].exists)
      {
	int n=0;
	int m=0;
	for (int j=0;j<n_bond;j++)
	  if (bond[j].exists && (bond[j].a==i || bond[j].b==i))
	    {
	      n+=bond[j].type;
	      if (bond[j].type>1) m++;
	    }
	atom[i].charge=0;
	bool cont=true;
	while (cont)
	  {
	    string::size_type pos=atom[i].label.find_first_of('-');
	    if (pos!=string::npos)
	      {
		atom[i].label.erase(pos,1);
		if (atom[i].label.length()>0 && isalpha(atom[i].label.at(0))) 
		  atom[i].charge--;
	      }
	    else
	      {
		pos=atom[i].label.find_first_of('+');
		if (pos!=string::npos)
		  {
		    atom[i].label.erase(pos,1);
		    if (atom[i].label.length()>0 && isalpha(atom[i].label.at(0))) 
		      atom[i].charge++;
		  }
		else cont=false;
	      }
	  }
        for (int j=0;j<n_bond;j++)
          if (bond[j].exists && bond[j].hash && bond[j].b==i)
            atom[i].charge=0;
                                                                                  
	atom[i].label=fix_atom_name(atom[i].label,n);
      }
}

Color getBgColor(Image image,bool inv)
{
  ColorGray c,r;
  r=image.pixelColor(1,1);
  for (int i=0;i<BG_PICK_POINTS;i++)
    {
      int x=(image.columns()*rand())/RAND_MAX;
      int y=(image.rows()*rand())/RAND_MAX;
      c=image.pixelColor(x,y);
      if ((!inv) && ((c.shade())>(r.shade()))) r=c;
      else if ((inv) && ((c.shade())<(r.shade()))) r=c;
    }
  return(r);
}

void debug(Image image,atom_t *atom, int n_atom,bond_t *bond,int n_bond, string fname)
{
  image.modifyImage();
  image.type(TrueColorType);
  image.strokeWidth(1);

 int max_x=image.columns();
 int max_y=image.rows();

 for (int i=0;i<n_bond;i++)
   {
     if ((bond[i].exists) && (atom[bond[i].a].exists) && (atom[bond[i].b].exists))
       {
	 if (bond[i].type==1)
	   {
	     image.strokeColor("green");
	   }
	 else if (bond[i].type==2)
	   {
	     image.strokeColor("yellow");
	   }
	 else if (bond[i].type>=3)
	   {
	     image.strokeColor("red");
	   }
	 if (bond[i].hash)
	   {
	     image.strokeColor("blue");
	   }
	 else if (bond[i].wedge)
	   {
	     image.strokeColor("purple");
	   }
	 image.draw( DrawableLine(atom[bond[i].a].x,atom[bond[i].a].y,atom[bond[i].b].x,atom[bond[i].b].y));
       }
   }


 for (int i=0;i<n_atom;i++)
      {
	if (atom[i].exists)
	  {
	    if ((int(atom[i].x)<max_x) && (int(atom[i].y<max_y)))
	    image.pixelColor(int(atom[i].x), int(atom[i].y),"blue");
	  }
      }
      

  image.write(fname);
}

void draw_square(Image *image,int x1, int y1, int x2, int y2,string color)
{
  image->strokeWidth(1);
  image->strokeColor(color);
  image->draw( DrawableLine(x1,y1,x2,y1));
  image->draw( DrawableLine(x1,y2,x2,y2));
  image->draw( DrawableLine(x1,y1,x1,y2));
  image->draw( DrawableLine(x2,y1,x2,y2));
}

void draw_box(Image image,box_t *boxes,int n_boxes, string fname)
{
  image.modifyImage();
  image.type(TrueColorType);
  for (int i=0;i<n_boxes;i++)
    {
      draw_square(&image,boxes[i].x1,boxes[i].y1,boxes[i].x2,boxes[i].y2,"green");
    }
  image.write(fname);
}

int next_atom(int cur, int begin, int total)
{
  int n=cur+1;
  if (n>total-1) {n=begin;}
  return(n);
}

bool dir_change(int n, int last,int begin, int total, atom_t *atom)
{ 
  int m=next_atom(n,begin,total);
  while (distance(atom[m].x,atom[m].y,atom[n].x,atom[n].y)<V_DISPLACEMENT && m!=n)
    m=next_atom(m,begin,total);
  if (m==n) return(false);
  double s=fabs(distance_from_bond_y(atom[n].x,atom[n].y,atom[last].x,atom[last].y,
				     atom[m].x,atom[m].y));
  if (s>V_DISPLACEMENT) return(true);
  return(false);
}

bool smaller_distance(int n, int last,int begin, int total, atom_t *atom)
{ 
  int m=next_atom(n,begin,total);
  double d1=distance(atom[n].x,atom[n].y,atom[last].x,atom[last].y);
  double d2=distance(atom[m].x,atom[m].y,atom[last].x,atom[last].y);
  if (d1>d2) {return(true);}
  return(false);
}


int find_bonds(atom_t *atom, bond_t *bond, int b_atom, int n_atom, int n_bond,potrace_path_t * p)
{
  int i=b_atom+1;
  int last=b_atom;
  while (i<n_atom)
    {
      if (atom[i].corner) 
      	{
      	  atom[i].exists=true;
      	  last=i;
      	  i++;
      	}
      else if (dir_change(i,last,b_atom,n_atom,atom)) 
      	{
      	  atom[i].exists=true;
      	  last=i;
      	  i++;
      	}
      else if (smaller_distance(i,last,b_atom,n_atom,atom))
      	{
	  atom[i].exists=true;
	  last=i;
	  i++;
      	}
      else {i++;}
    }
  for(i=b_atom;i<n_atom;i++)
    if (atom[i].exists)
      {
	bond[n_bond].a=i;
	bond[n_bond].exists=true;
	bond[n_bond].type=1;
	int j=next_atom(i,b_atom,n_atom);
	while (!atom[j].exists) {j=next_atom(j,b_atom,n_atom);}
	bond[n_bond].b=j;
	bond[n_bond].curve=p;
	bond[n_bond].hash=false;
	bond[n_bond].wedge=false;
	bond[n_bond].up=false;
	bond[n_bond].down=false;
	bond[n_bond].Small=false;
      	n_bond++;
	if (n_bond>=MAX_ATOMS) n_bond--;
      }
    return(n_bond);
}

void adjust_box(Image image,double THRESHOLD_BOND,ColorGray bgColor,
		int width,int height,int boundary,int allowed,
		int &top,int &left, int &bottom, int &right,int maxh,int maxw)
{
  int oldtop=top-1;
  int oldbottom=bottom-1;
  int oldleft=left-1;
  int oldright=right-1;
  while (((oldtop!=top) || (oldbottom!=bottom) || (oldleft!=left) || (oldright!=right)) &&
	 ((right-left)<=maxw) && ((bottom-top)<=maxh))
    {
      oldtop=top;
      oldbottom=bottom;
      oldleft=left;
      oldright=right;
      int s=1;
      while((top>0) && (s>allowed)  && ((right-left)<=maxw) && ((bottom-top)<=maxh))
	{
	  s=0;
	  for (int i=top;i<top+boundary;i++)
	    {
	      for (int j=left;j<=right;j++) s+=getPixel(image,bgColor,j,i,THRESHOLD_BOND);
	    }
	  if (s>allowed) top--;
	}
      s=1;
      while((bottom<height) && (s>allowed) && ((right-left)<=maxw) && ((bottom-top)<=maxh))
	{
	  s=0;
	  for (int i=bottom;i>bottom-boundary;i--)
	    {
	      for (int j=left;j<=right;j++) s+=getPixel(image,bgColor,j,i,THRESHOLD_BOND);
	    }
	  if (s>allowed) bottom++;
	}
      s=1;
      while((left>0) && (s>allowed) && ((right-left)<=maxw) && ((bottom-top)<=maxh))
	{
	  s=0;
	  for (int i=top;i<=bottom;i++)
	    {
	      for (int j=left;j<left+boundary;j++) s+=getPixel(image,bgColor,j,i,THRESHOLD_BOND);
	    }
	  if (s>allowed) left--;
	}
      s=1;
      while((right<width) && (s>allowed) && ((right-left)<=maxw) && ((bottom-top)<=maxh))
	{
	  s=0;
	  for (int i=top;i<=bottom;i++)
	    {
	      for (int j=right;j>right-boundary;j--) s+=getPixel(image,bgColor,j,i,THRESHOLD_BOND);
	    }
	  if (s>allowed) right++;
	}
    }

}


int find_chars(potrace_path_t *p,Image orig,letters_t *letters,
	       atom_t *atom,bond_t *bond,int n_atom,int n_bond,int height,int width,
	       ColorGray bgColor, double THRESHOLD, 
	       int max_font_width, int max_font_height,
	       int &real_font_width, int &real_font_height)
{
  int n, *tag,n_letters=0;
  potrace_dpoint_t (*c)[3];
  real_font_width=0;
  real_font_height=0;

  while (p != NULL) 
      {
	if ((p->sign == int('+')))
	  {
	    n = p->curve.n;
	    tag = p->curve.tag;
	    c = p->curve.c;
	    int top=height;
	    int x1=0;
	    int left=width;
	    int y1=0;
	    int bottom=0;
	    int x2=0;
	    int right=0;
	    int y2=0;
	    for (int i=0; i<n; i++) 
	      {
		switch (tag[i]) 
		  {
		  case POTRACE_CORNER:
		    if (c[i][1].x<left) {left=int(c[i][1].x);y1=int(c[i][1].y);}
		    if (c[i][1].x>right) {right=int(c[i][1].x);y2=int(c[i][1].y);}
		    if (c[i][1].y<top) {top=int(c[i][1].y);x1=int(c[i][1].x);}
		    if (c[i][1].y>bottom) {bottom=int(c[i][1].y);x2=int(c[i][1].x);}
		    break;
		  case POTRACE_CURVETO:
		    if (c[i][0].x<left) {left=int(c[i][0].x);y1=int(c[i][0].y);}
		    if (c[i][0].x>right) {right=int(c[i][0].x);y2=int(c[i][0].y);}
		    if (c[i][0].y<top) {top=int(c[i][0].y);x1=int(c[i][0].x);}
		    if (c[i][0].y>bottom) {bottom=int(c[i][0].y);x2=int(c[i][0].x);}
		    if (c[i][1].x<left) {left=int(c[i][1].x);y1=int(c[i][1].y);}
		    if (c[i][1].x>right) {right=int(c[i][1].x);y2=int(c[i][1].y);}
		    if (c[i][1].y<top) {top=int(c[i][1].y);x1=int(c[i][1].x);}
		    if (c[i][1].y>bottom) {bottom=int(c[i][1].y);x2=int(c[i][1].x);}
		    break;
		  }
		if (c[i][2].x<left) {left=int(c[i][2].x);y1=int(c[i][2].y);}
		if (c[i][2].x>right) {right=int(c[i][2].x);y2=int(c[i][2].y);}
		if (c[i][2].y<top) {top=int(c[i][2].y);x1=int(c[i][2].x);}
		if (c[i][2].y>bottom) {bottom=int(c[i][2].y);x2=int(c[i][2].x);}
	      }

	    if (((bottom-top)<=2*max_font_height) && 
		((right-left)<=2*max_font_width) && (right-left>V_DISPLACEMENT) 
		&& (bottom-top>MIN_FONT_HEIGHT))
	      {
		int s=1;
		while((top>0) && (s>0))
		  {
		    s=0;
		    s=getPixel(orig,bgColor,x1,top,THRESHOLD);
		    if (s>0) top--;
		  }
		s=1;
		while((bottom<height) && (s>0))
		  {
		    s=0;
		    s=getPixel(orig,bgColor,x2,bottom,THRESHOLD);
		    if (s>0) bottom++;
		  }
		s=1;
		while((left>0) && (s>0))
		  {
		    s=0;
		    s=getPixel(orig,bgColor,left,y1,THRESHOLD);
		    if (s>0) left--;
		  }
		s=1;
		while((right<width) && (s>0))
		  {
		    s=0;
		    s=getPixel(orig,bgColor,right,y2,THRESHOLD);
		    if (s>0) right++;
		  }
	      }

	    bool found=false;
	    if (((bottom-top)<=max_font_height) && 
		((right-left)<=max_font_width) && (right-left>V_DISPLACEMENT) 
		&& (bottom-top>MIN_FONT_HEIGHT))
	    {
	      
	      char label=0;
	      label=get_atom_label(orig,bgColor,left,top,right,bottom,THRESHOLD);
	      if (label !=0)
		{
		    //cout<<label<<endl;
		    letters[n_letters].a=label;
		    letters[n_letters].x=(left+right)/2;
		    letters[n_letters].y=(top+bottom)/2;
		    letters[n_letters].r=distance(left,top,right,bottom)/2;
		    if (right-left>real_font_width)
		      real_font_width=right-left;
		    if (bottom-top>real_font_height)
		      real_font_height=bottom-top;
		    letters[n_letters].free=true;
		    n_letters++;
		    if (n_letters>=MAX_ATOMS) n_letters--;
		    delete_bonds_in_char(bond,n_bond,atom,left,top,right,bottom);
		    delete_curve_with_children(atom,bond,n_atom,n_bond,p);
		    found=true;
		  }
	    }
	    if (((bottom-top)<=2*max_font_height) && 
		((right-left)<=max_font_width) && (right-left>V_DISPLACEMENT) 
		&& (bottom-top>MIN_FONT_HEIGHT) && !found)
	      {
	      
	      char label1=0;
	      int newtop=(top+bottom)/2;
	      label1=get_atom_label(orig,bgColor,left,newtop,right,bottom,THRESHOLD);
	      char label2=0;
	      int newbottom=(top+bottom)/2;
	      label2=get_atom_label(orig,bgColor,left,top,right,newbottom,THRESHOLD);
	      if ((label1 !=0) && (label2 != 0))
		  {
		    //cout<<label1<<label2<<endl;
		    letters[n_letters].a=label1;
		    letters[n_letters].x=(left+right)/2;
		    letters[n_letters].y=(newtop+bottom)/2;
		    letters[n_letters].r=distance(left,newtop,right,bottom)/2;
		    if (right-left>real_font_width)
		      real_font_width=right-left;
		    if (bottom-newtop>real_font_height)
		      real_font_height=bottom-newtop;
		    letters[n_letters].free=true;
		    n_letters++;
		    if (n_letters>=MAX_ATOMS) n_letters--;
		    letters[n_letters].a=label2;
		    letters[n_letters].x=(left+right)/2;
		    letters[n_letters].y=(top+newbottom)/2;
		    letters[n_letters].r=distance(left,top,right,newbottom)/2;
		    if (newbottom-top>real_font_height)
		      real_font_height=newbottom-top;
		    letters[n_letters].free=true;
		    n_letters++;
		    if (n_letters>=MAX_ATOMS) n_letters--;
		    delete_bonds_in_char(bond,n_bond,atom,left,top,right,bottom);
		    delete_curve_with_children(atom,bond,n_atom,n_bond,p);
		    found=true;
		  }
	    }
	    if (((bottom-top)<=max_font_height) && 
		((right-left)<=2*max_font_width) && (right-left>V_DISPLACEMENT) 
		&& (bottom-top>MIN_FONT_HEIGHT) && !found)
	    {
	      
	      char label1=0;
	      int newright=(left+right)/2;
	      label1=get_atom_label(orig,bgColor,left,top,newright,bottom,THRESHOLD);
	      char label2=0;
	      int newleft=(left+right)/2;
	      label2=get_atom_label(orig,bgColor,newleft,top,right,bottom,THRESHOLD);
	      if ((label1 !=0) && (label2 != 0))
		  {
		    //cout<<label1<<label2<<endl;
		    letters[n_letters].a=label1;
		    letters[n_letters].x=(left+newright)/2;
		    letters[n_letters].y=(top+bottom)/2;
		    letters[n_letters].r=distance(left,top,newright,bottom)/2;
		    if (newright-left>real_font_width)
		      real_font_width=newright-left;
		    if (bottom-top>real_font_height)
		      real_font_height=bottom-top;
		    letters[n_letters].free=true;
		    n_letters++;
		    if (n_letters>=MAX_ATOMS) n_letters--;
		    letters[n_letters].a=label2;
		    letters[n_letters].x=(newleft+right)/2;
		    letters[n_letters].y=(top+bottom)/2;
		    letters[n_letters].r=distance(newleft,top,right,bottom)/2;
		    if (right-newleft>real_font_width)
		      real_font_width=right-newleft;
		    letters[n_letters].free=true;
		    n_letters++;
		    if (n_letters>=MAX_ATOMS) n_letters--;
		    delete_bonds_in_char(bond,n_bond,atom,left,top,right,bottom);
		    delete_curve_with_children(atom,bond,n_atom,n_bond,p);
		    found=true;
		  }
	    }

	  }
	p = p->next;
      }
  if (real_font_width<1) real_font_width=max_font_width;
  else real_font_width++;
  if (real_font_height<1) real_font_height=max_font_height;
  else real_font_height++;
  return(n_letters);
}

int find_atoms(potrace_path_t *p, atom_t *atom,bond_t *bond,int *n_bond)
{
  int *tag,n_atom=0;
  potrace_dpoint_t (*c)[3];
  long n;

 while (p != NULL) 
      {
	    n = p->curve.n;
	    tag = p->curve.tag;
	    c = p->curve.c;
	    int b_atom=n_atom;
	    atom[n_atom].x=c[n-1][2].x;
	    atom[n_atom].y=c[n-1][2].y;
	    atom[n_atom].label=" ";
	    atom[n_atom].exists=false;
	    atom[n_atom].curve=p;
	    atom[n_atom].n=0;
	    atom[n_atom].corner=false;
	    n_atom++;
	    if (n_atom>=MAX_ATOMS) n_atom--;
	    for (long i=0; i<n; i++) 
	      {
		switch (tag[i]) 
		  {
		  case POTRACE_CORNER:
		    atom[n_atom].x=c[i][1].x;
		    atom[n_atom].y=c[i][1].y;
		    atom[n_atom].label=" ";
		    atom[n_atom].exists=false;
		    atom[n_atom].curve=p;
		    atom[n_atom].n=0;
		    atom[n_atom].corner=true;
		    n_atom++;
		    if (n_atom>=MAX_ATOMS) n_atom--;
		    break;
		  case POTRACE_CURVETO:
		    atom[n_atom].x=c[i][0].x;
		    atom[n_atom].y=c[i][0].y;
		    atom[n_atom].label=" ";
		    atom[n_atom].exists=false;
		    atom[n_atom].curve=p;
		    atom[n_atom].n=0;
		    atom[n_atom].corner=false;
		    n_atom++;
		    if (n_atom>=MAX_ATOMS) n_atom--;
		    atom[n_atom].x=c[i][1].x;
		    atom[n_atom].y=c[i][1].y;
		    atom[n_atom].label=" ";
		    atom[n_atom].exists=false;
		    atom[n_atom].curve=p;
		    atom[n_atom].n=0;
		    atom[n_atom].corner=false;
		    n_atom++;
		    if (n_atom>=MAX_ATOMS) n_atom--;
		    break;
		  }
		if (i!=n-1)
		  {
		    atom[n_atom].x=c[i][2].x;
		    atom[n_atom].y=c[i][2].y;
		    atom[n_atom].label=" ";
		    atom[n_atom].exists=false;
		    atom[n_atom].curve=p;
		    atom[n_atom].n=0;
		    atom[n_atom].corner=false;
		    n_atom++;
		    if (n_atom>=MAX_ATOMS) n_atom--;
		  }
	      }
	    *n_bond=find_bonds(atom,bond,b_atom,n_atom,*n_bond,p);
	    p = p->next;
      }
 return(n_atom);
}


bool check_boxes(int left,int top,int right,int bottom,box_t *boxes,int n_boxes)
{
  for (int i=0;i<n_boxes;i++)
    {
      if ((left>boxes[i].x1) && (left<boxes[i].x2) && (top>boxes[i].y1) && (top<boxes[i].y2)) 
	return(true);
      else if ((right>boxes[i].x1) && (right<boxes[i].x2) && (top>boxes[i].y1) && (top<boxes[i].y2)) 
	return(true);
      else if ((left>boxes[i].x1) && (left<boxes[i].x2) && (bottom>boxes[i].y1) && (bottom<boxes[i].y2)) 
	return(true);
      else if ((right>boxes[i].x1) && (right<boxes[i].x2) && (bottom>boxes[i].y1) && (bottom<boxes[i].y2)) 
	return(true);
    }
  return(false);
}

int calculate_area(potrace_path_t *p)
{
  int area=p->area;
  potrace_path_t *child=p->childlist;
  while (child !=NULL)
    {
      area-=child->area;
      child=child->sibling;
    }
  return(area);
}



int find_boxes(box_t *boxes,Image image,double THRESHOLD_BOND,ColorGray bgColor,
	       int width,int height,int res,int boundary, int working_resolution)
{
  potrace_bitmap_t *bm;
  potrace_param_t *param;
  potrace_path_t *p;
  potrace_state_t *st;
  int n_boxes=0;

  param = potrace_param_default();
  param->alphamax=0.;
  //param->turnpolicy=POTRACE_TURNPOLICY_MINORITY;
  bm = bm_new(width,height);
  param->turdsize=res;
  //param->turdsize=1;

    for(int i=0;i<width;i++)
      for(int j=0;j<height;j++)
	BM_PUT(bm,i,j,getPixel(image,bgColor,i,j,THRESHOLD_BOND));
    
    st = potrace_trace(param, bm);
    p = st->plist;
    while (p != NULL) 
      {
	int top=0;
	int left=0;
	int bottom=0;
	int right=0;
	if ((p->sign == int('+')))
	  {
	    long n = p->curve.n;
	    int *tag = p->curve.tag;
	    potrace_dpoint_t (*c)[3];
	    c = p->curve.c;
	    top=height;
	    left=width;
	    bottom=0;
	    right=0;
	    for (int i=0; i<n; i++) 
	      {
		switch (tag[i]) 
		  {
		  case POTRACE_CORNER:
		    if (c[i][1].x<left) {left=int(c[i][1].x);}
		    if (c[i][1].x>right) {right=int(c[i][1].x);}
		    if (c[i][1].y<top) {top=int(c[i][1].y);}
		    if (c[i][1].y>bottom) {bottom=int(c[i][1].y);}
		    break;
		  case POTRACE_CURVETO:
		    if (c[i][0].x<left) {left=int(c[i][0].x);}
		    if (c[i][0].x>right) {right=int(c[i][0].x);}
		    if (c[i][0].y<top) {top=int(c[i][0].y);}
		    if (c[i][0].y>bottom) {bottom=int(c[i][0].y);}
		    if (c[i][1].x<left) {left=int(c[i][1].x);}
		    if (c[i][1].x>right) {right=int(c[i][1].x);}
		    if (c[i][1].y<top) {top=int(c[i][1].y);}
		    if (c[i][1].y>bottom) {bottom=int(c[i][1].y);}
		    break;
		  }
		if (c[i][2].x<left) {left=int(c[i][2].x);}
		if (c[i][2].x>right) {right=int(c[i][2].x);}
		if (c[i][2].y<top) {top=int(c[i][2].y);}
		if (c[i][2].y>bottom) {bottom=int(c[i][2].y);}
		
		if (left<0) left=0;
		if (top<0) top=0;
		if (right>width-1) right=width-1;
		if (bottom>height-1) bottom=height-1;
	      }
	    int area=calculate_area(p);
	    double ratio=0,aspect=0;
	    if ((bottom!=top) && (right!=left)) ratio=1.*area/((bottom-top)*(right-left));
	    if (right!=left)  aspect=1.*(bottom-top)/(right-left);

	    if ((ratio<MAX_RATIO) && (ratio>0) && (aspect>MIN_ASPECT) && 
		(aspect<MAX_ASPECT) &&
		(!check_boxes(left,top,right,bottom,boxes,n_boxes)))
	      {
		adjust_box(image,THRESHOLD_BOND,bgColor,width,height,boundary, 
			   0,top,left,bottom,right,height,width);
		//left=0;top=0;right=width;bottom=height;
		if (left<0) left=0;
		if (top<0) top=0;
		if (right>width) right=width;
		if (bottom>height) bottom=height;
		if ((right-left)*300/working_resolution<MAX_WIDTH 
		    && (bottom-top)*300/working_resolution<MAX_HEIGHT 
		    && (right-left)>MIN_WIDTH && (bottom-top)>MIN_HEIGHT
		    || working_resolution<150)
		  {
		    boxes[n_boxes].x1=left;
		    boxes[n_boxes].x2=right;
		    boxes[n_boxes].y1=top;
		    boxes[n_boxes].y2=bottom;
		    n_boxes++;
		    if (n_boxes>=NUM_BOXES) n_boxes--;
		  }
	      }
	  }

	p = p->next;
      }
    //draw_box(image,boxes,n_boxes,"tmp.gif");
    potrace_state_free(st);
    potrace_param_free(param);
    free(bm);
    return(n_boxes);
}

int count_pages(string input)
{
  list<Image> imageList;
  readImages( &imageList,input);
  return(imageList.size());
}

string  image_type(string input)
{
  Image image;
  image.ping(input);
  return(image.magick());
}

int count_atoms(atom_t *atom,int n_atom)
{
  int r=0;
  for (int i=0;i<n_atom;i++)
    if (atom[i].exists) r++;
  return(r);
}

int count_bonds(bond_t *bond,int n_bond)
{
  int r=0;
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists) r++;
  return(r);
}

/*------------------- ThinImage - Thin binary image. --------------------------- * 
 *                                                            
 *    Description:                                                    
 *        Thins the supplied binary image using Rosenfeld's parallel   
 *        thinning algorithm.                                         
 *                                                                     
 *    On Entry:                                                        
 *        image = Image to thin.                                       
 *                                                                     
 *------------------------------------------------------------------------------- */ 
 
 
/* Direction masks:                  */ 
/*   N     S     W        E            */ 
static        unsigned int     masks[]         = { 0200, 0002, 0040, 0010 }; 
 
/*    True if pixel neighbor map indicates the pixel is 8-simple and  */ 
/*    not an end point and thus can be deleted.  The neighborhood     */ 
/*    map is defined as an integer of bits abcdefghi with a non-zero  */ 
/*    bit representing a non-zero pixel.  The bit assignment for the  */ 
/*    neighborhood is:                                                */ 
/*                                                                    */ 
/*                            a b c                                   */ 
/*                            d e f                                   */ 
/*                            g h i                                   */ 
 
static        unsigned char   todelete[512] = { 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; 


void thin1(unsigned char *ptr, unsigned int xsize,unsigned int ysize) 
{ 
      unsigned char *y_ptr, *y1_ptr;
      unsigned char bg_color=0,colour=1;
      unsigned int    x, y;           /* Pixel location               */ 
      unsigned int    i;              /* Pass index           */ 
      unsigned int    pc      = 0;    /* Pass count           */ 
      unsigned int    count   = 1;    /* Deleted pixel count          */ 
      unsigned int    p, q;           /* Neighborhood maps of adjacent*/ 
                                      /* cells                        */ 
      unsigned char   *qb;            /* Neighborhood maps of previous*/ 
                                      /* scanline                     */ 
      unsigned int    m;              /* Deletion direction mask      */ 

      qb=(unsigned char*)malloc(xsize*sizeof(unsigned char)); 
      qb[xsize-1] = 0;                /* Used for lower-right pixel   */ 
 
      while ( count ) {               /* Scan image while deletions   */ 
          pc++; 
          count = 0; 
 
          for ( i = 0 ; i < 4 ; i++ ) { 
 
              m = masks[i]; 
 
              /* Build initial previous scan buffer.                  */ 
              p = (ptr[0] == colour); 
              for ( x = 0 ; x < xsize-1 ; x++ ) 
                  qb[x] = (unsigned char) (p = ((p<<1)&0006) | (unsigned int)(ptr[x+1] == colour)); 
 
              /* Scan image for pixel deletion candidates.            */ 
	      y_ptr = ptr; y1_ptr = ptr + xsize; 
              for (y = 0; y < ysize - 1; y++, y_ptr += xsize, y1_ptr += xsize)
	      { 
                  q = qb[0]; 
                  p = ((q<<2)&0330) | (y1_ptr[0] == colour); 
 
                  for ( x = 0 ; x < xsize-1 ; x++ ) { 
                      q = qb[x]; 
                      p = ((p<<1)&0666) | ((q<<3)&0110) | (unsigned int) (y1_ptr[x+1]==colour); 
                      qb[x] = (unsigned char) p; 
                      if  ( ((p&m) == 0) && todelete[p] ) { 
                          count++; 
			  y_ptr[x] = bg_color;  /* delete the pixel */ 
                      } 
                  } 
 
                  /* Process right edge pixel.                        */ 
                  p = (p<<1)&0666; 
                  if  ( (p&m) == 0 && todelete[p] ) { 
                      count++; 
                      y_ptr[xsize-1] = bg_color;
                  } 
              } 
 
              /* Process bottom scan line.                            */ 
	      q = qb[0]; 
	      p = ((q<<2)&0330); 
 
	      y_ptr = ptr + xsize * (ysize - 1);
              for ( x = 0 ; x < xsize ; x++ ) { 
                  q = qb[x]; 
                  p = ((p<<1)&0666) | ((q<<3)&0110); 
                  if  ( (p&m) == 0 && todelete[p] ) { 
                      count++; 
                      y_ptr[x] = bg_color;
                  } 
              } 
          } 
      } 
      free (qb); 
} 

Image thin_image(Image box,double THRESHOLD_BOND,ColorGray bgColor)
{
  Image image( Geometry(box.columns(),box.rows()), "white" );
  image.type( GrayscaleType );
  unsigned int xsize=box.columns();
  unsigned int ysize=box.rows();
  unsigned char *ptr=(unsigned char*)malloc(xsize*ysize*sizeof(unsigned char));
  for (unsigned int i=0; i<xsize;i++)
    for (unsigned int j=0; j<ysize;j++) 
	ptr[i+j*xsize]=getPixel(box,bgColor,i,j,THRESHOLD_BOND);
  thin1(ptr,xsize,ysize);
  for (unsigned int i=0; i<xsize;i++)
    for (unsigned int j=0; j<ysize;j++) 
      if (ptr[i+j*xsize]==1)
	image.pixelColor(i, j, "black" );
  free(ptr);
  return(image);
}


int comp_dashes_x(const void *a,const void *b)
{
  dash_t *aa=(dash_t *) a;
  dash_t *bb=(dash_t *) b;
  if (aa->x<bb->x) return(-1);
  if (aa->x==bb->x) return(0);
  if (aa->x>bb->x) return(1);
  return(0);
}

int comp_dashes_y(const void *a,const void *b)
{
  dash_t *aa=(dash_t *) a;
  dash_t *bb=(dash_t *) b;
  if (aa->y<bb->y) return(-1);
  if (aa->y==bb->y) return(0);
  if (aa->y>bb->y) return(1);
  return(0);
}


void extend_dashed_bond(int a,int b,int n,atom_t *atom,double avg)
{
  double x0=atom[a].x;
  double y0=atom[a].y;
  double x1=atom[b].x;
  double y1=atom[b].y;
  double l=distance(x0,y0,x1,y1);
  double kx=(x1-x0)/l;
  double ky=(y1-y0)/l;
  //double e=max(avg,l);
  //atom[b].x=kx*e+x0;
  //atom[b].y=ky*e+y0;
  atom[a].x=kx*(-1.*l/(n-1))+x0;
  atom[a].y=ky*(-1.*l/(n-1))+y0;
  atom[b].x=kx*l/(n-1)+x1;
  atom[b].y=ky*l/(n-1)+y1;
}

int count_area(vector < vector<int> > *box, double &x0, double &y0)
{
  int a=0;
  int w=(*box).size();
  int h=(*box)[0].size();
  int x=int(x0);
  int y=int(y0);
  int xm=0,ym=0;
  if ((*box)[x][y]==1)
    {
      (*box)[x][y]=2;
      list<int> cx;
      list<int> cy;
      cx.push_back(x);
      cy.push_back(y);
      while(!cx.empty())
	{
	  x=cx.front();
	  y=cy.front();
	  cx.pop_front();
	  cy.pop_front();
	  (*box)[x][y]=0;
	  a++;
	  xm+=x;
	  ym+=y;
	  for(int i=x-1;i<x+2;i++)
	    for (int j=y-1;j<y+2;j++)
	      if (i<w && j<h && i>=0 && j>=0 && (*box)[i][j]==1)
		{
		  cx.push_back(i);
		  cy.push_back(j);
		  (*box)[i][j]=2;
		}
	}
    }
  else return(0);
  x0=1.*xm/a;
  y0=1.*ym/a;
  return(a);
}

int find_dashed_bonds(potrace_path_t *p, atom_t *atom,bond_t *bond,int n_atom,
		      int *n_bond,int max,double avg,Image img,ColorGray bg,  
		      double THRESHOLD, bool thick)
{
  int n,n_dot=0;
  potrace_dpoint_t (*c)[3];
  dash_t dot[100];
  vector < vector<int> > box(img.columns());

  for (unsigned int i=0;i<img.columns();i++)
    for (unsigned int j=0;j<img.rows();j++)
      box[i].push_back(getPixel(img,bg,i,j,THRESHOLD));

  while (p != NULL) 
      {
	if ((p->sign == int('+')) && (p->area<max))
	  {
	    n = p->curve.n;
	    c = p->curve.c;
	    int *tag = p->curve.tag;
	    dot[n_dot].x=c[n-1][2].x;
	    dot[n_dot].y=c[n-1][2].y;
	    double l=c[n-1][2].x;
	    double r=c[n-1][2].x;
	    double t=c[n-1][2].y;
	    double b=c[n-1][2].y;
	    dot[n_dot].curve=p;
	    dot[n_dot].free=true;
	    int tot=1;
	    for (long i=0; i<n; i++) 
	      {
		switch (tag[i]) 
		  {
		  case POTRACE_CORNER:
		    dot[n_dot].x+=c[i][1].x;
		    dot[n_dot].y+=c[i][1].y;
		    if (c[i][1].x<l) l=c[i][1].x;
		    if (c[i][1].x>r) r=c[i][1].x;
		    if (c[i][1].y<t) t=c[i][1].y;
		    if (c[i][1].x>b) b=c[i][1].y;
		    tot++;
		    break;
		  case POTRACE_CURVETO:
		    dot[n_dot].x+=c[i][0].x;
		    dot[n_dot].y+=c[i][0].y;
		    if (c[i][0].x<l) l=c[i][0].x;
		    if (c[i][0].x>r) r=c[i][0].x;
		    if (c[i][0].y<t) t=c[i][0].y;
		    if (c[i][0].x>b) b=c[i][0].y;
		    dot[n_dot].x+=c[i][1].x;
		    dot[n_dot].y+=c[i][1].y;
		    if (c[i][1].x<l) l=c[i][1].x;
		    if (c[i][1].x>r) r=c[i][1].x;
		    if (c[i][1].y<t) t=c[i][1].y;
		    if (c[i][1].x>b) b=c[i][1].y;
		    tot+=2;
		    break;
		  }
		if (i!=n-1)
		  {
		    dot[n_dot].x+=c[i][2].x;
		    dot[n_dot].y+=c[i][2].y;
		    if (c[i][2].x<l) l=c[i][2].x;
		    if (c[i][2].x>r) r=c[i][2].x;
		    if (c[i][2].y<t) t=c[i][2].y;
		    if (c[i][2].x>b) b=c[i][2].y;
		    tot++;
		  }
	      }
	    dot[n_dot].x/=tot;
	    dot[n_dot].y/=tot;
	    if (thick)
	      dot[n_dot].area=count_area(&box,dot[n_dot].x,dot[n_dot].y);
	    else
	      dot[n_dot].area=p->area;
	    if (distance(l,t,r,b)<avg/3) 
	      n_dot++;
	    if (n_dot>=100) n_dot--;
	  }
	p = p->next;
      }
  for(int i=0;i<n_dot;i++)
    if (dot[i].free)
      {
	dash_t dash[100];
	dash[0]=dot[i];
	dot[i].free=false;
	double l=dot[i].x;
	double r=dot[i].x;
	double t=dot[i].y;
	double b=dot[i].y;
	double mx=l;
	double my=t;
	double dist_next=FLT_MAX;
	int next_dot=i;
	for(int j=i+1;j<n_dot;j++)
	  if (dot[j].free
	      && distance(dash[0].x,dash[0].y,dot[j].x,dot[j].y)<=1.8*avg
	      && distance(dash[0].x,dash[0].y,dot[j].x,dot[j].y)<dist_next)
	    {
	      dash[1]=dot[j];
	      dist_next=distance(dash[0].x,dash[0].y,dot[j].x,dot[j].y);
	      next_dot=j;
	    }

	int n=1;
	if (next_dot!=i)
	  {
	    dot[next_dot].free=false;
	    if (dash[1].x<l) l=dash[1].x;
	    if (dash[1].x>r) r=dash[1].x;
	    if (dash[1].y<t) t=dash[1].y;
	    if (dash[1].y>b) b=dash[1].y;
	    mx=(mx+dash[1].x)/2;
	    my=(my+dash[1].y)/2;
	    n=2;
	  }
	bool found=true;
	while (n>1 && found)
	  {
	    dist_next=FLT_MAX;
	    found=false;
	    int minj=next_dot;
	    for(int j=next_dot+1;j<n_dot;j++)
	      if (dot[j].free && distance(mx,my,dot[j].x,dot[j].y)<=1.8*avg
		  && distance(mx,my,dot[j].x,dot[j].y)<dist_next
		  && fabs(angle4(dash[0].x,dash[0].y,dash[n-1].x,dash[n-1].y,
				 dash[0].x,dash[0].y,dot[j].x,dot[j].y))>D_T_TOLERANCE)
		{
		  dash[n]=dot[j];
		  dist_next=distance(mx,my,dot[j].x,dot[j].y);
		  found=true;
		  minj=j;
		}
	    if (found) 
	      {
		dot[minj].free=false;
		if (dash[n].x<l) l=dash[n].x;
		if (dash[n].x>r) r=dash[n].x;
		if (dash[n].y<t) t=dash[n].y;
		if (dash[n].y>b) b=dash[n].y;
		mx=(mx+dash[n].x)/2;
		my=(my+dash[n].y)/2;
		n++;
	      }
	  }
      
	if (n>2) 
	  {
	    if((r-l)>(b-t))
	      {
		qsort(dash,n,sizeof(dash_t),comp_dashes_x);
	      }
	    else
	      {
		qsort(dash,n,sizeof(dash_t),comp_dashes_y);
	      }
	    bool one_line=true;
	    double dx=dash[n-1].x-dash[0].x;
	    double dy=dash[n-1].y-dash[0].y;
	    double k=0;
	    if (fabs(dx)>fabs(dy)) k=dy/dx;
	    else k=dx/dy;
	    for(int j=1;j<n-1;j++)
	      {
		double nx=dash[j].x-dash[0].x;
		double ny=dash[j].y-dash[0].y;
		double diff=0;
		if (fabs(dx)>fabs(dy))
		  diff=k*nx-ny;
		else
		  diff=k*ny-nx;
		if (fabs(diff)>V_DISPLACEMENT) one_line=false;
	      }
	    if (one_line)
	      {
		for(int j=0;j<n;j++)   
		  delete_curve(atom,bond,n_atom,*n_bond,dash[j].curve);
		atom[n_atom].x=dash[0].x;
		atom[n_atom].y=dash[0].y;
		atom[n_atom].label=" ";
		atom[n_atom].exists=true;
		atom[n_atom].curve=dash[0].curve;
		atom[n_atom].n=0;
		atom[n_atom].corner=false;
		n_atom++;
		if (n_atom>=MAX_ATOMS) n_atom--;
		atom[n_atom].x=dash[n-1].x;
		atom[n_atom].y=dash[n-1].y;
		atom[n_atom].label=" ";
		atom[n_atom].exists=true;
		atom[n_atom].curve=dash[n-1].curve;
		atom[n_atom].n=0;
		atom[n_atom].corner=false;
		n_atom++;
		if (n_atom>=MAX_ATOMS) n_atom--;
		bond[*n_bond].a=n_atom-2;
		bond[*n_bond].exists=true;
		bond[*n_bond].type=1;
		bond[*n_bond].b=n_atom-1;
		bond[*n_bond].curve=dash[0].curve;
		if (dash[0].area>dash[n-1].area)
		  bond_end_swap(bond,*n_bond);
		bond[*n_bond].hash=true;
		bond[*n_bond].wedge=false;
		bond[*n_bond].up=false;
		bond[*n_bond].down=false;
		bond[*n_bond].Small=false;
		extend_dashed_bond(bond[*n_bond].a,bond[*n_bond].b,n,atom,avg);
		(*n_bond)++;
		if ((*n_bond)>=MAX_ATOMS) (*n_bond)--;
	      }
	  }
      }
	     
  return(n_atom);
}

int find_small_bonds(potrace_path_t *p, atom_t *atom,bond_t *bond,int n_atom,
		     int *n_bond,double max_area,double Small,double thickness)
{

  while (p != NULL) 
      {
	if ((p->sign == int('+')) && (p->area<=max_area))
	  {
	    int n_dot=0;
	    dash_t dot[20];
	    for (int i=0;i<n_atom;i++)
	      if ((atom[i].exists) && (atom[i].curve==p) && (n_dot<20))
		{
		  dot[n_dot].x=atom[i].x;
		  dot[n_dot].y=atom[i].y;
		  dot[n_dot].curve=p;
		  dot[n_dot].free=true;
		  n_dot++;
		  if (n_dot>=20) n_dot--;
		}
      
	    if ((n_dot>2))
	      {
		double l=dot[0].x;
		double r=dot[0].x;
		double t=dot[0].y;
		double b=dot[0].y;
		for(int i=1;i<n_dot;i++)
		  {
		    if (dot[i].x<l) l=dot[i].x;
		    if (dot[i].x>r) r=dot[i].x;
		    if (dot[i].y<t) t=dot[i].y;
		    if (dot[i].y>b) b=dot[i].y;
		  }
		if((r-l)>(b-t))
		  {
		    qsort(dot,n_dot,sizeof(dash_t),comp_dashes_x);
		  }
		else
		  {
		    qsort(dot,n_dot,sizeof(dash_t),comp_dashes_y);
		  }
		double d=0;
		for (int i=1;i<n_dot-1;i++)
		  d=max(d,fabs(distance_from_bond_y(dot[0].x,dot[0].y,dot[n_dot-1].x,
						    dot[n_dot-1].y,dot[i].x,dot[i].y)));
		if (d<thickness || p->area<Small)
		  {
		    delete_curve(atom,bond,n_atom,*n_bond,p);
		    atom[n_atom].x=dot[0].x;
		    atom[n_atom].y=dot[0].y;
		    atom[n_atom].label=" ";
		    atom[n_atom].exists=true;
		    atom[n_atom].curve=p;
		    atom[n_atom].n=0;
		    atom[n_atom].corner=false;
		    n_atom++;
		    if (n_atom>=MAX_ATOMS) n_atom--;
		    atom[n_atom].x=dot[n_dot-1].x;
		    atom[n_atom].y=dot[n_dot-1].y;
		    atom[n_atom].label=" ";
		    atom[n_atom].exists=true;
		    atom[n_atom].curve=p;
		    atom[n_atom].n=0;
		    atom[n_atom].corner=false;
		    n_atom++;
		    if (n_atom>=MAX_ATOMS) n_atom--;
		    bond[*n_bond].a=n_atom-2;
		    bond[*n_bond].exists=true;
		    bond[*n_bond].type=1;
		    bond[*n_bond].b=n_atom-1;
		    bond[*n_bond].curve=p;
		    bond[*n_bond].hash=false;
		    bond[*n_bond].wedge=false;
		    bond[*n_bond].up=false;
		    bond[*n_bond].down=false;
		    bond[*n_bond].Small=true;
		    (*n_bond)++;
		    if ((*n_bond)>=MAX_ATOMS) (*n_bond)--;
		  }
	      }
	  }
	p = p->next;
      }
  return(n_atom);
}

int resolve_bridge_bonds(atom_t* atom,int n_atom,bond_t* bond,int n_bond)
{
  int rotors1,rotors2,f1,f2,rings1,rings2;
  double confidence;
  string smiles1=get_smiles(atom,bond,n_bond,rotors1,confidence,f1,rings1);
  for (int i=0;i<n_atom;i++)
    if ((atom[i].exists) && (atom[i].label==" "))
      {
	list<int> con;
	for (int j=0;j<n_bond;j++)
	  if ((bond[j].exists) && (bond[j].a==i || bond[j].b==i))
	    con.push_back(j);
	if (con.size()==4)
	  {
	    int a=con.front();
	    con.pop_front();
	    int b=0;
	    int e=0;
	    while ((con.size()>2) && (e++<3))
	      {
		b=con.front();
		con.pop_front();
		if (angle_between_connected_bonds(bond,a,b,atom)<FLAT_TOLERANCE)
		  con.push_back(b);
	      }
	    if (con.size()==2)
	      {
		int c=con.front();
		con.pop_front();
		int d=con.front();
		con.pop_front();
		vector<int> term;
		term.push_back(a);term.push_back(b);term.push_back(c);term.push_back(d);
		bool terminal=false;
		for (unsigned int k=0;k<term.size();k++)
		  {
		    bool terminal_a=terminal_bond(bond[term[k]].a,term[k],bond,n_bond);
		    bool terminal_b=terminal_bond(bond[term[k]].b,term[k],bond,n_bond);
		    if (terminal_a || terminal_b) terminal=true;
		  }
		if (bond[a].type==1 && bond[b].type==1 &&
		    bond[c].type==1 && bond[d].type==1 &&
		    angle_between_connected_bonds(bond,c,d,atom)>FLAT_TOLERANCE
		    && !terminal
		    )
		  {
		    bond[b].exists=false;
		    bond[d].exists=false;
		    atom[i].exists=false;
		    if (bond[a].a==bond[b].a) bond[a].a=bond[b].b;
		    else if (bond[a].a==bond[b].b) bond[a].a=bond[b].a;
		    else if (bond[a].b==bond[b].a) bond[a].b=bond[b].b;
		    else if (bond[a].b==bond[b].b) bond[a].b=bond[b].a;
		    if (bond[c].a==bond[d].a) bond[c].a=bond[d].b;
		    else if (bond[c].a==bond[d].b) bond[c].a=bond[d].a;
		    else if (bond[c].b==bond[d].a) bond[c].b=bond[d].b;
		    else if (bond[c].b==bond[d].b) bond[c].b=bond[d].a;
		    string smiles2=get_smiles(atom,bond,n_bond,rotors2,
					      confidence,f2,rings2);
		    if (f1!=f2 || rotors1!=rotors2 || rings1-rings2==2)
		      {
			bond[b].exists=true;
			bond[d].exists=true;
			atom[i].exists=true;
			if (bond[a].a==bond[b].a) bond[a].a=bond[b].b;
			else if (bond[a].a==bond[b].b) bond[a].a=bond[b].a;
			else if (bond[a].b==bond[b].a) bond[a].b=bond[b].b;
			else if (bond[a].b==bond[b].b) bond[a].b=bond[b].a;
			if (bond[c].a==bond[d].a) bond[c].a=bond[d].b;
			else if (bond[c].a==bond[d].b) bond[c].a=bond[d].a;
			else if (bond[c].b==bond[d].a) bond[c].b=bond[d].b;
			else if (bond[c].b==bond[d].b) bond[c].b=bond[d].a;
		      }
		  }
	      }
	  }
      }
  return(f1);
}

void collapse_atoms(atom_t *atom, bond_t *bond, int n_atom,int n_bond,
		    double dist)
{
  bool found=true;
  while (found)
    {
      found=false;
      for (int i=0;i<n_atom;i++)
	if (atom[i].exists)
	  for (int j=0;j<n_atom;j++)
	    if (atom[j].exists && j!=i &&
		distance(atom[i].x,atom[i].y,atom[j].x,atom[j].y)<dist)
	      {
		atom[j].exists=false;
		atom[i].x=(atom[i].x+atom[j].x)/2;
		atom[i].y=(atom[i].y+atom[j].y)/2;
		if (atom[j].label!=" " && atom[i].label==" ")
		  atom[i].label=atom[j].label;
		for (int k=0;k<n_bond;k++)
		  if (bond[k].exists)
		    if (bond[k].a==j) {bond[k].a=i;}
		    else if (bond[k].b==j) {bond[k].b=i;}
		found=true;
	      }
    }
}


int fix_one_sided_bonds(bond_t *bond,int n_bond,atom_t *atom, double thickness)
{
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists && bond[i].type==1)
      for (int j=0;j<n_bond;j++)
	if (bond[j].exists && j!=i && bond[j].type==1)
	  {
	    double d1=fabs(distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,
						atom[bond[i].b].x,atom[bond[i].b].y,
						atom[bond[j].a].x,atom[bond[j].a].y));
	    double d2=fabs(distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,
						atom[bond[i].b].x,atom[bond[i].b].y,
						atom[bond[j].b].x,atom[bond[j].b].y));
	    double l=bond_length(bond,i,atom);
	    if (d1<thickness && !(bond[j].a==bond[i].b || bond[j].a==bond[i].a))
	      {
		double l1=distance(atom[bond[j].a].x,atom[bond[j].a].y,
				   atom[bond[i].a].x,atom[bond[i].a].y);
		double l2=distance(atom[bond[j].a].x,atom[bond[j].a].y,
				   atom[bond[i].b].x,atom[bond[i].b].y);
		if (l1<l && l2<l)
		  {
		    if (bond[j].b==bond[i].b || bond[j].b==bond[i].a)
		      {
			bond[j].exists=false;
		      }
		    else
		      {
			bond[n_bond].b=bond[i].b;
			bond[n_bond].exists=true;
			bond[n_bond].type=1;
			bond[n_bond].a=bond[j].a;
			bond[n_bond].curve=bond[i].curve;
			bond[n_bond].hash=false;
			if (bond[i].wedge) bond[n_bond].wedge=true;
			else bond[n_bond].wedge=false;
			bond[n_bond].Small=false;
			bond[n_bond].up=false;
			bond[n_bond].down=false;
			if (bond[i].arom) bond[n_bond].arom=true;
			else bond[n_bond].arom=false;
			n_bond++;
			if (n_bond>=MAX_ATOMS) n_bond--;
			bond[i].b=bond[j].a;
			bond[i].wedge=false;
		      }
		  }
	      }
	    else  if (d2<thickness && !(bond[j].b==bond[i].b || bond[j].b==bond[i].a))
	      {
		double l1=distance(atom[bond[j].b].x,atom[bond[j].b].y,
				   atom[bond[i].a].x,atom[bond[i].a].y);
		double l2=distance(atom[bond[j].b].x,atom[bond[j].b].y,
				   atom[bond[i].b].x,atom[bond[i].b].y);
		if (l1<l && l2<l)
		  {
		    if (bond[j].a==bond[i].b || bond[j].a==bond[i].a)
		      {
			bond[j].exists=false;
		      }
		    else
		      {
			bond[n_bond].b=bond[i].b;
			bond[n_bond].exists=true;
			bond[n_bond].type=1;
			bond[n_bond].a=bond[j].b;
			bond[n_bond].curve=bond[i].curve;
			bond[n_bond].hash=false;
			if (bond[i].wedge) bond[n_bond].wedge=true;
			else bond[n_bond].wedge=false;
			bond[n_bond].Small=false;
			bond[n_bond].up=false;
			bond[n_bond].down=false;
			if (bond[i].arom) bond[n_bond].arom=true;
			else bond[n_bond].arom=false;
			n_bond++;
			if (n_bond>=MAX_ATOMS) n_bond--;
			bond[i].b=bond[j].b;
			bond[i].wedge=false;
		      }
		  }
	      }
	  }
  return(n_bond);
}

int find_fused_chars(bond_t *bond,int n_bond,atom_t *atom,
		     letters_t *letters,int n_letters,
		     int max_font_height,int max_font_width,
		     double r, Image orig,  ColorGray bgColor, double THRESHOLD)
{
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists && bond_length(bond,i,atom)<r)
      {
	list<int> t;
	t.push_back(i);
	for (int j=i+1;j<n_bond;j++)
	  if (bond[j].exists && bond_length(bond,j,atom)<r)
	    {
	      double dx=max(max(fabs(atom[bond[j].a].x-atom[bond[i].a].x),
				fabs(atom[bond[j].a].x-atom[bond[i].b].x)),
			    max(fabs(atom[bond[j].b].x-atom[bond[i].a].x),
				fabs(atom[bond[j].b].x-atom[bond[i].b].x)));
	      double dy=max(max(fabs(atom[bond[j].a].y-atom[bond[i].a].y),
				fabs(atom[bond[j].a].y-atom[bond[i].b].y)),
			    max(fabs(atom[bond[j].b].y-atom[bond[i].a].y),
				fabs(atom[bond[j].b].y-atom[bond[i].b].y)));
	      if (dx<max_font_width && dy<max_font_height)
		t.push_back(j);
	    }
	if (t.size()>2)
	  {
	    double cx=0;
	    double cy=0;
	    int n=0;
	    while (!t.empty())
	      {
		int k=t.front();
		t.pop_front();
		cx+=atom[bond[k].a].x+atom[bond[k].b].x;
		cy+=atom[bond[k].a].y+atom[bond[k].b].y;
		n+=2;
	      }
	    cx/=n;
	    cy/=n;
	    int left=int(cx-max_font_width/2);
	    int right=int(cx+max_font_width/2);
	    int top=int(cy-max_font_height/2);
	    int bottom=int(cy+max_font_height/2);
	    

	    char label=0;
	    label=get_atom_label(orig,bgColor,left,top,right,bottom,THRESHOLD);
	    if (label !=0 
		&& label!='P' && label!='p' && label!='F' 
		&& label!='X' && label!='Y'
		&& label!='n' && label!='F' && label!='U' && label!='u'
		&& label!='h'
		)
	      {
		bool overlap=false;
		for (int j=0;j<n_letters;j++)
		  {
		    if (distance((left+right)/2,(top+bottom)/2,
				 letters[j].x,letters[j].y)<letters[j].r)
		      overlap=true;
		  }
		if (!overlap)
		  {
		    letters[n_letters].a=label;
		    letters[n_letters].x=(left+right)/2;
		    letters[n_letters].y=(top+bottom)/2;
		    letters[n_letters].r=distance(left,top,right,bottom)/2;
		    letters[n_letters].free=true;
		    n_letters++;
		    if (n_letters>=MAX_ATOMS) n_letters--;
		    delete_bonds_in_char(bond,n_bond,atom,left,top,right,bottom);
		  }

	      }
	  }
      }
  return(n_letters);
}

int comp_boxes(const void *a,const void *b)
{
  box_t *aa=(box_t *) a;
  box_t *bb=(box_t *) b;
  if (aa->y2<bb->y1) return(-1);
  if (aa->y1>bb->y2) return(1);
  if (aa->x1>bb->x1) return(1);
  if (aa->x1<bb->x1) return(-1);
  return(0);
}


double noise_factor(Image image, int width, int height, ColorGray bgColor, 
		     double THRESHOLD_BOND)
{
  int n1=0,n2=0,n3=0;
  for(int i=0;i<width;i++)
    {
      int j=0;
      while(j<height)
	{
	  while(!getPixel(image,bgColor,i,j,THRESHOLD_BOND) && j<height) j++;
	  int l=0;
	  while(getPixel(image,bgColor,i,j,THRESHOLD_BOND) && j<height) 
	    {
	      l++;
	      j++;
	    }
	  if (l==1) n1++;
	  else if (l==2) n2++;
	  else if (l==3) n3++;
	}
    }
  for(int i=0;i<height;i++)
    {
      int j=0;
      while(j<width)
	{
	  while(!getPixel(image,bgColor,j,i,THRESHOLD_BOND) && j<width) j++;
	  int l=0;
	  while(getPixel(image,bgColor,j,i,THRESHOLD_BOND) && j<width) 
	    {
	      l++;
	      j++;
	    }
	  if (l==1) n1++;
	  else if (l==2) n2++;
	  else if (l==3) n3++;
	}
    }
  //  cout<<n1<<" "<<n2<<" "<<n3<<endl;
  return(1.*n2/n3);
}

double thickness_hor(Image image,int x1,int y1, ColorGray bgColor, 
		     double THRESHOLD_BOND)
{
  int i=0,s=0,w=0;
  int width=image.columns();
  s=getPixel(image,bgColor,x1,y1,THRESHOLD_BOND);
  if (s==0 && x1+1<width)
    {
      x1++;
      s=getPixel(image,bgColor,x1,y1,THRESHOLD_BOND);
    }
  if (s==0 && x1-2>=0)
    {
      x1-=2;
      s=getPixel(image,bgColor,x1,y1,THRESHOLD_BOND);
    }
  if (s==1)
    {
      while(x1+i<width && s==1)
	s=getPixel(image,bgColor,x1+i++,y1,THRESHOLD_BOND);
      w=i-1;
      i=1;s=1;
      while(x1-i>=0 && s==1)
	s=getPixel(image,bgColor,x1-i++,y1,THRESHOLD_BOND);
      w+=i-1;
    }
  return(1.*w);
}

double thickness_ver(Image image,int x1,int y1, ColorGray bgColor, 
		     double THRESHOLD_BOND)
{
  int i=0,s=0,w=0;
  int height=image.rows();
  s=getPixel(image,bgColor,x1,y1,THRESHOLD_BOND);
  if (s==0 && y1+1<height)
    {
      y1++;
      s=getPixel(image,bgColor,x1,y1,THRESHOLD_BOND);
    }
  if (s==0 && y1-2>=0)
    {
      y1-=2;
      s=getPixel(image,bgColor,x1,y1,THRESHOLD_BOND);
    }
  if (s==1)
    {
      while(y1+i<height && s==1)
	s=getPixel(image,bgColor,x1,y1+i++,THRESHOLD_BOND);
      w=i-1;
      i=1;s=1;
      while(y1-i>=0 && s==1)
	s=getPixel(image,bgColor,x1,y1-i++,THRESHOLD_BOND);
      w+=i-1;
    }
  return(1.*w);
}



double find_wedge_bonds(Image image,atom_t* atom, int n_atom,bond_t* bond,int n_bond, 
			ColorGray bgColor,double THRESHOLD_BOND, 
			double max_dist_double_bond, double avg)
{
  double l,t=0;
  double a[MAX_ATOMS];
  int n=0;

  for (int i=0;i<n_bond;i++)
    if (bond[i].exists && !bond[i].hash && bond[i].type==1 
	&& (l=bond_length(bond,i,atom))>max_dist_double_bond)
      {
	double d=l/16;
	double x1=atom[bond[i].a].x+(atom[bond[i].b].x-atom[bond[i].a].x)*d/l;
	double y1=atom[bond[i].a].y+(atom[bond[i].b].y-atom[bond[i].a].y)*d/l;
	double x2=atom[bond[i].b].x+(atom[bond[i].a].x-atom[bond[i].b].x)*d/l;
	double y2=atom[bond[i].b].y+(atom[bond[i].a].y-atom[bond[i].b].y)*d/l;
	
	double w1,w2,w3;
	double w3_ver=thickness_ver(image,int((x1+x2)/2),int((y1+y2)/2),bgColor,
				    THRESHOLD_BOND);
	double w3_hor=thickness_hor(image,int((x1+x2)/2),int((y1+y2)/2),bgColor,
				    THRESHOLD_BOND);
	if (w3_ver==0 && w3_hor==0) continue;
	if ((w3_ver<w3_hor && w3_ver>0) || w3_hor==0)
	  {
	    w1=thickness_ver(image,int(x1),int(y1),bgColor,THRESHOLD_BOND);
	    w2=thickness_ver(image,int(x2),int(y2),bgColor,THRESHOLD_BOND);
	    w3=w3_ver;
	  }
	else
	  {
	    w1=thickness_hor(image,int(x1),int(y1),bgColor,THRESHOLD_BOND);
	    w2=thickness_hor(image,int(x2),int(y2),bgColor,THRESHOLD_BOND);
	    w3=w3_hor;
	  }
	if (w2>w3 && w3>w1 && w2<2*MAX_BOND_THICKNESS && w2<avg/3 && w1>0)
	  {
	    bond[i].wedge=true;
	  }
	if (w1>w3 && w3>w2 && w1<2*MAX_BOND_THICKNESS && w1<avg/3 && w2>0)
	  {
	    bond[i].wedge=true;
	    bond_end_swap(bond,i);
	  }
	if (bond[i].wedge)
	  {
	    double w=max(w1,w2);
	    if (w<2*MAX_BOND_THICKNESS && w<avg/3) 
	      for (int j=0;j<n_atom;j++)
		if (atom[j].exists && j!=bond[i].b &&
		    distance(atom[bond[i].b].x,atom[bond[i].b].y,
			     atom[j].x,atom[j].y)<=w)
		  {
		    atom[j].exists=false;
		    atom[bond[i].b].x=(atom[bond[i].b].x+atom[j].x)/2;
		    atom[bond[i].b].y=(atom[bond[i].b].y+atom[j].y)/2;
		    for (int k=0;k<n_bond;k++)
		      if (bond[k].exists)
			if (bond[k].a==j) {bond[k].a=bond[i].b;}
			else if (bond[k].b==j) {bond[k].b=bond[i].b;}
		  }
	  }
	if (!bond[i].wedge) a[n++]=w3;
      }
  qsort(a,n,sizeof(double),num_comp);
  t=a[n/2];
  return(t);
}

void collapse_double_bonds(bond_t *bond,int n_bond,atom_t *atom,int n_atom,double dist)
{
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists && bond[i].type==2 && bond[i].conjoined)
      for (int j=0;j<n_bond;j++)
	if (bond[j].exists && j!=i && bond[j].type==1 && bond_length(bond,j,atom)<=dist)
	  if (bond[j].a==bond[i].a)
	    {
	      bond[j].exists=false;
	      atom[bond[i].a].x=(atom[bond[i].a].x+atom[bond[j].b].x)/2;
	      atom[bond[i].a].y=(atom[bond[i].a].y+atom[bond[j].b].y)/2;
	      for (int k=0;k<n_bond;k++)
		if (bond[k].exists)
		  if (bond[k].a==bond[j].b) {bond[k].a=bond[i].a;}
		  else if (bond[k].b==bond[j].b) {bond[k].b=bond[i].a;}
	    }
	  else if (bond[j].b==bond[i].a)
	    {
	      bond[j].exists=false;
	      atom[bond[i].a].x=(atom[bond[i].a].x+atom[bond[j].a].x)/2;
	      atom[bond[i].a].y=(atom[bond[i].a].y+atom[bond[j].a].y)/2;
	      for (int k=0;k<n_bond;k++)
		if (bond[k].exists)
		  if (bond[k].a==bond[j].a) {bond[k].a=bond[i].a;}
		  else if (bond[k].b==bond[j].a) {bond[k].b=bond[i].a;}
	    }
	  else if (bond[j].a==bond[i].b)
	    {
	      bond[j].exists=false;
	      atom[bond[i].b].x=(atom[bond[i].b].x+atom[bond[j].b].x)/2;
	      atom[bond[i].b].y=(atom[bond[i].b].y+atom[bond[j].b].y)/2;
	      for (int k=0;k<n_bond;k++)
		if (bond[k].exists)
		  if (bond[k].a==bond[j].b) {bond[k].a=bond[i].b;}
		  else if (bond[k].b==bond[j].b) {bond[k].b=bond[i].b;}
	    }
	  else if (bond[j].b==bond[i].b)
	    {
	      bond[j].exists=false;
	      atom[bond[i].b].x=(atom[bond[i].b].x+atom[bond[j].a].x)/2;
	      atom[bond[i].b].y=(atom[bond[i].b].y+atom[bond[j].a].y)/2;
	      for (int k=0;k<n_bond;k++)
		if (bond[k].exists)
		  if (bond[k].a==bond[j].a) {bond[k].a=bond[i].b;}
		  else if (bond[k].b==bond[j].a) {bond[k].b=bond[i].b;}
	    }

}

void find_up_down_bonds(bond_t* bond,int n_bond,atom_t* atom, double thickness)
{
  for(int i=0;i<n_bond;i++)
    if(bond[i].exists && bond[i].type==2)
      {
	if (atom[bond[i].a].x>atom[bond[i].b].x)
	  bond_end_swap(bond,i);
	if (atom[bond[i].a].x==atom[bond[i].b].x && atom[bond[i].a].y>atom[bond[i].b].y)
	  bond_end_swap(bond,i);

	for(int j=0;j<n_bond;j++)
	  if (bond[j].exists && bond[j].type==1 && !bond[j].wedge && !bond[j].hash)
	    {
	      bond[j].down=false;
	      bond[j].up=false;
	      if (bond[j].b==bond[i].a)
		{
		  double h=distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,
						atom[bond[i].b].x,atom[bond[i].b].y,
						atom[bond[j].a].x,atom[bond[j].a].y);
		  if (h>thickness)  bond[j].down=true;
		  else if (h<-thickness) bond[j].up=true;
		}
	      else if (bond[j].a==bond[i].a)
		{
		  bond_end_swap(bond,j);
		  double h=distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,
						atom[bond[i].b].x,atom[bond[i].b].y,
						atom[bond[j].a].x,atom[bond[j].a].y);
		  if (h>thickness)  bond[j].down=true;
		  else if (h<-thickness) bond[j].up=true;
		}
	      else if (bond[j].a==bond[i].b)
		{
		  double h=distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,
						atom[bond[i].b].x,atom[bond[i].b].y,
						atom[bond[j].b].x,atom[bond[j].b].y);
		  if (h>thickness)  bond[j].up=true;
		  else if (h<-thickness) bond[j].down=true;
		}
	      else if (bond[j].b==bond[i].b)
		{
		  bond_end_swap(bond,j);
		  double h=distance_from_bond_y(atom[bond[i].a].x,atom[bond[i].a].y,
						atom[bond[i].b].x,atom[bond[i].b].y,
						atom[bond[j].b].x,atom[bond[j].b].y);
		  if (h>thickness)  bond[j].up=true;
		  else if (h<-thickness) bond[j].down=true;
		}
	    }
      }
}

bool detect_curve(bond_t *bond,int n_bond, potrace_path_t *curve)
{
  bool res=false;
  for(int i=0;i<n_bond;i++)
    if (bond[i].exists && bond[i].curve==curve && bond[i].type==1 
	&& !bond[i].wedge && !bond[i].hash) 
      res=true;
  return(res);
}

int find_plus_minus(potrace_path_t *p,letters_t *letters,
		    atom_t *atom,bond_t *bond,int n_atom,int n_bond,int height,
		    int width, int max_font_height, int max_font_width,
		    int n_letters, double avg)
{
  int n, *tag;
  potrace_dpoint_t (*c)[3];

  while (p != NULL) 
      {
	if ((p->sign == int('+')) && detect_curve(bond,n_bond,p))
	  {
	    n = p->curve.n;
	    tag = p->curve.tag;
	    c = p->curve.c;
	    int top=height;
	    int x1=0;
	    int left=width;
	    int y1=0;
	    int bottom=0;
	    int x2=0;
	    int right=0;
	    int y2=0;
	    for (int i=0; i<n; i++) 
	      {
		switch (tag[i]) 
		  {
		  case POTRACE_CORNER:
		    if (c[i][1].x<left) {left=int(c[i][1].x);y1=int(c[i][1].y);}
		    if (c[i][1].x>right) {right=int(c[i][1].x);y2=int(c[i][1].y);}
		    if (c[i][1].y<top) {top=int(c[i][1].y);x1=int(c[i][1].x);}
		    if (c[i][1].y>bottom) {bottom=int(c[i][1].y);x2=int(c[i][1].x);}
		    break;
		  case POTRACE_CURVETO:
		    if (c[i][0].x<left) {left=int(c[i][0].x);y1=int(c[i][0].y);}
		    if (c[i][0].x>right) {right=int(c[i][0].x);y2=int(c[i][0].y);}
		    if (c[i][0].y<top) {top=int(c[i][0].y);x1=int(c[i][0].x);}
		    if (c[i][0].y>bottom) {bottom=int(c[i][0].y);x2=int(c[i][0].x);}
		    if (c[i][1].x<left) {left=int(c[i][1].x);y1=int(c[i][1].y);}
		    if (c[i][1].x>right) {right=int(c[i][1].x);y2=int(c[i][1].y);}
		    if (c[i][1].y<top) {top=int(c[i][1].y);x1=int(c[i][1].x);}
		    if (c[i][1].y>bottom) {bottom=int(c[i][1].y);x2=int(c[i][1].x);}
		    break;
		  }
		if (c[i][2].x<left) {left=int(c[i][2].x);y1=int(c[i][2].y);}
		if (c[i][2].x>right) {right=int(c[i][2].x);y2=int(c[i][2].y);}
		if (c[i][2].y<top) {top=int(c[i][2].y);x1=int(c[i][2].x);}
		if (c[i][2].y>bottom) {bottom=int(c[i][2].y);x2=int(c[i][2].x);}
	      }


	    if (((bottom-top)<=max_font_height) && 
		((right-left)<=max_font_width) && (right-left>1)
		&& (right-left)<avg/2
		)
	    {
	      double aspect=1.*(bottom-top)/(right-left);
	      double fill=1.*p->area/((bottom-top)*(right-left));
	      char c=' ';
	      bool char_to_right=false;
	      bool inside_char=false;
	      for(int j=0;j<n_letters;j++)
		{
		  if (letters[j].x>right && (top+bottom)/2>letters[j].y-letters[j].r
		      && (top+bottom)/2<letters[j].y+letters[j].r
		      && right>letters[j].x-2*letters[j].r
		      && letters[j].a!='-' && letters[j].a!='+')
		    char_to_right=true;
		  if (letters[j].x-letters[j].r<=left 
		      && letters[j].x+letters[j].r>=right
		      && letters[j].y-letters[j].r<=top 
		      && letters[j].y+letters[j].r>=bottom)
		    inside_char=true;
		}
	      //cout<<left<<","<<y1<<" "<<right<<","<<y2<<" "<<top<<","<<x1<<" "<<bottom<<","<<x2<<endl;

	      if (aspect<0.7 && fill>0.9 && !char_to_right && !inside_char)  c='-';
	      else if (aspect>0.7 && aspect<1./0.7 
		       && abs(y1-y2)<3 && abs(y1+y2-bottom-top)/2<3
		       && abs(x1-x2)<3 && abs(x1+x2-right-left)/2<3
		       && !inside_char
		       //&& !char_to_right
		       )
		c='+';
	      if (c!=' ')
		{
		  letters[n_letters].a=c;
		  letters[n_letters].x=(left+right)/2;
		  letters[n_letters].y=(top+bottom)/2;
		  letters[n_letters].r=distance(left,top,right,bottom)/2;
		  letters[n_letters].free=true;
		  n_letters++;
		  if (n_letters>=MAX_ATOMS) n_letters--;
		  delete_curve_with_children(atom,bond,n_atom,n_bond,p);
		}
	    }
	  }
	p = p->next;
      }
  return(n_letters);	 
}

void  find_old_aromatic_bonds(potrace_path_t *p,bond_t *bond,int n_bond,
			      atom_t *atom,int n_atom, double avg)
{
  potrace_path_t *p1=p;
  for(int i=0;i<n_bond;i++)
    if (bond[i].exists)
      bond[i].arom=false;
  while (p != NULL) 
    {
      if ((p->sign == int('-')) && detect_curve(bond,n_bond,p))
	{
	  potrace_path_t *child=p->childlist;
	  if (child != NULL && child->sign == int('+'))
	    {
	      potrace_path_t *gchild=child->childlist;
	      if (gchild != NULL && gchild->sign == int('-'))
		{
		  for(int i=0;i<n_bond;i++)
		    if (bond[i].exists && bond[i].curve==p)
		      bond[i].arom=true;
		  delete_curve_with_children(atom,bond,n_atom,n_bond,child);
		}
	    }
	}
      p = p->next;
    }
    
  while (p1 != NULL) 
    {
      if (p1->sign == int('+') && detect_curve(bond,n_bond,p1))
	{
	  potrace_path_t *child=p1->childlist;
	  if (child != NULL && child->sign == int('-'))
	    {
	      vector<int> vert;
	      double circum=0;
	      for(int i=0;i<n_bond;i++)
		if (bond[i].exists && bond[i].curve==p1)
		  circum+=bond_length(bond,i,atom);
	      for(int i=0;i<n_atom;i++)
		if (atom[i].exists && atom[i].curve==p1)
		  vert.push_back(i);
	      if (vert.size()>4)
		{
		  double diameter=0,center_x=0,center_y=0;
		  int num=0;
		  for(unsigned int i=0;i<vert.size();i++)
		    {
		      for(unsigned int j=i+1;j<vert.size();j++)
			{
			  double dist=distance(atom[vert[i]].x,atom[vert[i]].y,
					       atom[vert[j]].x,atom[vert[j]].y);
			  if(dist>diameter)
			    diameter=dist;
			}
		      center_x+=atom[vert[i]].x;
		      center_y+=atom[vert[i]].y;
		      num++;
		    }
		  center_x/=num;
		  center_y/=num;
		  bool centered=true;
		  for(unsigned int i=0;i<vert.size();i++)
		    {
		      double dist=distance(atom[vert[i]].x,atom[vert[i]].y,
					   center_x,center_y);
		      if(fabs(dist-diameter/2)>V_DISPLACEMENT)
			centered=false;
		    }

		  if (circum<PI*diameter && diameter>avg/2 && diameter<3*avg
		      && centered)
		    {
		      delete_curve_with_children(atom,bond,n_atom,n_bond,p1);
		      for(int i=0;i<n_bond;i++)
			if (bond[i].exists)
			  {
			    double dist=distance(
					 (atom[bond[i].a].x+atom[bond[i].b].x)/2,
					 (atom[bond[i].a].y+atom[bond[i].b].y)/2,
					 center_x,center_y);
			    double ang=angle4(atom[bond[i].b].x,atom[bond[i].b].y,
					      atom[bond[i].a].x,atom[bond[i].a].y,
					      center_x,center_y,
					      atom[bond[i].a].x,atom[bond[i].a].y);
			    ang=acos(ang)* 180.0 / PI;
			    if (ang<90 && dist<(avg/3+diameter/2))
			      {
				bond[i].arom=true;
			      }
			  }
		    }
		
		}
	    }
	}
      p1 = p1->next;
    }
  
}

void flatten_bonds(bond_t *bond,int n_bond,atom_t *atom,int n_atom,double maxh)
{
  bool found=true;


  while (found)
    {
      found=false;
      for (int i=0;i<n_bond;i++)
	if (bond[i].exists && bond[i].type<3)
	  {
	    int n=0;
	    int f=i;
	    for (int j=0;j<n_bond;j++)
	      if (j!=i && bond[j].exists && bond[j].type<3 &&
		  (bond[i].a==bond[j].a || bond[i].a==bond[j].b))
		{
		  n++;
		  f=j;
		}
	    if (n==1)
	      {
		if (bond[i].a==bond[f].b)
		  {
		    double h=fabs(distance_from_bond_y(atom[bond[i].a].x,
						       atom[bond[i].a].y,
						       atom[bond[i].b].x,
						       atom[bond[i].b].y,
						       atom[bond[f].a].x,
						       atom[bond[f].a].y));
		    if (h<=maxh)
		      {
			bond[f].exists=false;
			atom[bond[f].b].exists=false;
			bond[i].a=bond[f].a;
			bond[i].type=max(bond[i].type,bond[f].type);
			if (bond[f].arom) bond[i].arom=true;
			if (bond[f].hash) bond[i].hash=true;
			if (bond[f].wedge) bond[i].wedge=true;
			found=true;
		      }
		  }
		else 
		  {
		    double h=fabs(distance_from_bond_y(atom[bond[i].a].x,
						       atom[bond[i].a].y,
						       atom[bond[i].b].x,
						       atom[bond[i].b].y,
						       atom[bond[f].b].x,
						       atom[bond[f].b].y));
		    if (h<=maxh)
		      {

			bond[f].exists=false;
			atom[bond[f].a].exists=false;
			
			if (bond[f].hash || bond[f].wedge)
			  {
			    bond[i].a=bond[i].b;
			    bond[i].b=bond[f].b;
			  }
			else bond[i].a=bond[f].b;
			bond[i].type=max(bond[i].type,bond[f].type);
			if (bond[f].arom) bond[i].arom=true;
			if (bond[f].hash) bond[i].hash=true;
			if (bond[f].wedge) bond[i].wedge=true;
			found=true;
		      }
		  }
	      }

	    n=0;
	    f=i;
	    for (int j=0;j<n_bond;j++)
	      if (j!=i && bond[j].exists && bond[j].type<3 &&
		  (bond[i].b==bond[j].a || bond[i].b==bond[j].b))
		{
		  n++;
		  f=j;
		}
	    if (n==1)
	      {
		if (bond[i].b==bond[f].b)
		  {
		    double h=fabs(distance_from_bond_y(atom[bond[i].a].x,
						       atom[bond[i].a].y,
						       atom[bond[i].b].x,
						       atom[bond[i].b].y,
						       atom[bond[f].a].x,
						       atom[bond[f].a].y));
		    if (h<=maxh)
		      {
			bond[f].exists=false;
			atom[bond[f].b].exists=false;
			if (bond[f].hash || bond[f].wedge)
			  {
			    bond[i].b=bond[i].a;
			    bond[i].a=bond[f].a;
			  }
			else bond[i].b=bond[f].a;
			bond[i].type=max(bond[i].type,bond[f].type);
			if (bond[f].arom) bond[i].arom=true;
			if (bond[f].hash) bond[i].hash=true;
			if (bond[f].wedge) bond[i].wedge=true;
			found=true;
		      }
		  }
		else 
		  {
		    double h=fabs(distance_from_bond_y(atom[bond[i].a].x,
						       atom[bond[i].a].y,
						       atom[bond[i].b].x,
						       atom[bond[i].b].y,
						       atom[bond[f].b].x,
						       atom[bond[f].b].y));
		    if (h<=maxh)
		      {
			bond[f].exists=false;
			atom[bond[f].a].exists=false;
			bond[i].b=bond[f].b;
			bond[i].type=max(bond[i].type,bond[f].type);
			if (bond[f].arom) bond[i].arom=true;
			if (bond[f].hash) bond[i].hash=true;
			if (bond[f].wedge) bond[i].wedge=true;
			found=true;
		      }
		  }
	      }

	  
	  }
    }
  
}

void clean_unrecognized_characters(bond_t *bond,int n_bond,atom_t *atom,
				   int real_font_height, int real_font_width,
				   unsigned int size)
{
  vector<int> all_bonds(n_bond,0);
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists) 
      all_bonds[i]=1;

  for (int i=0;i<n_bond;i++)
    if (all_bonds[i]==1)
      {
	list<int> bag;
	list<int> trash;
	all_bonds[i]=2;
	bag.push_back(i);
	while (!bag.empty())
	  {
	    int k=bag.front();
	    bag.pop_front();
	    all_bonds[k]=0;
	    for (int j=0;j<n_bond;j++)
	      if (j!=k && all_bonds[j]==1 && 
		  (bond[k].a==bond[j].a || bond[k].a==bond[j].b ||
		   bond[k].b==bond[j].a || bond[k].b==bond[j].b))
		{
		  all_bonds[j]=2;
		  bag.push_back(j);
		}
	  }
      double t=FLT_MAX,b=0,l=FLT_MAX,r=0;
      while (!bag.empty())
	{
	  int k=bag.front();
	  bag.pop_front();
	  trash.push_back(k);

	  if (atom[bond[k].a].x<l) l=atom[bond[k].a].x;
	  if (atom[bond[k].b].x<l) l=atom[bond[k].b].x;
	  if (atom[bond[k].a].x>r) r=atom[bond[k].a].x;
	  if (atom[bond[k].b].x>r) r=atom[bond[k].b].x;
	  if (atom[bond[k].a].y<t) t=atom[bond[k].a].y;
	  if (atom[bond[k].b].y<t) t=atom[bond[k].b].y;
	  if (atom[bond[k].a].y>b) b=atom[bond[k].a].y;
	  if (atom[bond[k].b].y>b) b=atom[bond[k].b].y;
	}
      if ((r-l)<real_font_width && (b-t)<real_font_height && trash.size()>size)
	while (!trash.empty())
	  {
	    int k=trash.front();
	    trash.pop_front();
	    bond[k].exists=false;
	  }
      }
}

job_t *JOB;

int main(int argc,char **argv)
{
    fclose(stderr);
    srand(1);
    TCLAP::CmdLine cmd("OSRA: Optical Structure Recognition, created by Igor Filippov, 2007-2008",' ',OSRA_VERSION);
    TCLAP::UnlabeledValueArg<string>  input( "in", "input file",true,"", "filename"  );
    cmd.add(input);
    TCLAP::ValueArg<double> threshold("t","threshold","Gray level threshold",
				      false,0,"0.2..0.8");
    cmd.add(threshold);
    TCLAP::ValueArg<string> output("o","output","Write out images to files",false,"","filename prefix");
    cmd.add(output);
    TCLAP::ValueArg<int> resolution_param("r","resolution","Resolution in dots per inch",false,0,"default: auto");
    cmd.add(resolution_param);
    TCLAP::SwitchArg inv("n","negate","Invert color (white on black)",false);
    cmd.add(inv);
    TCLAP::ValueArg<string> resize("s","size","Resize image on output",false,"","dimensions, 300x400");
    cmd.add(resize);
    TCLAP::SwitchArg conf("p","print","Print out confidence estimate",false);
    cmd.add(conf);
    TCLAP::SwitchArg guess("g","guess","Print out resolution guess",false);
    cmd.add(guess);
    cmd.parse( argc, argv );

    int input_resolution=resolution_param.getValue();
    string type=image_type(input.getValue());
    bool invert=inv.getValue();
    if ((type=="PDF") || (type=="PS")) input_resolution=150;
    int page=count_pages(input.getValue());

    int image_count=0;
    int ttt=1;

    for(int l=0;l<page;l++)
      {
	Image image;
	image.density("150x150");
	stringstream pname;
	pname<<input.getValue()<<"["<<l<<"]";
	image.read(pname.str());
	
	image.modifyImage();
	image.type( TrueColorType );
	if (!invert)
	  {
	    double a=0;
	    ColorRGB c;
	    for (int i=0;i<BG_PICK_POINTS;i++)
	      {
		int x=(image.columns()*rand())/RAND_MAX;
		int y=(image.rows()*rand())/RAND_MAX;
		c=image.pixelColor(x,y);
		a+=(c.red()+c.green()+c.blue())/3;
	      }
	    a/=BG_PICK_POINTS;
	    if (a<0.5) invert=true;
	  }
	for (unsigned int i=0;i<image.columns();i++)
	  for (unsigned int j=0;j<image.rows();j++)
	    {
	      ColorRGB c,b;
	      b=image.pixelColor(i,j);
	      double a=min(b.red(),min(b.green(),b.blue()));
	      if (invert)
		a=max(b.red(),max(b.green(),b.blue()));
	      c.red(a);c.green(a);c.blue(a);
	      image.pixelColor(i,j,c);
	    }
	image.contrast(2);
	image.type( GrayscaleType );
	
	int num_resolutions=NUM_RESOLUTIONS;
	if (input_resolution!=0) num_resolutions=1;
	vector<int> select_resolution(num_resolutions,input_resolution);
	vector < vector <string> > array_of_smiles(num_resolutions);
	vector<double> array_of_confidence(num_resolutions,-FLT_MAX);
	vector< vector <Image> >  array_of_images(num_resolutions);
	    
	if (input_resolution==0)
	  {
	    select_resolution[0]=72;
	    select_resolution[1]=150;
	    select_resolution[2]=300;
	  }
	int res_iter;
#pragma omp parallel for default(shared) shared(threshold,invert,output,resize,type,page,l,num_resolutions,select_resolution,array_of_smiles,array_of_confidence,array_of_images,image,image_count,conf,guess) private(res_iter,JOB)
    for (res_iter=0;res_iter<num_resolutions;res_iter++)
      {
	int n_boxes=0,total_boxes=0;
	double total_confidence=0;
	box_t boxes[NUM_BOXES];

	int resolution=select_resolution[res_iter];
	int working_resolution=resolution;

	potrace_param_t *param;
	param = potrace_param_default();
	param->alphamax=0.;
	//    param->turnpolicy=POTRACE_TURNPOLICY_MINORITY;
	param->turdsize=0;

	double THRESHOLD_BOND,THRESHOLD_CHAR;
	THRESHOLD_BOND=threshold.getValue();
	if (THRESHOLD_BOND<0.0001)
	  {
	    if (resolution>=150)
	      {
		THRESHOLD_BOND=THRESHOLD_GLOBAL;
	      }
	    else 
	      {
		THRESHOLD_BOND=0.2;
	      }
	  }
	THRESHOLD_CHAR=THRESHOLD_BOND;

	    if (resolution>300)
	      {
		int percent=(100*300)/resolution;
		stringstream scale;
		scale<<percent<<"%";
		image.scale(scale.str());
		working_resolution=300;
	      }

	    ColorGray bgColor=getBgColor(image,invert);


	    int width=image.columns();
	    int height=image.rows();
	    int max_font_height=MAX_FONT_HEIGHT*resolution/150;
	    int max_font_width=MAX_FONT_WIDTH*resolution/150;
	    int boundary=2*BOUNDARY;
	    bool thick=true;
	    if (resolution<300)	boundary=BOUNDARY;
	    if (resolution<=150) thick=false;

	    n_boxes=find_boxes(boxes,image,THRESHOLD_BOND,bgColor,width,height,
			       COARSE_GRAIN,boundary,working_resolution);
	    qsort(boxes,n_boxes,sizeof(box_t),comp_boxes);
	    for (int k=0;k<n_boxes;k++)
	      {
		int n_atom=0,n_bond=0,n_letters=0,n_label=0;
		atom_t atom[MAX_ATOMS];
		bond_t bond[MAX_ATOMS];
		letters_t letters[MAX_ATOMS];
		label_t label[MAX_ATOMS];
		potrace_bitmap_t *bm;
		potrace_path_t *p;
		potrace_state_t *st;

		Image orig_box=image;

		width=orig_box.columns();
		height=orig_box.rows();


		if (boxes[k].x1>MIN_FONT_HEIGHT 
		    || boxes[k].y1>MIN_FONT_HEIGHT 
		    || width-boxes[k].x2>MIN_FONT_HEIGHT 
		    || height-boxes[k].y2>MIN_FONT_HEIGHT)
		  orig_box.crop(Geometry(boxes[k].x2-boxes[k].x1,
					 boxes[k].y2-boxes[k].y1,
					 boxes[k].x1,boxes[k].y1));
		width=orig_box.columns();
		height=orig_box.rows();
		Image thick_box;


		if (resolution>=300)
		  {
		    double nf=noise_factor(orig_box,width,height,bgColor,
					   THRESHOLD_BOND);
		    if (nf>0.5 && nf<1.)
		      {
			thick_box=anisotropic_smoothing(orig_box,width,height);
		      }
		      else thick_box=orig_box;
		  }
		else if (resolution<300 && resolution>150)
		  {
		    int nw=width*300/resolution;
		    int nh=height*300/resolution;
		    thick_box=anisotropic_scaling(orig_box,width,height,nw,nh);
		    width=thick_box.columns();
		    height=thick_box.rows();
		    int percent=(100*300)/resolution;
		    stringstream scale;
		    scale<<percent<<"%";
		    orig_box.scale(scale.str());
		    working_resolution=300;
		  }
		else 
		    thick_box=orig_box;
		
		    
		param->turnpolicy=POTRACE_TURNPOLICY_MINORITY;
		double c_width=1.*width*72/working_resolution;
		double c_height=1.*height*72/working_resolution;
		if (c_height*c_width<SMALL_PICTURE_AREA)
		  param->turnpolicy=POTRACE_TURNPOLICY_BLACK;
	

		Image box;
		if (thick)
		  box=thin_image(thick_box,THRESHOLD_BOND,bgColor);
		else  box=thick_box;
	   
	    
		bm = bm_new(width,height);
		for(int i=0;i<width;i++)
		  for(int j=0;j<height;j++)
		    BM_PUT(bm,i,j,getPixel(box,bgColor,i,j,THRESHOLD_BOND));
	
		st = potrace_trace(param, bm);
		p = st->plist;
		n_atom=find_atoms(p,atom,bond,&n_bond);
		int real_font_width,real_font_height;
		
		n_letters=find_chars(p,orig_box,letters,atom,bond,n_atom,n_bond,
				     height,width,bgColor,THRESHOLD_CHAR,
				     max_font_width,max_font_height,
				     real_font_width,real_font_height);
	
		double avg_bond=percentile75(bond,n_bond,atom);
		if (working_resolution==300)
		  {
		    n_letters=find_fused_chars(bond,n_bond,atom,letters,n_letters,
					       real_font_height,real_font_width,
					       avg_bond/3,orig_box,bgColor,
					       THRESHOLD_CHAR);
		  }
		double max_area=avg_bond*5;
		if (thick) max_area=avg_bond;
		n_letters=find_plus_minus(p,letters,atom,bond,n_atom,n_bond,
					  height,width,real_font_height,
					  real_font_width,n_letters,avg_bond);

	
		n_atom=find_small_bonds(p,atom,bond,n_atom,&n_bond,
					max_area,avg_bond/2,5);
		find_old_aromatic_bonds(p,bond,n_bond,atom,n_atom,avg_bond);
		
		
		double thickness=skeletize(atom,bond,n_bond,box,THRESHOLD_BOND,bgColor);
		remove_disconnected_atoms(atom,bond,n_atom,n_bond);
		n_bond=fix_one_sided_bonds(bond,n_bond,atom,2);
		collapse_atoms(atom,bond,n_atom,n_bond,3);
		remove_zero_bonds(bond,n_bond,atom);
		flatten_bonds(bond,n_bond,atom,n_atom,3);
		remove_zero_bonds(bond,n_bond,atom);
		avg_bond=percentile75(bond,n_bond,atom);

		avg_bond=percentile75(bond,n_bond,atom);
		double max_dist_double_bond=0;
		n_bond=double_triple_bonds(atom,bond,n_bond,avg_bond,n_atom,
					   5,max_dist_double_bond);


		n_atom=find_dashed_bonds(p,atom,bond,n_atom,&n_bond,
					 max(MAX_DASH,int(avg_bond/3)),
					 avg_bond,orig_box,bgColor,
					 THRESHOLD_BOND,thick);



		n_letters=remove_small_bonds(bond,n_bond,atom,letters,n_letters,
					     real_font_height,MIN_FONT_HEIGHT,avg_bond);

		clean_unrecognized_characters(bond,n_bond,atom,
					      real_font_height,real_font_width,1);

		thickness=find_wedge_bonds(thick_box,atom,n_atom,bond,n_bond,bgColor,
					   THRESHOLD_BOND,max_dist_double_bond,avg_bond);

		n_label=assemble_labels(letters,n_letters,label);


		remove_disconnected_atoms(atom,bond,n_atom,n_bond);		
		collapse_atoms(atom,bond,n_atom,n_bond,thickness);
		remove_zero_bonds(bond,n_bond,atom);
		flatten_bonds(bond,n_bond,atom,n_atom,2*thickness);
		remove_zero_bonds(bond,n_bond,atom);


		avg_bond=percentile75(bond,n_bond,atom);

		collapse_double_bonds(bond,n_bond,atom,n_atom,max_dist_double_bond);
		extend_terminal_bond_to_label(atom,letters,n_letters,bond,n_bond,
					      label,n_label,avg_bond,
					      thickness,max_dist_double_bond);

		remove_disconnected_atoms(atom,bond,n_atom,n_bond);
		collapse_atoms(atom,bond,n_atom,n_bond,thickness);
		remove_zero_bonds(bond,n_bond,atom);
		flatten_bonds(bond,n_bond,atom,n_atom,thickness);
		remove_zero_bonds(bond,n_bond,atom);
		remove_disconnected_atoms(atom,bond,n_atom,n_bond);


		extend_terminal_bond_to_bonds(atom,n_atom,bond,n_bond,avg_bond,
					      2*thickness,max_dist_double_bond);


		collapse_atoms(atom,bond,n_atom,n_bond,3);
		remove_zero_bonds(bond,n_bond,atom);
		flatten_bonds(bond,n_bond,atom,n_atom,3);
		remove_zero_bonds(bond,n_bond,atom);
		clean_unrecognized_characters(bond,n_bond,atom,
					      real_font_height,real_font_width,0);
					


		assign_charge(atom,bond,n_atom,n_bond);
		find_up_down_bonds(bond,n_bond,atom,thickness);
		int real_atoms=count_atoms(atom,n_atom);
		if ((real_atoms>MIN_A_COUNT) && (real_atoms<MAX_A_COUNT))
		  {
		    int f=resolve_bridge_bonds(atom,n_atom,bond,n_bond);
		    
		    flatten_bonds(bond,n_bond,atom,n_atom,min(real_font_width,
							      real_font_height));
		    remove_zero_bonds(bond,n_bond,atom);
		   debug(thick_box,atom,n_atom,bond,n_bond,"tmp.png");
		    int rotors,rings;
		    double confidence=0;
		    string smiles=get_smiles(atom,bond,n_bond,rotors,confidence,f,rings);
		    if (f<5 && smiles!="")
		      {
			array_of_smiles[res_iter].push_back(smiles);
			total_boxes++;
			total_confidence+=confidence;
			array_of_images[res_iter].push_back(orig_box);
		      }
		  }

		potrace_state_free(st);
		free(bm);
	      }
	    if (total_boxes>0) 
	      array_of_confidence[res_iter]=total_confidence/total_boxes;
	    potrace_param_free(param); 
      }
    double max_conf=-FLT_MAX;
    int max_res=0;
    for (int i=0;i<num_resolutions;i++)
      {
	if (array_of_confidence[i]>max_conf && array_of_smiles[i].size()>0)
	  {
	    max_conf=array_of_confidence[i];
	    max_res=i;
	  }
      }

    for (unsigned int i=0;i<array_of_smiles[max_res].size();i++)
      {
	cout<<array_of_smiles[max_res][i];
	if (guess.getValue()) cout<<" "<<select_resolution[max_res];
	if (conf.getValue()) cout<<" "<<max_conf;
	cout<<endl;
	stringstream fname;
	if (output.getValue()!="") fname<<output.getValue()<<image_count<<".png";
	image_count++;
	if (fname.str()!="")
	  {
	    Image tmp=array_of_images[max_res][i];
	    if (resize.getValue()!="")
	      {
		tmp.scale(resize.getValue());
	      }
	    tmp.write(fname.str());
	  }
      }
   
   }

   

  return 0;
}
