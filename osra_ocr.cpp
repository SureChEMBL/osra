/*********************************************************************
  OSRA: Optical Structure Recognition

  This is a U.S. Government work (year) and is therefore not subject to copyright.  
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


extern "C" {
#include "pgm2asc.h"
}

#include <algorithm>
#include <cstdio>
#include <vector>

#include "common.h"
#include "rectangle.h"
#include "ucs.h"
#include "bitmap.h"
//#include "block.h"
#include "blob.h"
#include "character.h"

#include "osra.h"


char get_atom_label(Magick::Image image, Magick::ColorGray bg, int x1, int y1, int x2, int y2, double THRESHOLD, int dropx,int dropy)
{
  Control control;


  char c=0,c1=0;
  unsigned char* tmp;
  job_t job;
  double f=1.;
  JOB=&job;
  job_init(&job);
  job.cfg.cfilter="oOcCnNHFsSBuUgMeEXYZRPp23568";

  //job.cfg.cs=160;
  //job.cfg.certainty=80;
  //job.cfg.dust_size=1;
//  if ((y2-y1)>MAX_FONT_HEIGHT) f=1.*MAX_FONT_HEIGHT/(y2-y1);


  job.src.p.x=int(f*(x2-x1+1));
  job.src.p.y=int(f*(y2-y1+1));
  job.src.p.bpp=1;
  job.src.p.p = (unsigned char *)malloc(job.src.p.x*job.src.p.y);

  Blob *b=new Blob(0,0,job.src.p.x,job.src.p.y);

  tmp=(unsigned char *)malloc(int((x2-x1+1)*(y2-y1+1)));

  for(int i=0;i<job.src.p.x*job.src.p.y;i++) job.src.p.p[i]=255;

  for (int i=y1;i<=y2;i++)
    for (int j=x1;j<=x2;j++)
      tmp[(i-y1)*(x2-x1+1)+j-x1]=(unsigned char)(255-255*getPixel(image,bg,j,i,THRESHOLD));


  int t=1;
  int y=dropy-y1+1;
  int x=dropx-x1;
  //int y=0;
  //int x=int((x2-x1+1)/2);
  while ((t!=0) && (y<int(y2-y1+1)))
    {
      t=tmp[y*(x2-x1+1)+x];
      y++;
    }
  y--;
  if (t==0)
    {
      tmp[y*(x2-x1+1)+x]=2;
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
	  tmp[y*(x2-x1+1)+x]=1;
	  for(int i=x-1;i<x+2;i++)
	    for (int j=y-1;j<y+2;j++)
	      if ((i<(x2-x1+1)) && (j<(y2-y1+1)) && (i>=0) && (j>=0) &&
		  (tmp[j*(x2-x1+1)+i]==0))
		{
		  cx.push_back(i);
		  cy.push_back(j);
		  tmp[j*(x2-x1+1)+i]=2;
		}
	}

      for (int i=0;i<(y2-y1+1);i++)
	for(int j=0;j<(x2-x1+1);j++)
	  if (tmp[i*(x2-x1+1)+j]==1) 
	    {
	      tmp[i*(x2-x1+1)+j]=0;
	    }
	  else tmp[i*(x2-x1+1)+j]=255;

      int count=0;
      int zeros=0;
      for (int i=y1;i<=y2;i++)
	{
	  for (int j=x1;j<=x2;j++)
	    {
	      int x=int(f*(j-x1));
	      int y=int(f*(i-y1));
	      if ((x<job.src.p.x) && (y<job.src.p.y) && (job.src.p.p[y*job.src.p.x+x]==255))
		{
		  job.src.p.p[y*job.src.p.x+x]= tmp[(i-y1)*(x2-x1+1)+j-x1];
		  if (tmp[(i-y1)*(x2-x1+1)+j-x1]==0) 
		    {
		      b->set_bit(y,x,true);
		      if(x>0 && x<job.src.p.x-1 && y>0 && y<job.src.p.y-1) count++;
		    }
		  else 
		    if(x>0 && x<job.src.p.x-1 && y>0 && y<job.src.p.y-1)
		      zeros++;
		}

	    }
	}
/*      cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;           
      for (int i=0;i<job.src.p.y;i++)
	{
	  for(int j=0;j<job.src.p.x;j++)
	    cout<<job.src.p.p[i*job.src.p.x+j]/255;
	  cout<<endl;
	  }
  */    
      if (count>MIN_CHAR_POINTS && zeros>MIN_CHAR_POINTS)
	{
	  try {
	    pgm2asc(&job);
	  }
	  catch(...){}
	  char *l;
	  l=(char *)job.res.linelist.start.next->data;
	  if (l!=NULL)  c1=l[0];
	  //cout<<"c1="<<c1<<endl;
	  if (isalnum(c1)) c=c1;
	  else
	    {
	      char c2=0;
	      b->find_holes();
	      Character a(b);
	      a.recognize1(control.charset,Rectangle::Rectangle( a.left(), a.top(), a.right(), a.bottom()));
	      c2=a.result();
	      //cout<<"c2="<<c2<<endl;
	      string patern=job.cfg.cfilter;
	      if (patern.find(c2,0)==string::npos) c2='_';
	      if (isalnum(c2)) c=c2;
	    }


	}
      //cout<<c<<endl;//<<"=========================="<<endl;
    }
      job_free(&job);
    
  if (isalnum(c))
    {
      return (c);
    }
  else
    {
      return(0);
    }
}



string fix_atom_name(string s,int n)
{
  string r=s;
  //cout<<s<<"-->";
  if (s.length()==1) r=toupper(s.at(0));
  if (s=="Ci" || s=="Cf" || s=="Cll") r="Cl";
  else if (s=="H" && n>1) r="N";
  else if (s=="HN" || s=="NH" || s=="M" || s=="Hm" || s=="MN" || s=="N2"
	   || s=="NM" || s=="NH2" || s=="H2N" || s=="NHZ" || s=="HZN" || s=="NH3"
	   || s=="nu")   r="N";
  else if (s=="OH" || s=="oH" || s=="Ho" || s=="HO" || s=="ol"
	   || s=="On" || s=="on" || s=="no" || s=="nO") r="O";
  else if (s=="Meo" || s=="oMe" || s=="oMg" || s=="omg" || s=="Mgo"
	   || s=="leo" || s=="ohle" || s=="lleo" || s=="olllle" 
	   || s=="OMe" || s=="OM8" || s=="OMo" || s=="OMB" || s=="OCH3" 
	   || s=="OCHS")   r="MeO";
  else if (s=="NC" || s=="YC")  r="CN";
  else if ((s=="nBU") || (s=="neU") ||(s=="ngU")) r="nBu";
  else if ((s=="Eto") || (s=="oEt") || (s=="Elo") || (s=="oEl")
      || s=="ElO") r="EtO";
  else if ((s=="olgU") || (s=="oleU") || s=="OlBU") r="OiBu";
  else if ((s=="npr") || (s=="llpll") || (s=="lpl") || (s=="npl")
      || s=="lPl") r="iPr";
  else if ((s=="tBU") || (s=="BU") || (s=="llBU") || (s=="lBU")) r="tBu";
  else if (s=="CooH" || s=="HooC" || s=="Co2H" || s=="CO2H" || s=="HOOC" || s=="CO2n") r="COOH";
  else if (s=="AC" || s=="pC" || s=="pc") r="Ac";
  else if (s=="ACo" || s=="opC" || s=="pcO" || s=="ACO" || s=="oCO" 
	   || s=="OoC" || s=="OpC" || s=="pCO" || s=="RCO" || s=="ORC") r="AcO";
  else if (s=="Bl" || s=="el") r="Br";
  else if (s=="CH3" || s=="H3C") r="C";
  else if (s=="R" || s=="Rl" || s=="Rlo" || s=="R2" || s=="R3" || s=="Rg"
      || s=="R4" || s=="R5" || s=="R6" || s=="R7" || s=="R8" || s=="Z"  || s=="Y" 
      || s=="2" || s=="RlO") r="X";
  else if (s=="pl" || s=="nl") r="Ar";
  else if (s=="oX") r="Ox";
  else if (s=="NoZ" || s=="o2N" || s=="No2" || s=="No" || s=="O2N"
      || s=="NOZ" || s=="MO2") r="NO2";
//  else if (s=="ph" || s=="Pl" || s=="pl") r="Ph";
  else if (s=="F3C" || s=="CF" || s=="FC" || s=="Co" || s=="F8l" || s=="CFS"
      || s=="FSC") r="CF3";
  else if (s=="F3Co") r="F3CN";
  else if (s=="S3" || s=="Se" || s=="lS" || s=="8") r="S";
  else if (s=="lH") r="H";
  else if (s=="NHnC") r="NHAc";
  else if (s=="OlHP" || s=="lHPO") r="THPO";
  else if (s=="NlOHCH3") r="NOHCH3";
  else if (s=="HO3S") r="SO3H";
  else if (s=="NMe" || s=="NHMe") r="MeN";
  else if (s=="RO") r="OR";
  else if (s=="lHPO" || s=="OlHP") r="THPO";
  else if (s=="NCOlRlH3") r="N(OH)CH3";
  else if (s=="pZO" || s=="p2O" || s=="OBX" || s=="BZO" || s=="B2O" || s=="OB2") r="BzO";
  //cout<<r<<endl;
  return(r);
}

unsigned char Character::result() const throw()
{
   if( guesses() )
    {
      const unsigned char ch = UCS::map_to_byte( gv[0].code );
      if( ch ) return ch;
    }
  return '_';
}
                          
