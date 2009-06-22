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

bool detect_bracket(int x, int y,unsigned char *pic)
{
  Control control;
  char c1=0;
  job_t job;
  JOB=&job;
  job_init(&job);
  job.cfg.cfilter="([{";

  //job.cfg.cs=160;
  //job.cfg.certainty=80;
  //job.cfg.dust_size=1;

  bool res=false;

  job.src.p.x=x;
  job.src.p.y=y;
  job.src.p.bpp=1;
  job.src.p.p = pic;

  Blob *b=new Blob(0,0,job.src.p.x,job.src.p.y);

  int count=0;
  int zeros=0;
  for (int i=0;i<=y;i++)
    for (int j=0;j<=x;j++)
      {
	if (pic[i*x+j]==0) 
	  {
	    b->set_bit(y,x,true);
	    count++;
	  }
	else 
	  zeros++;
      }

  /* for (int i=0;i<job.src.p.y;i++)
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
	  if (c1=='(' || c1=='[' || c1=='{') res=true;
	  else
	    {
	      char c2=0;
	      b->find_holes();
	      Character a(b);
	      a.recognize1(control.charset,Rectangle::Rectangle( a.left(), a.top(), a.right(), a.bottom()));
	      c2=a.result();
	      if (c2=='(' || c2=='[' || c2=='{') res=true;
	    }
	}

      job_free(&job);
      
      return(res);
}


string fix_atom_name(string s,int n,map<string,string> fix, bool debug)
{
  string r=s;
  if (s.length()==1) r=toupper(s.at(0));
  if (s=="H" && n>1) r="N";

  map<string,string>::iterator it=fix.find(s);
  string mapped=" ";
  if (it!=fix.end())   
    {
      r=it->second;
      mapped=r;
    }


  if (debug && s!=" " && s!="") 
    cout<<s<<" --> "<<mapped<<" --> "<<r<<endl;

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
                          
