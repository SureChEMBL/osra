extern "C" {
#include "pgm2asc.h"
}

#include <algorithm>
#include <cstdio>
#include <vector>

#include "./common.h"
#include "./rectangle.h"
#include "./ucs.h"
#include "./bitmap.h"
#include "./block.h"
#include "./character.h"

#include "osra.h"


char get_atom_label(Image image, ColorGray bg, int x1, int y1, int x2, int y2, double THRESHOLD)
{
  Control control;


  char c=0,c1=0;
  unsigned char* tmp;
  job_t job;
  double f=1.;
  JOB=&job;
  job_init(&job);
  job.cfg.cfilter="oOcCnNHFsSBuUgMeEXYZRPp23568h";

  //job.cfg.cs=160;
  //job.cfg.certainty=80;
  //job.cfg.dust_size=1;
  if ((y2-y1)>MAX_FONT_HEIGHT) f=1.*MAX_FONT_HEIGHT/(y2-y1);


  job.src.p.x=int(f*(x2-x1+1));
  job.src.p.y=int(f*(y2-y1+1));
  job.src.p.bpp=1;
  job.src.p.p = (unsigned char *)malloc(job.src.p.x*job.src.p.y);

  Block *b=new Block(0,0,job.src.p.x,job.src.p.y);

  tmp=(unsigned char *)malloc(int((x2-x1+1)*(y2-y1+1)));

  for(int i=0;i<job.src.p.x*job.src.p.y;i++) job.src.p.p[i]=255;

  for (int i=y1;i<=y2;i++)
    for (int j=x1;j<=x2;j++)
      tmp[(i-y1)*(x2-x1+1)+j-x1]=(unsigned char)(255-255*getPixel(image,bg,j,i,THRESHOLD));


  int t=1;
  int y=0;
  int x=int((x2-x1+1)/2);
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
	      c2=a.result(control);
	      //cout<<"c2="<<c2<<endl;
	      string patern=job.cfg.cfilter;
	      if (patern.find(c2,0)==string::npos) c2='_';
	      if (isalnum(c2)) c=c2;
	    }


	}
      //cout<<c<<endl<<"=========================="<<endl;
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
