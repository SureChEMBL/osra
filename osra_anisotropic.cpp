#define cimg_use_magick
#define cimg_plugin "greycstoration.h"
#include "CImg.h"
using namespace cimg_library;
using namespace Magick;

Image anisotropic_smoothing(Image image,int width,int height)
{

  CImg<unsigned char> source(width,height,1,1,0);
  unsigned char color[1]={0};
  unsigned char cc;
  ColorGray c;
  Image res( Geometry(width,height), "white" );
  res.type( GrayscaleType );

  for(int i=0;i<width;i++)
    for(int j=0;j<height;j++)
      {
	c=image.pixelColor(i,j);
	color[0]=(unsigned char)(255*c.shade());
	source.draw_point(i,j,color);
      }
  CImg<unsigned char> dest(source);  
  //  const float gfact = (sizeof(T)==2)?1.0f/256:1.0f;
  const float gfact=1.;
  const float amplitude=20.; // 40 20!
  const float sharpness=0.3; // 0.2! 0.3
  const float anisotropy=1.;
  const float alpha=.6; //0.6! 0.8
  const float sigma=2.; //1.1 2.!
  const float dl=0.8;
  const float da=30.;
  const float gauss_prec=2.;
  const unsigned int interp=0;
  const bool fast_approx=true;
  const unsigned int tile=512; // 512 0
  const unsigned int btile=4;
  const unsigned int threads=2; // 2 1

  dest.greycstoration_run(amplitude,sharpness,anisotropy,alpha,sigma,gfact,dl,da,gauss_prec,interp,fast_approx,tile,btile,threads);
  do {
    cimg::wait(200);
  } while (dest.greycstoration_is_running());

  for(int i=0;i<width;i++)
    for(int j=0;j<height;j++)
      {
	cc=dest(i,j);
	c.shade(1.*cc/255);
	res.pixelColor(i,j,c);
      }
  return(res);     
}

Image anisotropic_scaling(Image image,int width,int height, int nw, int nh)
{

  CImg<unsigned char> source(width,height,1,1,0);
  unsigned char color[1]={0};
  unsigned char cc;
  ColorGray c;
  Image res( Geometry(nw,nh), "white" );
  res.type( GrayscaleType );

  for(int i=0;i<width;i++)
    for(int j=0;j<height;j++)
      {
	c=image.pixelColor(i,j);
	color[0]=(unsigned char)(255*c.shade());
	source.draw_point(i,j,color);
      }

  //  const float gfact = (sizeof(T)==2)?1.0f/256:1.0f;
  const float gfact=1.;
  const float amplitude=20.; // 40 20!
  const float sharpness=0.2; // 0.2! 0.3
  const float anisotropy=1.;
  const float alpha=.6; //0.6! 0.8
  const float sigma=2.; //1.1 2.!
  const float dl=0.8;
  const float da=30.;
  const float gauss_prec=2.;
  const unsigned int interp=0;
  const bool fast_approx=true;
  const unsigned int tile=512; // 512 0
  const unsigned int btile=4;
  const unsigned int threads=2; // 2 1

  const unsigned int init=5;
  CImg<unsigned char> mask;

  mask.assign(source.dimx(),source.dimy(),1,1,255);
  mask = !mask.resize(nw,nh,1,1,4);
  source.resize(nw,nh,1,-100,init);
  CImg<unsigned char> dest(source);  

  dest.greycstoration_run(mask,amplitude,sharpness,anisotropy,alpha,sigma,gfact,dl,da,gauss_prec,interp,fast_approx,tile,btile,threads);
  do {
    cimg::wait(200);
  } while (dest.greycstoration_is_running());

  for(int i=0;i<nw;i++)
    for(int j=0;j<nh;j++)
      {
	cc=dest(i,j);
	c.shade(1.*cc/255);
	res.pixelColor(i,j,c);
      }
  return(res);     
}
