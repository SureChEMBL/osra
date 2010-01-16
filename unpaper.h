#include <Magick++.h>

void fromImageToStruct(Magick::Image source,struct IMAGE* image, int* type)
{
 *type = PGM;
 image->bitdepth = 8;
 image->color = FALSE;
 image->width=source.columns();
 image->height=source.rows();
 int maxColorIndex=255;
 int bytesPerLine = image->width;
 int inputSize = bytesPerLine * image->height;
 image->buffer = (unsigned char*)malloc(inputSize);
 Magick::ColorGray c;
 for (int j=0;j<image->height;j++)
   for (int i=0;i<image->width;i++)
     {
       c=source.pixelColor(i,j);
       image->buffer[j*bytesPerLine+i]=(unsigned char)(255*c.shade());
     }
 image->bufferGrayscale = image->buffer;
 image->bufferLightness = image->buffer;
 image->bufferDarknessInverse = image->buffer;
}

void fromStructToImage(Magick::Image &target,struct IMAGE* image)
{
  int bytesPerLine = image->width;
  Magick::ColorGray c;
  target.modifyImage();
  target.erase();
  for (int j=0;j<image->height;j++)
    for (int i=0;i<image->width;i++)
      {
	//	printf("%d ",image->buffer[j*bytesPerLine+i]);
	c.shade((1.*image->buffer[j*bytesPerLine+i])/255);
	target.pixelColor(i,j,c);
      }
  //  target.write("debug.png");
}
