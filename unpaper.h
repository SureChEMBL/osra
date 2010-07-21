/******************************************************************************
 OSRA: Optical Structure Recognition

 This is a U.S. Government work (2007-2010) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2010 (igorf@helix.nih.gov)

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

#include <Magick++.h>

void fromImageToStruct(Magick::Image source, struct IMAGE* image, int* type) {
	*type = PGM;

	image->bitdepth = 8;
	image->color = FALSE;
	image->width = source.columns();
	image->height = source.rows();

	int maxColorIndex = 255;
	int bytesPerLine = image->width;
	int inputSize = bytesPerLine * image->height;

	image->buffer = (unsigned char*) malloc(inputSize);
	Magick::ColorGray c;

	for (int j = 0; j < image->height; j++)
		for (int i = 0; i < image->width; i++) {
			c = source.pixelColor(i, j);
			image->buffer[j * bytesPerLine + i] = (unsigned char) (255 * c.shade());
		}

	image->bufferGrayscale = image->buffer;
	image->bufferLightness = image->buffer;
	image->bufferDarknessInverse = image->buffer;
}

void fromStructToImage(Magick::Image &target, struct IMAGE* image) {
	int bytesPerLine = image->width;
	Magick::ColorGray c;

	target.modifyImage();
	target.erase();

	for (int j = 0; j < image->height; j++)
		for (int i = 0; i < image->width; i++) {
			//printf("%d ", image->buffer[j * bytesPerLine + i]);
			c.shade((1. * image->buffer[j * bytesPerLine + i]) / 255);
			target.pixelColor(i, j, c);
		}

	//target.write("debug.png");
}
