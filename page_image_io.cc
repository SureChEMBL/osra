/*  GNU Ocrad - Optical Character Recognition program
    Copyright (C) 2003, 2004, 2005, 2006 Antonio Diaz Diaz.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
*/

#include <algorithm>
#include <cctype>
#include <climits>
#include <cstdio>
#include <string>
#include <vector>

#include "common.h"
#include "rational.h"
#include "rectangle.h"
#include "page_image.h"


namespace {

// binarization by Otsu's method based on maximization of inter-class variance
//
int otsu_th( const std::vector< std::vector< unsigned char > > & data,
             const Rectangle & re, const int maxval ) throw()
  {
  if( maxval == 1 ) return 0;

  std::vector< int > hist( maxval + 1, 0 );	// histogram of image data
  for( int row = re.top(); row <= re.bottom(); ++row )
    for( int col = re.left(); col <= re.right(); ++col )
      ++hist[data[row][col]];

  std::vector< int > chist;		// cumulative histogram
  chist.reserve( maxval + 1 );
  chist.push_back( hist[0] );
  std::vector< long long > cmom;	// cumulative moment
  cmom.reserve( maxval + 1 );
  cmom.push_back( 0 );			// 0 times hist[0] equals zero
  for( int i = 1; i <= maxval; ++i )
    {
    chist.push_back( chist[i-1] + hist[i] );
    cmom.push_back( cmom[i-1] + ( i * hist[i] ) );
    }

  const double cmom_max = cmom[maxval];
  double bvar_max = 0;
  int threshold = 0;			// threshold for binarization
  for( int i = 0; i < maxval; ++i )
    if( chist[i] > 0 && chist[i] < re.size() )
      {
      double bvar = (double)cmom[i] / chist[i];
      bvar -= ( cmom_max - cmom[i] ) / ( re.size() - chist[i] );
      bvar *= bvar; bvar *= chist[i]; bvar *= ( re.size() - chist[i] );
      if( bvar > bvar_max ) { bvar_max = bvar; threshold = i; }
      }

  if( Ocrad::verbose )
    std::fprintf( stderr, "maxval = %d, automatic threshold = %d (%s)\n",
                  maxval, threshold,
                  Rational( threshold, maxval ).to_decimal( 1, -3 ).c_str() );
  return threshold;
  }


unsigned char pnm_getrawbyte( FILE * f ) throw( Page_image::Error )
  {
  int ch = std::fgetc( f );

  if( ch == EOF )
    throw Page_image::Error( "end-of-file reading pnm file." );

  return static_cast< unsigned char > (ch);
  }


unsigned char pnm_getc( FILE * f ) throw( Page_image::Error )
  {
  unsigned char ch;
  bool comment = false;

  do {
    ch = pnm_getrawbyte( f );
    if( ch == '#' ) comment = true;
    else if( ch == '\n' ) comment = false;
    }
  while( comment );
  return ch;
  }


int pnm_getint( FILE * f ) throw( Page_image::Error )
  {
  unsigned char ch;
  int i = 0;

  do ch = pnm_getc( f ); while( std::isspace( ch ) );
  if( !std::isdigit( ch ) )
    throw Page_image::Error( "junk in pnm file where an integer should be." );
  do {
    if( ( INT_MAX - (ch - '0') ) / 10 < i )
      throw Page_image::Error( "number too big in pnm file." );
    i = (i * 10) + (ch - '0');
    ch = pnm_getc( f );
    }
  while( std::isdigit( ch ) );
  return i;
  }


unsigned char pbm_getbit( FILE * f ) throw( Page_image::Error )
  {
  unsigned char ch;

  do ch = pnm_getc( f ); while( std::isspace( ch ) );

  if( ch == '0' ) return 0;
  if( ch == '1' ) return 1;
  throw Page_image::Error( "junk in pbm file where bits should be." );
  }

} // end namespace


void Page_image::read_p1( FILE * f, const bool invert ) throw( Page_image::Error )
  {
  _maxval = 1; _threshold = 0;
  const int rows = height(), cols = width();
  if( !invert )
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; ++col )
        data[row].push_back( 1 - pbm_getbit( f ) );
  else
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; ++col )
        data[row].push_back( pbm_getbit( f ) );
  }


void Page_image::read_p4( FILE * f, const bool invert ) throw( Page_image::Error )
  {
  _maxval = 1; _threshold = 0;
  const int rows = height(), cols = width();
  if( !invert )
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; )
        {
        unsigned char byte = pnm_getrawbyte( f );
        for( unsigned char mask = 0x80; mask > 0 && col < cols; mask >>= 1, ++col )
          data[row].push_back( ( byte & mask ) ? 0 : 1 );
        }
  else
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; )
        {
        unsigned char byte = pnm_getrawbyte( f );
        for( unsigned char mask = 0x80; mask > 0 && col < cols; mask >>= 1, ++col )
          data[row].push_back( ( byte & mask ) ? 1 : 0 );
        }
  }


void Page_image::read_p2( FILE * f, const Rational & th, const bool invert ) throw( Page_image::Error )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Page_image::Error( "zero maxval in pgm file." );
  _maxval = std::min( maxval, 255 );
  const int rows = height(), cols = width();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      int val = pnm_getint( f );
      if( val > maxval ) throw Page_image::Error( "value > maxval in pgm file." );
      if( invert ) val = maxval - val;
      if( maxval > 255 ) { val *= 255; val /= maxval; }
      data[row].push_back( val );
      }

  if( th >= 0 && th <= 1 ) _threshold = ( th * _maxval ).trunc();
  else _threshold = otsu_th( data, *this, _maxval );
  }


void Page_image::read_p5( FILE * f, const Rational & th, const bool invert ) throw( Page_image::Error )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Page_image::Error( "zero maxval in pgm file." );
  if( maxval > 255 ) throw Page_image::Error( "maxval > 255 in pgm \"P5\" file." );
  _maxval = maxval;
  const int rows = height(), cols = width();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      unsigned char val = pnm_getrawbyte( f );
      if( val > _maxval ) throw Page_image::Error( "value > maxval in pgm file." );
      if( invert ) val = _maxval - val;
      data[row].push_back( val );
      }

  if( th >= 0 && th <= 1 ) _threshold = ( th * _maxval ).trunc();
  else _threshold = otsu_th( data, *this, _maxval );
  }


void Page_image::read_p3( FILE * f, const Rational & th, const bool invert ) throw( Page_image::Error )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Page_image::Error( "zero maxval in ppm file." );
  _maxval = std::min( maxval, 255 );
  const int rows = height(), cols = width();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      const int r = pnm_getint( f );			// Red value
      const int g = pnm_getint( f );			// Green value
      const int b = pnm_getint( f );			// Blue value
      if( r > maxval || g > maxval || b > maxval )
        throw Page_image::Error( "value > maxval in ppm file." );
      int val;
      if( !invert ) val = std::min( r, std::min( g, b ) );
      else val = maxval - std::max( r, std::max( g, b ) );
      if( maxval > 255 ) { val *= 255; val /= maxval; }
      data[row].push_back( val );
      }

  if( th >= 0 && th <= 1 ) _threshold = ( th * _maxval ).trunc();
  else _threshold = otsu_th( data, *this, _maxval );
  }


void Page_image::read_p6( FILE * f, const Rational & th, const bool invert ) throw( Page_image::Error )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Page_image::Error( "zero maxval in ppm file." );
  if( maxval > 255 ) throw Page_image::Error( "maxval > 255 in ppm \"P6\" file." );
  _maxval = maxval;
  const int rows = height(), cols = width();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      const unsigned char r = pnm_getrawbyte( f );	// Red value
      const unsigned char g = pnm_getrawbyte( f );	// Green value
      const unsigned char b = pnm_getrawbyte( f );	// Blue value
      if( r > _maxval || g > _maxval || b > _maxval )
        throw Page_image::Error( "value > maxval in ppm file." );
      unsigned char val;
      if( !invert ) val = std::min( r, std::min( g, b ) );
      else val = _maxval - std::max( r, std::max( g, b ) );
      data[row].push_back( val );
      }

  if( th >= 0 && th <= 1 ) _threshold = ( th * _maxval ).trunc();
  else _threshold = otsu_th( data, *this, _maxval );
  }


// Creates a Page_image from a pbm, pgm or ppm file
// "P1" (pbm), "P4" (pbm RAWBITS), "P2" (pgm), "P5" (pgm RAWBITS),
// "P3" (ppm), "P6" (ppm RAWBITS) file formats are recognized.
//
Page_image::Page_image( FILE * f, const Rational & th, const bool invert ) throw( Page_image::Error )
  : Rectangle( 0, 0, 0, 0 )
  {
  unsigned char filetype = 0;

  if( pnm_getrawbyte( f ) == 'P' )
    {
    unsigned char ch = pnm_getrawbyte( f );
    if( ch >= '1' && ch <= '6' ) filetype = ch;
    }
  if( filetype == 0 )
    throw Error( "bad magic number - not a pbm, pgm or ppm file." );

  {
  int tmp = pnm_getint( f );
  if( tmp == 0 ) throw Error( "zero width in pnm file." );
  Rectangle::width( tmp );
  tmp = pnm_getint( f );
  if( tmp == 0 ) throw Error( "zero height in pnm file." );
  Rectangle::height( tmp );
  if( width() < 16 || height() < 16 )
    throw Error( "image too small. Minimum size is 16x16." );
  if( (long long)width() * height() > (long long)INT_MAX )
    throw Error( "image too big. `int' will overflow." );
  }

  data.resize( height() );
  for( unsigned int row = 0; row < data.size(); ++row )
    data[row].reserve( width() );

  if( Ocrad::verbose )
    {
    std::fprintf( stderr, "file type is P%c\n", filetype );
    std::fprintf( stderr, "file size is %dw x %dh\n", width(), height() );
    }

  switch( filetype )
    {
    case '1': read_p1( f, invert ); break;
    case '4': read_p4( f, invert ); break;
    case '2': read_p2( f, th, invert ); break;
    case '5': read_p5( f, th, invert ); break;
    case '3': read_p3( f, th, invert ); break;
    case '6': read_p6( f, th, invert ); break;
    }
  if( Ocrad::verbose && th >= 0 && th <= 1 )
    std::fprintf( stderr, "maxval = %d, manual threshold = %d (%s)\n",
                  _maxval, _threshold, th.to_decimal( 1, -3 ).c_str() );
  }


void Page_image::adapt_thresholds() throw()
  {
  if( _maxval > 1 && ( rv.size() > 1 || ( rv.size() == 1 && rv[0] != *this ) ) )
    {
    for( unsigned int i = 0; i < rv.size(); ++i )
      tv[i] = otsu_th( data, rv[i], _maxval );
    }
  }


bool Page_image::crop( const Rational ltrb[4], const Rational & th ) throw()
  {
  if( ltrb[0] < 0 || ltrb[1] < 0 || ltrb[2] < 0 || ltrb[3] < 0 ) return false;

  Rectangle re = *this;
  int a;

  if( ltrb[0] <= 1 ) a = left() + ( ltrb[0] * ( width() - 1 ) ).trunc();
  else a = ltrb[0].round();
  if( a > re.left() ) re.left( a );

  if( ltrb[1] <= 1 ) a = top() + ( ltrb[1] * ( height() - 1 ) ).trunc();
  else a = ltrb[1].round();
  if( a > re.top() ) re.top( a );

  if( ltrb[2] <= 1 ) a = left() + ( ltrb[2] * ( width() - 1 ) ).trunc();
  else a = ltrb[2].round();
  if( a < re.right() ) { if( a < re.left() ) return false; re.right( a ); }

  if( ltrb[3] <= 1 ) a = top() + ( ltrb[3] * ( height() - 1 ) ).trunc();
  else a = ltrb[3].round();
  if( a < re.bottom() ) { if( a < re.top() ) return false; re.bottom( a ); }

  if( re.width() < 16 || re.height() < 16 ) return false;

  // cropping is performed here
  if( re.bottom() < bottom() ) data.resize( re.bottom() - top() + 1 );
  if( re.right() < right() )
    {
    const int w = re.right() - left() + 1;
    for( int row = data.size() - 1; row >= 0 ; --row ) data[row].resize( w );
    }
  if( re.top() > top() )
    data.erase( data.begin(), data.begin() + re.top() - top() );
  if( re.left() > left() )
    {
    const int d = re.left() - left();
    for( int row = data.size() - 1; row >= 0 ; --row )
      data[row].erase( data[row].begin(), data[row].begin() + d );
    }
  Rectangle::left( 0 );
  Rectangle::top( 0 );
  Rectangle::right( data[0].size() - 1 );
  Rectangle::bottom( data.size() - 1 );

  if( Ocrad::verbose )
    std::fprintf( stderr, "file cropped to %dw x %dh\n", width(), height() );
  if( _maxval > 1 && ( th < 0 || th > 1 ) )
    _threshold = otsu_th( data, *this, _maxval );
  rv.clear(); tv.clear();
  return true;
  }


bool Page_image::save( FILE * f, const char filetype, const int i ) const throw()
  {
  if( filetype < '1' || filetype > '6' || i >= (int)rv.size() ) return false;
  const Rectangle & re = ( i < 0 ) ? *this : rv[i];
  if( !this->includes( re ) ) return false;
  std::fprintf( f, "P%c\n%d %d\n", filetype, re.width(), re.height() );

  if( filetype == '1' )					// pbm
    for( int row = re.top(); row <= re.bottom(); ++row )
      {
      for( int col = re.left(); col <= re.right(); ++col )
        std::putc( get_bit( row, col ) ? '1' : '0', f );
      std::putc( '\n', f );
      }
  else if( filetype == '4' )				// pbm RAWBITS
    for( int row = re.top(); row <= re.bottom(); ++row )
      {
      unsigned char byte = 0, mask = 0x80;
      for( int col = re.left(); col <= re.right(); ++col )
        {
        if( get_bit( row, col ) ) byte |= mask;
        mask >>= 1;
        if( mask == 0 ) { std::putc( byte, f ); byte = 0; mask = 0x80; }
        }
      if( mask != 0x80 ) std::putc( byte, f ); // incomplete byte at end of row
      }
  else if( filetype == '2' )				// pgm
    for( int row = re.top(); row <= re.bottom(); ++row )
      {
      for( int col = re.left(); col < re.right(); ++col )
        std::fprintf( f, "%d ", data[row][col] );
      std::fprintf( f, "%d\n", data[row][re.right()] );
      }
  else if( filetype == '5' )				// pgm RAWBITS
    for( int row = re.top(); row <= re.bottom(); ++row )
      for( int col = re.left(); col <= re.right(); ++col )
        std::fprintf( f, "%c", data[row][col] );
  else if( filetype == '3' )				// ppm
    for( int row = re.top(); row <= re.bottom(); ++row )
      {
      for( int col = re.left(); col < re.right(); ++col )
        {
        const unsigned char d = data[row][col];
        std::fprintf( f, "%d %d %d ", d, d, d );
        }
      const unsigned char d = data[row][re.right()];
      std::fprintf( f, "%d %d %d\n", d, d, d );
      }
  else if( filetype == '6' )				// ppm RAWBITS
    for( int row = re.top(); row <= re.bottom(); ++row )
      for( int col = re.left(); col <= re.right(); ++col )
        {
        const unsigned char d = data[row][col];
        std::fprintf( f, "%c %c %c ", d, d, d );
        }
  return true;
  }
