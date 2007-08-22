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
#include <cstdio>
#include <vector>

#include "common.h"
#include "rational.h"
#include "rectangle.h"
#include "track.h"
#include "page_image.h"


void Page_image::find_columns( const Rectangle & rin, bool recursive ) throw()
  {
  if( !this->includes( rin ) ) return;
  const int colmin = ( width() > 100 ) ? width() / 10 : 10;
  const int gapmin = ( width() > 300 ) ? width() / 100 : 3;
  const int min_width = ( 2 * colmin ) + gapmin;
  if( rin.width() < min_width ) { rv.push_back( rin ); return; }
  const unsigned int rv_size_orig = rv.size();
  const int ldiff = rin.left() - left();
  const int tdiff = rin.top() - top();
  std::vector< int > h_outline( rin.width(), 0 );

  for( int row = 0; row < rin.height(); ++row )
    {
    const std::vector< unsigned char > & datarow = data[row+tdiff];
    for( int col = 0; col < rin.width(); ++col )
      if( datarow[col+ldiff] == 0 ) ++h_outline[col];
    }
  int total_dots = 0;
  for( int col = 0; col < rin.width(); ++col ) total_dots += h_outline[col];
  if( 10 * total_dots > 8 * rin.size() ) return;	// eliminates images

  const int threshold_col = std::max( total_dots / ( rin.width() * 4 ), 10 );
  const int threshold_gap = std::max( total_dots / ( rin.width() * 20 ), 1 );
  int ileft = 0, iright = rin.width() - 1;

  while( ileft < iright && h_outline[ileft] <= 0 ) ++ileft; // cut left border
  if( rin.left() == left() && 10 * ileft < iright )	// cut left border noise
    {
    int l = ileft + gapmin;
    while( l < iright && h_outline[l] <= 0 ) ++l;
    if( l > ileft + ( 2 * gapmin ) ) ileft = l;
    }
  ileft = std::max( 0, ileft - gapmin );

  while( iright > ileft && h_outline[iright] <= 0 ) --iright; // cut right border
  iright = std::min( iright + gapmin, rin.width() - 1 );

  while( iright - ileft >= 2 * colmin )
    {
    int l, r;
    for( l = r = ileft; r < iright; ++r )		// look for text
      if( h_outline[r] < threshold_col )
        { if( r - l >= colmin ) break; else l = r; }
    if( r - l < colmin ) break;

    for( l = r; r < iright; ++r )			// look for gap
      if( h_outline[r] > threshold_gap )
        { if( r - l >= gapmin ) break; else l = r; }
    if( r - l < gapmin ) break;

    if( r < iright )			// cut by a minimum near the center
      {
      int mid = ( r + l ) / 2, half = ( r - l ) / 2;
      r = mid;
      for( int i = 1; i <= half && h_outline[r] > 0; ++i )
        {
        if( h_outline[mid+i] < h_outline[r] ) r = mid + i;
        if( h_outline[mid-i] < h_outline[r] ) r = mid - i;
        }
      }
    Rectangle re( rin.left() + ileft, rin.top(), rin.left() + r, rin.bottom() );
    if( recursive ) find_rows( re, re.width() >= min_width );
    else rv.push_back( re );
    ileft = r;
    }
  if( iright - ileft > gapmin )
    {
    Rectangle re( rin.left() + ileft, rin.top(), rin.left() + iright, rin.bottom() );
    if( recursive && ( rv.size() > rv_size_orig || rv.size() == 0 ) )
      find_rows( re, re.width() >= min_width );
    else rv.push_back( re );
    }
  }


void Page_image::find_rows( const Rectangle & rin, bool recursive ) throw()
  {
  if( !this->includes( rin ) ) return;
  const int rowmin = ( height() > 100 ) ? height() / 10 : 10;
  const int gapmin = ( height() > 300 ) ? height() / 100 : 3;
  const int min_height = ( 2 * rowmin ) + gapmin;
  if( rin.height() < min_height ) { rv.push_back( rin ); return; }
  const unsigned int rv_size_orig = rv.size();
  const int ldiff = rin.left() - left();
  const int tdiff = rin.top() - top();
  std::vector< int > v_outline( rin.height(), 0 );

  int total_dots = 0;
  for( int row = 0; row < rin.height(); ++row )
    {
    const std::vector< unsigned char > & datarow = data[row+tdiff];
    for( int col = 0; col < rin.width(); ++col )
      if( datarow[col+ldiff] == 0 ) ++v_outline[row];
    total_dots += v_outline[row];
    }
  if( 10 * total_dots > 8 * rin.size() ) return;	// eliminates images

  const int threshold_gap = ( total_dots / ( rin.height() * 20 ) ) + 1;
  int itop = 0, ibottom = rin.height() - 1;

  while( itop < ibottom && v_outline[itop] <= 0 ) ++itop;	// cut top border
  itop = std::max( 0, itop - gapmin );
  while( ibottom > itop && v_outline[ibottom] <= 0 ) --ibottom; // cut bottom border
  ibottom = std::min( ibottom + gapmin, rin.height() - 1 );

  while( ibottom - itop >= min_height )
    {
    int t, b;					// top and bottom of gap
    for( t = b = itop + gapmin; t < ibottom - gapmin; ++t )
      if( v_outline[t] < threshold_gap )
        {
        for( b = t + 1; b < ibottom && v_outline[b] < threshold_gap; ++b );
        if( b - t >= gapmin ) break; else t = b;
        }
    if( b - t < gapmin ) break;

    if( b < ibottom )			// cut by a minimum near the center
      {
      int mid = ( b + t ) / 2, half = ( b - t ) / 2;
      b = mid;
      for( int i = 1; i <= half && v_outline[b] > 0; ++i )
        {
        if( v_outline[mid+i] < v_outline[b] ) b = mid + i;
        if( v_outline[mid-i] < v_outline[b] ) b = mid - i;
        }
      }
    Rectangle re( rin.left(), rin.top() + itop, rin.right(), rin.top() + b );
    if( recursive ) find_columns( re, re.height() >= min_height && ibottom - b > gapmin );
    else rv.push_back( re );
    itop = b;
    }
  if( ibottom - itop > gapmin )
    {
    Rectangle re( rin.left(), rin.top() + itop, rin.right(), rin.top() + ibottom );
    if( recursive && rv.size() > rv_size_orig )
      find_columns( re, re.height() >= min_height );
    else rv.push_back( re );
    }
  }


// Creates a reduced Page_image
//
Page_image::Page_image( const Page_image & source, const int scale ) throw()
  : Rectangle( source ), _maxval( source._maxval ), _threshold( source._threshold )
  {
  if( scale < 2 || scale > source.width() || scale > source.height() )
    Ocrad::internal_error( "bad parameter building a reduced Page_image" );

  const int scale2 = scale * scale;
  Rectangle::height( source.height() / scale );
  Rectangle::width( source.width() / scale );

  data.resize( height() );
  for( int row = 0; row < height(); ++row )
    {
    const int srow = ( row * scale ) + scale;
    data[row].reserve( width() );
    std::vector< unsigned char > & datarow = data[row];
    for( int col = 0; col < width(); ++col )
      {
      const int scol = ( col * scale ) + scale;
      int sum = 0;
      for( int i = srow - scale; i < srow; ++i )
        {
        const std::vector< unsigned char > & sdatarow = source.data[i];
        for( int j = scol - scale; j < scol; ++j )
          sum += sdatarow[j];
        }
      datarow.push_back( sum / scale2 );
      }
    }
  }


// Creates a reduced, b/w Page_image
//
Page_image::Page_image( const Page_image & source, const int scale, const Rational & th ) throw()
  : Rectangle( source ), _maxval( 1 ), _threshold( 0 )
  {
  if( scale < 2 || scale > source.width() || scale > source.height() )
    Ocrad::internal_error( "bad parameter building a reduced Page_image" );

  const int ith = ( th >= 0 && th <= 1 ) ?
                  ( th * ( scale * scale ) ).trunc() :
                  ( ( scale * scale ) - 1 ) / 2;
  Rectangle::height( source.height() / scale );
  Rectangle::width( source.width() / scale );

  data.resize( height() );
  for( int row = 0; row < height(); ++row )
    {
    const int srow = ( row * scale ) + scale;
    data[row].reserve( width() );
    std::vector< unsigned char > & datarow = data[row];
    for( int col = 0; col < width(); ++col )
      {
      const int scol = ( col * scale ) + scale;
      int counter = 0;
      for( int i = srow - scale; i < srow; ++i )
        {
        const std::vector< unsigned char > & sdatarow = source.data[i];
        for( int j = scol - scale; j < scol; ++j )
          if( sdatarow[j] <= source._threshold && ++counter > ith ) goto L1;
        }
      L1: datarow.push_back( ( counter > ith ) ? 0 : 1 );
      }
    }
  }


int Page_image::analyse_layout( const int layout_level ) throw()
  {
  rv.clear();

  if( layout_level >= 1 && layout_level <= 2 &&
      left() == 0 && top() == 0 && width() > 200 && height() > 200 )
    {
    Page_image reduced( *this, 10, Rational( 9, 100 ) );
    reduced.find_columns( reduced, layout_level >= 2 );
    rv = reduced.rv;
    for( unsigned int i = 0; i < rv.size(); ++i )
      rv[i].enlarge( 10 );
    }
  if( rv.size() == 0 ) rv.push_back( *this );
  tv.clear(); tv.insert( tv.end(), rv.size(), _threshold );
  return rv.size();
  }


void Page_image::draw_rectangle( const Rectangle & re ) throw()
  {
  int l = std::max( left(), re.left() );
  int t = std::max( top(), re.top() );
  int r = std::min( right(), re.right() );
  int b = std::min( bottom(), re.bottom() );
  if( l == re.left() )
    for( int row = t; row <= b; ++row ) set_bit( row, l, true );
  if( t == re.top() )
    for( int col = l; col <= r; ++col ) set_bit( t, col, true );
  if( r == re.right() )
    for( int row = t; row <= b; ++row ) set_bit( row, r, true );
  if( b == re.bottom() )
    for( int col = l; col <= r; ++col ) set_bit( b, col, true );
  }


void Page_image::draw_track( const Track & tr ) throw()
  {
  int l = std::max( left(), tr.left() );
  int r = std::min( right(), tr.right() );
  if( l == tr.left() )
    for( int row = tr.top( l ); row <= tr.bottom( l ); ++row )
      if( row >= top() && row <= bottom() ) set_bit( row, l, true );
  if( r == tr.right() )
    for( int row = tr.top( r ); row <= tr.bottom( r ); ++row )
      if( row >= top() && row <= bottom() ) set_bit( row, r, true );
  for( int col = l; col <= r; ++col )
    {
    int row = tr.top( col );
    if( row >= top() && row <= bottom() ) set_bit( row, col, true );
    row = tr.bottom( col );
    if( row >= top() && row <= bottom() ) set_bit( row, col, true );
    }
  }


void Page_image::histogramize( const bool vertical ) throw()
  {
  if( vertical )
    {
    for( int col = left(); col <= right(); ++col )
      for( int row = top(), y = bottom(); row <= bottom(); ++row )
        if( get_bit( row, col ) )
          { set_bit( row, col, false ); set_bit( y--, col, true ); }
    }
  else
    {
    for( int row = top(); row <= bottom(); ++row )
      for( int col = left(), x = left(); col <= right(); ++col )
        if( get_bit( row, col ) )
          { set_bit( row, col, false ); set_bit( row, x++, true ); }
    }
  }
