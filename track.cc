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
#include <vector>

#include "common.h"
#include "rectangle.h"
#include "track.h"


namespace {

void error( const char * s ) throw() __attribute__ ((noreturn));
void error( const char * s ) throw()
  { Ocrad::internal_error( s ); }


int good_reference( const Rectangle & r1, const Rectangle & r2, int & val,
                    int mean_height, int mean_width ) throw()
  {
  if( 4 * r1.height() >= 3 * mean_height &&
      4 * r2.height() >= 3 * mean_height &&
      ( r1.width() >= mean_width || r2.width() >= mean_width ) && val > 0 )
    {
    if( 4 * r1.height() <= 5 * mean_height &&
        4 * r2.height() <= 5 * mean_height )
      {
      if( 9 * r1.height() <= 10 * mean_height &&
          9 * r2.height() <= 10 * mean_height &&
          10 * std::abs( r1.bottom() - r2.bottom() ) <= mean_height )
        { val = 0; return ( r1.height() <= r2.height() ) ? 0 : 1; }
      if( val > 1 && 10 * std::abs( r1.vcenter() - r2.vcenter() ) <= mean_height )
        { val = 1; return ( r1.bottom() <= r2.bottom() ) ? 0 : 1; }
      }
    if( val > 2 && 10 * std::abs( r1.vcenter() - r2.vcenter() ) <= mean_height )
      { val = 2; return ( r1.bottom() <= r2.bottom() ) ? 0 : 1; }
    }
  return -1;
  }


int set_l( const std::vector< Rectangle > & rectangle_vector,
           int mean_height, int mean_width ) throw()
  {
  const int rectangles = rectangle_vector.size();
  const int imax = rectangles / 4;
  int ibest = -1, val = 3;
  for( int i1 = 0; i1 < imax && val > 0; ++i1 )
    for( int i2 = i1 + 1; i2 <= imax && i2 <= i1 + 2; ++i2 )
      {
      int i = good_reference( rectangle_vector[i1], rectangle_vector[i2],
                              val, mean_height, mean_width );
      if( i >= 0 ) { ibest = (i == 0) ? i1 : i2; if( val == 0 ) break; }
      }
  return ibest;
  }


int set_r( const std::vector< Rectangle > & rectangle_vector,
           int mean_height, int mean_width ) throw()
  {
  const int rectangles = rectangle_vector.size();
  const int imin = rectangles - 1 - ( rectangles / 4 );
  int ibest = -1, val = 3;
  for( int i1 = rectangles - 1; i1 > imin && val > 0; --i1 )
    for( int i2 = i1 - 1; i2 >= imin && i2 >= i1 - 2; --i2 )
      {
      int i = good_reference( rectangle_vector[i1], rectangle_vector[i2],
                              val, mean_height, mean_width );
      if( i >= 0 ) { ibest = (i == 0) ? i1 : i2; if( val == 0 ) break; }
      }
  return ibest;
  }

} // end namespace


Vrhomboid::Vrhomboid( int l, int lc, int r, int rc, int h ) throw()
  {
  if( r < l || h <= 0 )
    {
    std::fprintf( stderr, "l = %d, lc = %d, r = %d, rc = %d, h = %d\n",
                  l, lc, r, rc, h );
    error( "bad parameter building a Vrhomboid" );
    }
  _left = l; _lvcenter = lc; _right = r; _rvcenter = rc; _height = h;
  }


void Vrhomboid::left( int l ) throw()
  {
  if( l > _right ) error( "left, bad parameter resizing a Vrhomboid" );
  _left = l;
  }


void Vrhomboid::right( int r ) throw()
  {
  if( r < _left ) error( "right, bad parameter resizing a Vrhomboid" );
  _right = r;
  }


void Vrhomboid::height( int h ) throw()
  {
  if( h <= 0 ) error( "height, bad parameter resizing a Vrhomboid" );
  _height = h;
  }


void Vrhomboid::extend_left( int l ) throw()
  {
  if( l > _right )
    error( "extend_left, bad parameter resizing a Vrhomboid" );
  _lvcenter = vcenter( l ); _left = l;
  }


void Vrhomboid::extend_right( int r ) throw()
  {
  if( r < _left )
    error( "extend_right, bad parameter resizing a Vrhomboid" );
  _rvcenter = vcenter( r ); _right = r;
  }


int Vrhomboid::vcenter( int col ) const throw()
  {
  int dx = _right - _left, dy = _rvcenter - _lvcenter;
  int vc = _lvcenter;
  if( dx && dy ) vc += ( dy * ( col - _left ) ) / dx;
  return vc;
  }


bool Vrhomboid::includes( const Rectangle & r ) const throw()
  {
  if( r.left() < _left || r.right() > _right ) return false;
  int tl = top( r.left() ), bl = tl +_height - 1;
  int tr = top( r.right() ), br = tr +_height - 1;
  int t = std::max( tl, tr ), b = std::min( bl, br );
  return ( t <= r.top() && b >= r.bottom() );
  }


bool Vrhomboid::includes( int row, int col ) const throw()
  {
  if( col < _left || col > _right ) return false;
  int t = top( col ), b = t +_height - 1;
  return ( t <= row && b >= row );
  }


// rectangle_vector must be ordered by increasing hcenter().
//
void Track::set_track( const std::vector< Rectangle > & rectangle_vector ) throw()
  {
  if( data.size() ) data.clear();
  if( !rectangle_vector.size() ) return;

  const int rectangles = rectangle_vector.size();
  int mean_vcenter = 0, mean_height = 0, mean_width = 0;

  for( int i = 0; i < rectangles; ++i )
    {
    mean_vcenter += rectangle_vector[i].vcenter();
    mean_height += rectangle_vector[i].height();
    mean_width += rectangle_vector[i].width();
    }
  if( rectangles )
    { mean_vcenter /= rectangles; mean_height /= rectangles; mean_width /= rectangles; }

  // short line
  if( rectangles < 8 )
    {
    const Rectangle & r1 = rectangle_vector[0];
    const Rectangle & r2 = rectangle_vector[rectangles-1];
    data.push_back( Vrhomboid( r1.hcenter(), mean_vcenter,
                               r2.hcenter(), mean_vcenter, mean_height ) );
    return;
    }

  // look for references
  int l = set_l( rectangle_vector, mean_height, mean_width );
  int r = set_r( rectangle_vector, mean_height, mean_width );

  int lcol, lvc, rcol, rvc;
  if( l >= 0 )
    {
    lcol = rectangle_vector[l].hcenter();
    lvc = Vrhomboid::vcenter( rectangle_vector[l].bottom(), mean_height );
    }
  else { lcol = rectangle_vector.front().hcenter(); lvc = mean_vcenter; }
  if( r >= 0 )
    {
    rcol = rectangle_vector[r].hcenter();
    rvc = Vrhomboid::vcenter( rectangle_vector[r].bottom(), mean_height );
    }
  else { rcol = rectangle_vector.back().hcenter(); rvc = mean_vcenter; }

  data.push_back( Vrhomboid( lcol, lvc, rcol, rvc, mean_height ) );
  if( l != 0 ) data.front().extend_left( rectangle_vector[0].hcenter() );
  if( r != rectangles - 1 )
    data.back().extend_right( rectangle_vector[rectangles - 1].hcenter() );
  }


int Track::vcenter( int col ) const throw()
  {
  for( unsigned int i = 0; i < data.size(); ++i )
    {
    const Vrhomboid & vr = data[i];
    if( col <= vr.right() || i >= data.size() - 1 ) return vr.vcenter( col );
    }
  return 0;
  }


int Track::bottom( int col ) const throw()
  {
  for( unsigned int i = 0; i < data.size(); ++i )
    {
    const Vrhomboid & vr = data[i];
    if( col <= vr.right() || i >= data.size() - 1 ) return vr.bottom( col );
    }
  return 0;
  }


int Track::top( int col ) const throw()
  {
  for( unsigned int i = 0; i < data.size(); ++i )
    {
    const Vrhomboid & vr = data[i];
    if( col <= vr.right() || i >= data.size() - 1 ) return vr.top( col );
    }
  return 0;
  }
