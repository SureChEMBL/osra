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
#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include "common.h"
#include "rectangle.h"


namespace {

void error( const char * s ) throw() __attribute__ ((noreturn));
void error( const char * s ) throw()
  { Ocrad::internal_error( s ); }


int hypoti( const int c1, const int c2 )
  {
  long long temp = c1; temp *= temp;
  long long target = c2; target *= target; target += temp;
  int lower = std::max( std::abs(c1), std::abs(c2) );
  int upper = std::abs(c1) + std::abs(c2);
  while( upper - lower > 1 )
    {
    int m = ( lower + upper ) / 2;
    temp = m; temp *= temp;
    if( temp < target ) lower = m; else upper = m;
    }
  temp = lower; temp *= temp; target *= 2; target -= temp;
  temp = upper; temp *= temp;
  if( target < temp ) return lower;
  else return upper;
  }

} // end namespace


Rectangle::Rectangle( const int l, const int t, const int r, const int b ) throw()
  {
  if( r < l || b < t )
    {
    std::fprintf( stderr, "l = %d, t = %d, r = %d, b = %d\n", l, t, r, b );
    error( "bad parameter building a Rectangle" );
    }
  _left = l; _top = t; _right = r; _bottom = b;
  }


void Rectangle::left( const int l ) throw()
  {
  if( l > _right ) error( "left, bad parameter resizing a Rectangle" );
  _left = l;
  }


void Rectangle::top( const int t ) throw()
  {
  if( t > _bottom ) error( "top, bad parameter resizing a Rectangle" );
  _top = t;
  }


void Rectangle::right( const int r ) throw()
  {
  if( r < _left ) error( "right, bad parameter resizing a Rectangle" );
  _right = r;
  }


void Rectangle::bottom( const int b ) throw()
  {
  if( b < _top ) error( "bottom, bad parameter resizing a Rectangle" );
  _bottom = b;
  }


void Rectangle::height( const int h ) throw()
  {
  if( h <= 0 ) error( "height, bad parameter resizing a Rectangle" );
  _bottom = _top + h - 1;
  }


void Rectangle::width( const int w ) throw()
  {
  if( w <= 0 ) error( "width, bad parameter resizing a Rectangle" );
  _right = _left + w - 1;
  }


void Rectangle::add_point( const int row, const int col ) throw()
  {
  if( row > _bottom ) _bottom = row; else if( row < _top ) _top = row;
  if( col > _right ) _right = col;   else if( col < _left ) _left = col;
  }


void Rectangle::add_rectangle( const Rectangle & re ) throw()
  {
  if( re._left < _left ) _left = re._left;
  if( re._top < _top ) _top = re._top;
  if( re._right > _right ) _right = re._right;
  if( re._bottom > _bottom ) _bottom = re._bottom;
  }


void Rectangle::enlarge( const int scale ) throw()
  {
  if( scale > 0 )
    { _left *= scale; _top *= scale; _right *= scale; _bottom *= scale; }
  }


void Rectangle::move( const int row, const int col ) throw()
  {
  int d = row - _top; if( d ) { _top += d; _bottom += d; }
  d = col - _left; if( d ) { _left += d; _right += d; }
  }


bool Rectangle::includes( const Rectangle & re ) const throw()
  {
  return ( _left  <= re._left  && _top    <= re._top &&
           _right >= re._right && _bottom >= re._bottom );
  }


bool Rectangle::includes( const int row, const int col ) const throw()
  {
  return ( _left <= col && _right >= col && _top <= row && _bottom >= row );
  }


bool Rectangle::strictly_includes( const Rectangle & re ) const throw()
  {
  return ( _left  < re._left  && _top    < re._top &&
           _right > re._right && _bottom > re._bottom );
  }


bool Rectangle::strictly_includes( const int row, const int col ) const throw()
  {
  return ( _left < col && _right > col && _top < row && _bottom > row );
  }


bool Rectangle::includes_hcenter( const Rectangle & re ) const throw()
  { return ( _left <= re.hcenter() && _right >= re.hcenter() ); }


bool Rectangle::includes_vcenter( const Rectangle & re ) const throw()
  { return ( _top <= re.vcenter() && _bottom >= re.vcenter() ); }


bool Rectangle::h_includes( const Rectangle & re ) const throw()
  { return ( _left <= re._left && _right >= re._right ); }


bool Rectangle::v_includes( const Rectangle & re ) const throw()
  { return ( _top <= re._top && _bottom >= re._bottom ); }


bool Rectangle::h_includes( const int col ) const throw()
  { return ( _left <= col && _right >= col ); }


bool Rectangle::v_includes( const int row ) const throw()
  { return ( _top <= row && _bottom >= row ); }


bool Rectangle::h_overlaps( const Rectangle & re ) const throw()
  { return ( _left <= re._right && _right >= re._left ); }


bool Rectangle::v_overlaps( const Rectangle & re ) const throw()
  { return ( _top <= re._bottom && _bottom >= re._top ); }


bool Rectangle::is_hcentred_in( const Rectangle & re ) const throw()
  {
  if( this->h_includes( re.hcenter() ) ) return true;
  int w = std::min( re.height(), re.width() ) / 2;
  if( width() < w )
    {
    int d = ( w + 1 ) / 2;
    if( hcenter() - d <= re.hcenter() && hcenter() + d >= re.hcenter() )
      return true;
    }
  return false;
  }


bool Rectangle::is_vcentred_in( const Rectangle & re ) const throw()
  {
  if( this->v_includes( re.vcenter() ) ) return true;
  int h = std::min( re.height(), re.width() ) / 2;
  if( height() < h )
    {
    int d = ( h + 1 ) / 2;
    if( vcenter() - d <= re.vcenter() && vcenter() + d >= re.vcenter() )
      return true;
    }
  return false;
  }


bool Rectangle::h_precedes( const Rectangle & re ) const throw()
  { return ( hcenter() < re.hcenter() ); }


bool Rectangle::v_precedes( const Rectangle & re ) const throw()
  {
  if( _bottom < re.vcenter() || vcenter() < re._top ) return true;
  if( this->includes_vcenter( re ) && re.includes_vcenter( *this ) )
    return this->h_precedes( re );
  return false;
  }


int Rectangle::distance( const Rectangle & re ) const throw()
  { return hypoti( h_distance( re ), v_distance( re ) ); }


int Rectangle::distance( const int row, const int col ) const throw()
  { return hypoti( h_distance( col ), v_distance( row ) ); }


int Rectangle::h_distance( const Rectangle & re ) const throw()
  {
  if( re._right <= _left ) return _left - re._right;
  if( re._left >= _right ) return re._left - _right;
  return 0;
  }

int Rectangle::h_distance( const int col ) const throw()
  {
  if( col <= _left ) return _left - col;
  if( col >= _right ) return col - _right;
  return 0;
  }

int Rectangle::v_distance( const Rectangle & re ) const throw()
  {
  if( re._bottom <= _top ) return _top - re._bottom;
  if( re._top >= _bottom ) return re._top - _bottom;
  return 0;
  }

int Rectangle::v_distance( const int row ) const throw()
  {
  if( row <= _top ) return _top - row;
  if( row >= _bottom ) return row - _bottom;
  return 0;
  }
