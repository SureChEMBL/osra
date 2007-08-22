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
#include <cstdio>
#include <vector>

#include "common.h"
#include "rational.h"
#include "rectangle.h"
#include "page_image.h"


namespace {

void enlarge_2( std::vector< std::vector< unsigned char > > & data ) throw()
  {
  const int height = data.size(), width = data[0].size();
  std::vector< std::vector< unsigned char > > new_data( 2 * height );

  for( unsigned int row = 0; row < new_data.size(); ++row )
    new_data[row].resize( 2 * width, 1 );

  for( int row = 0; row < height; ++row )
    {
    const std::vector< unsigned char > & datarow = data[row];
    std::vector< unsigned char > & new_datarow0 = new_data[2*row];
    std::vector< unsigned char > & new_datarow1 = new_data[2*row+1];
    for( int col = 0; col < width; ++col )
      {
      if( datarow[col] == 0 )
        {
        const bool l = col > 0 && datarow[col-1] == 0;
        const bool t = row > 0 && data[row-1][col] == 0;
        const bool r = col < width - 1 && datarow[col+1] == 0;
        const bool b = row < height - 1 && data[row+1][col] == 0;
        const bool lt = row > 0 && col > 0 && data[row-1][col-1] == 0;
        const bool rt = row > 0 && col < width - 1 && data[row-1][col+1] == 0;
        const bool lb = row < height - 1 && col > 0 && data[row+1][col-1] == 0;
        const bool rb = row < height - 1 && col < width - 1 && data[row+1][col+1] == 0;

        if( l || t || lt || ( !rt && !lb ) ) new_datarow0[2*col] = 0;
        if( r || t || rt || ( !lt && !rb ) ) new_datarow0[2*col+1] = 0;
        if( l || b || lb || ( !lt && !rb ) ) new_datarow1[2*col] = 0;
        if( r || b || rb || ( !rt && !lb ) ) new_datarow1[2*col+1] = 0;
        }
      }
    }
  data.swap( new_data );
  }


void enlarge_3( std::vector< std::vector< unsigned char > > & data ) throw()
  {
  const int height = data.size(), width = data[0].size();
  std::vector< std::vector< unsigned char > > new_data( 3 * height );

  for( unsigned int row = 0; row < new_data.size(); ++row )
    new_data[row].resize( 3 * width, 1 );

  for( int row = 0; row < height; ++row )
    {
    const int row3 = 3 * row;
    const std::vector< unsigned char > & datarow = data[row];
    std::vector< unsigned char > & new_datarow0 = new_data[row3];
    std::vector< unsigned char > & new_datarow1 = new_data[row3+1];
    std::vector< unsigned char > & new_datarow2 = new_data[row3+2];
    for( int col = 0; col < width; ++col )
      {
      const int col3 = 3 * col;
      const bool l = col > 0 && datarow[col-1] == 0;
      const bool t = row > 0 && data[row-1][col] == 0;
      const bool r = col < width - 1 && datarow[col+1] == 0;
      const bool b = row < height - 1 && data[row+1][col] == 0;
      const bool lt = row > 0 && col > 0 && data[row-1][col-1] == 0;
      const bool rt = row > 0 && col < width - 1 && data[row-1][col+1] == 0;
      const bool lb = row < height - 1 && col > 0 && data[row+1][col-1] == 0;
      const bool rb = row < height - 1 && col < width - 1 && data[row+1][col+1] == 0;
      if( datarow[col] == 0 )
        {
        if( l || t || lt || ( !rt && !lb ) ) new_datarow0[col3] = 0;
        new_datarow0[col3+1] = 0;
        if( r || t || rt || ( !lt && !rb ) ) new_datarow0[col3+2] = 0;
        new_datarow1[col3] = new_datarow1[col3+1] = new_datarow1[col3+2] = 0;
        if( l || b || lb || ( !lt && !rb ) ) new_datarow2[col3] = 0;
        new_datarow2[col3+1] = 0;
        if( r || b || rb || ( !rt && !lb ) ) new_datarow2[col3+2] = 0;
        }
      else
        {
        if( l && t && lt && ( !rt || !lb ) ) new_datarow0[col3] = 0;
        if( r && t && rt && ( !lt || !rb ) ) new_datarow0[col3+2] = 0;
        if( l && b && lb && ( !lt || !rb ) ) new_datarow2[col3] = 0;
        if( r && b && rb && ( !rt || !lb ) ) new_datarow2[col3+2] = 0;
        }
      }
    }
  data.swap( new_data );
  }


void enlarge_n( std::vector< std::vector< unsigned char > > & data, const int n ) throw()
  {
  if( n < 2 ) return;
  const int height = data.size(), width = data[0].size();
  std::vector< std::vector< unsigned char > > new_data;
  new_data.reserve( n * height );

  for( int row = 0; row < height; ++row )
    {
    const std::vector< unsigned char > & datarow = data[row];
    new_data.push_back( std::vector< unsigned char >() );
    for( int col = 0; col < width; ++col )
      {
      const unsigned char d = datarow[col];
      for( int i = 0; i < n; ++i ) new_data.back().push_back( d );
      }
    for( int i = 1; i < n; ++i ) new_data.push_back( new_data.back() );
    }
  data.swap( new_data );
  }


void mirror_left_right( std::vector< std::vector< unsigned char > > & data ) throw()
  {
  const int height = data.size();
  for( int row = 0; row < height; ++row )
    std::reverse( data[row].begin(), data[row].end() );
  }


void mirror_top_bottom( std::vector< std::vector< unsigned char > > & data ) throw()
  {
  for( int u = 0, d = data.size() - 1; u < d; ++u, --d )
    data[u].swap( data[d] );
  }


void mirror_diagonal( std::vector< std::vector< unsigned char > > & data,
                      Rectangle & re ) throw()
  {
  const int size = std::max( re.height(), re.width() );

  if( re.height() < size )
    {
    data.resize( size );
    for( int row = re.height(); row < size; ++row )
      data[row].resize( size );
    }
  else if( re.width() < size )
    for( int row = 0; row < re.height(); ++row )
      data[row].resize( size );

  for( int row = 0; row < size; ++row )
    {
    std::vector< unsigned char > & datarow = data[row];
    for( int col = 0; col < row; ++col )
      {
      unsigned char tmp = datarow[col];
      datarow[col] = data[col][row]; data[col][row] = tmp;
      }
    }

  const int h = re.height(), w = re.width();
  re.height( w ); re.width( h );
  if( re.height() < size ) data.resize( re.height() );
  else if( re.width() < size )
    for( int row = 0; row < re.height(); ++row )
      data[row].resize( re.width() );
  }

} // end namespace


bool Page_image::scale( int n ) throw( Page_image::Error )
  {
  if( n <= -2 )
    { Page_image reduced( *this, -n ); *this = reduced; return true; }
  if( n >= 2 )
    {
    if( (long long)width() * height() * n > (long long)INT_MAX )
      throw Error( "scale factor too big. `int' will overflow." );
    if( _maxval == 1 )
      {
      if( n && ( n % 2 ) == 0 ) { enlarge_2( data ); n /= 2; }
      else if( n && ( n % 3 ) == 0 ) { enlarge_3( data ); n /= 3; }
      }
    if( n > 1 ) enlarge_n( data, n );
    Rectangle::height( data.size() );
    Rectangle::width( data[0].size() );
    return true;
    }
  return false;
  }


void Page_image::transform( const Transformation & t ) throw()
  {
  switch( t.type() )
    {
    case Transformation::none:
      break;
    case Transformation::rotate90:
      mirror_diagonal( data, *this ); mirror_top_bottom( data ); break;
    case Transformation::rotate180:
      mirror_left_right( data ); mirror_top_bottom( data ); break;
    case Transformation::rotate270:
      mirror_diagonal( data, *this ); mirror_left_right( data ); break;
    case Transformation::mirror_lr:
      mirror_left_right( data ); break;
    case Transformation::mirror_tb:
      mirror_top_bottom( data ); break;
    case Transformation::mirror_d1:
      mirror_diagonal( data, *this ); break;
    case Transformation::mirror_d2:
      mirror_diagonal( data, *this );
      mirror_left_right( data ); mirror_top_bottom( data ); break;
    }
  }
