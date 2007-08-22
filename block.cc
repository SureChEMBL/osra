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
#include "rectangle.h"
#include "bitmap.h"
#include "block.h"


namespace {

void delete_hole( std::vector< Bitmap * > & holep_vector,
                  std::vector< Bitmap * > & v1,
                  std::vector< Bitmap * > & v2, Bitmap * p, int i ) throw()
  {
  std::replace( v1.begin() + i, v1.end(), p, (Bitmap *) 0 );
  std::replace( v2.begin(), v2.begin() + i, p, (Bitmap *) 0 );

  i = holep_vector.size();
  while( --i >= 0 && holep_vector[i] != p );
  if( i < 0 ) Ocrad::internal_error( "delete_hole, lost hole." );
  holep_vector.erase( holep_vector.begin() + i );
  delete p;
  }


inline void join_holes( std::vector< Bitmap * > & holep_vector,
                        std::vector< Bitmap * > & v1,
                        std::vector< Bitmap * > & v2,
                        Bitmap * p1, Bitmap * p2, int i ) throw()
  {
  if( p1->top() > p2->top() )
    {
    Bitmap *temp = p1; p1 = p2; p2 = temp;
    std::replace( v2.begin(), v2.begin() + i + 1, p2, p1 );
    }
  else std::replace( v1.begin() + i, v1.end(), p2, p1 );

  i = holep_vector.size();
  while( --i >= 0 && holep_vector[i] != p2 );
  if( i < 0 ) Ocrad::internal_error( "join_holes, lost hole" );
  holep_vector.erase( holep_vector.begin() + i );

  p1->add_bitmap( *p2 );
  delete p2;
  }


void delete_outer_holes( const Rectangle & re, std::vector< Bitmap * > & bpv ) throw()
  {
  for( int i = bpv.size() - 1; i >= 0; --i )
    {
    Bitmap & h = *bpv[i];
    if( !re.strictly_includes( h ) )
      { delete &h; bpv.erase( bpv.begin() + i ); }
    }
  }

} // end namespace


Block::Block( const Block & b ) throw()
  : Bitmap( b ), bpv( b.bpv )
  {
  for( unsigned int i = 0; i < bpv.size(); ++i )
    bpv[i] = new Bitmap( *b.bpv[i] );
  }


Block & Block::operator=( const Block & b ) throw()
  {
  if( this != &b )
    {
    Bitmap::operator=( b );
    for( unsigned int i = 0; i < bpv.size(); ++i ) delete bpv[i];
    bpv = b.bpv;
    for( unsigned int i = 0; i < bpv.size(); ++i )
      bpv[i] = new Bitmap( *b.bpv[i] );
    }
  return *this;
  }


Block::~Block() throw()
  {
  for( unsigned int i = 0; i < bpv.size(); ++i ) delete bpv[i];
  }


void Block::left( const int l ) throw()
  {
  const int d = l - left();
  if( d ) { Bitmap::left( l ); if( d > 0 ) delete_outer_holes( *this, bpv ); }
  }


void Block::top( const int t ) throw()
  {
  const int d = t - top();
  if( d ) { Bitmap::top( t ); if( d > 0 ) delete_outer_holes( *this, bpv ); }
  }


void Block::right( const int r ) throw()
  {
  const int d = r - right();
  if( d ) { Bitmap::right( r ); if( d < 0 ) delete_outer_holes( *this, bpv ); }
  }


void Block::bottom( const int b ) throw()
  {
  const int d = b - bottom();
  if( d ) { Bitmap::bottom( b ); if( d < 0 ) delete_outer_holes( *this, bpv ); }
  }


const Bitmap & Block::hole( const int i ) const throw()
  {
  if( i < 0 || i >= holes() )
    Ocrad::internal_error( "hole, index out of bounds" );
  return *bpv[i];
  }


int Block::id( const int row, const int col ) const throw()
  {
  if( this->includes( row, col ) )
    {
    if( get_bit( row, col ) ) return 1;
    for( int i = 0; i < holes(); ++i )
      if( bpv[i]->includes( row, col ) && bpv[i]->get_bit( row, col ) )
        return -( i + 1 );
    }
  return 0;
  }


void Block::print( FILE * outfile ) const throw()
  {
  for( int row = top(); row <= bottom(); ++row )
    {
    for( int col = left(); col <= right(); ++col )
      {
      if( get_bit( row, col ) ) std::fprintf( outfile, " O" );
      else std::fprintf( outfile, " ·" );
      }
    std::fputs( "\n", outfile );
    }
  std::fputs( "\n", outfile );
  }


void Block::fill_hole( const int i ) throw()
  {
  if( i < 0 || i >= holes() )
    Ocrad::internal_error( "fill_hole, index out of bounds" );
  add_bitmap( *bpv[i] );
  delete bpv[i];
  bpv.erase( bpv.begin() + i );
  }


void Block::find_holes() throw()
  {
  for( unsigned int i = 0; i < bpv.size(); ++i ) delete bpv[i];
  bpv.clear();
  if( height() < 3 || width() < 3 ) return;

  std::vector< Bitmap * > old_data( width(), (Bitmap *) 0 );
  std::vector< Bitmap * > new_data( width(), (Bitmap *) 0 );

  for( int row = top(); row <= bottom(); ++row )
    {
    old_data.swap( new_data );
    new_data[0] = get_bit( row, left() ) ? this : 0;
    for( int col = left() + 1; col < right(); ++col )
      {
      const int dcol = col - left();
      if( get_bit( row, col ) ) new_data[dcol] = this;	// black point
      else						// white point
        {
        Bitmap *p;
        Bitmap *lp = new_data[dcol-1];
        Bitmap *tp = old_data[dcol];
        if( lp == 0 || tp == 0 )
          {
          p = 0;
          if( lp && lp != this )
            delete_hole( bpv, old_data, new_data, lp, dcol );
          else if( tp && tp != this )
            delete_hole( bpv, old_data, new_data, tp, dcol );
          }
        else if( lp != this ) { p = lp; p->add_point( row, col ); }
        else if( tp != this ) { p = tp; p->add_point( row, col ); }
        else
          {
          p = new Bitmap( col, row, col, row );
          p->set_bit( row, col, true );
          bpv.push_back( p );
          }
        new_data[dcol] = p;
        if( p && lp != tp && lp != this && tp != this )
          join_holes( bpv, old_data, new_data, lp, tp, dcol );
        }
      }
    if( !get_bit( row, right() ) )
      {
      Bitmap *lp = new_data[width()-2];
      if( lp && lp != this )
        delete_hole( bpv, old_data, new_data, lp, width() - 1 );
      }
    }

  for( int i = bpv.size() - 1; i >= 0; --i )	// FIXME noise holes removal
    {
    Bitmap & h = *bpv[i];
    if( this->strictly_includes( h ) &&
        ( h.height() > 4 || h.width() > 4 ||
          ( ( h.height() > 2 || h.width() > 2 ) && h.area() > 3 ) ) )
      continue;
    delete &h; bpv.erase( bpv.begin() + i );
    }
/*
  for( int i = bpv.size() - 1; i >= 0; --i )
    {
    Bitmap & h = *bpv[i];
    if( !this->strictly_includes( h ) )
      { delete &h; bpv.erase( bpv.begin() + i ); }
    if( 20 * h.height() < height() && 16 * h.width() < width() ) fill_hole( i );
//    else if( h.height() < 2 && h.width() < 2 && h.area() < 2 )
//      { delete &h; bpv.erase( bpv.begin() + i ); }
    }
*/
  }
