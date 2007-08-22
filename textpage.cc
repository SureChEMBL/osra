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

#include <cstdio>
#include <string>
#include <vector>

#include "common.h"
#include "rectangle.h"
#include "track.h"
#include "bitmap.h"
#include "block.h"
#include "character.h"
#include "page_image.h"
#include "textline.h"
#include "textblock.h"
#include "textpage.h"


namespace {

struct Zone
  {
  Rectangle rectangle;
  std::vector< std::vector< Block * > > blockp_matrix;
  Zone( const Rectangle & r ) : rectangle( r ) {}
  };


int blocks_in_zone( const std::vector< Zone > & zone_vector, int i ) throw()
  {
  const std::vector< std::vector< Block * > > & blockp_matrix = zone_vector[i].blockp_matrix;
  int sum = 0;
  for( unsigned int cut = 0; cut < blockp_matrix.size(); ++cut )
    sum += blockp_matrix[cut].size();
  return sum;
  }


int blocks_in_page( const std::vector< Zone > & zone_vector ) throw()
  {
  int sum = 0;
  for( unsigned int i = 0; i < zone_vector.size(); ++i )
    sum += blocks_in_zone( zone_vector, i );
  return sum;
  }


void bprint( const std::vector< Zone > & zone_vector, FILE * outfile ) throw()
  {
//  std::fprintf( outfile, "page size %dw x %dh\n", width(), height() );
  std::fprintf( outfile, "total zones in page %d\n", (int)zone_vector.size() );
  std::fprintf( outfile, "total blocks in page %d\n\n", blocks_in_page( zone_vector ) );
  for( unsigned int zindex = 0; zindex < zone_vector.size(); ++zindex )
    {
    const Rectangle & r = zone_vector[zindex].rectangle;
    const std::vector< std::vector< Block * > > & blockp_matrix = zone_vector[zindex].blockp_matrix;

    std::fprintf( outfile, "zone %d of %d\n", zindex + 1, (int)zone_vector.size() );
    std::fprintf( outfile, "zone size %dw x %dh\n", r.width(), r.height() );
    std::fprintf( outfile, "total cuts in zone %d\n", (int)blockp_matrix.size() );
    std::fprintf( outfile, "total blocks in zone %d\n\n", blocks_in_zone( zone_vector, zindex ) );

    for( unsigned int cut = 0; cut < blockp_matrix.size(); ++cut )
      {
      std::fprintf( outfile, "cut %d blocks %d\n", cut + 1, (int)blockp_matrix[cut].size() );
      for( unsigned int i = 0; i < blockp_matrix[cut].size(); ++i )
      blockp_matrix[cut][i]->print( outfile );
      }
    }
  }


inline void join_blocks( std::vector< Block * > & blockp_vector,
                         std::vector< Block * > & v1,
                         std::vector< Block * > & v2,
                         Block * p1, Block * p2, int i ) throw()
  {
  if( p1->top() > p2->top() )
    {
    Block * temp = p1; p1 = p2; p2 = temp;
    std::replace( v2.begin(), v2.begin() + i + 1, p2, p1 );
    }
  else std::replace( v1.begin() + i, v1.end(), p2, p1 );

  i = blockp_vector.size();
  while( --i >= 0 && blockp_vector[i] != p2 );
  if( i < 0 ) Ocrad::internal_error( "join_blocks, lost block" );
  blockp_vector.erase( blockp_vector.begin() + i );

  p1->add_bitmap( *p2 );
  delete p2;
  }


void ignore_abnormal_blocks( std::vector< Block * > & blockp_vector ) throw()
  {
  for( int i = blockp_vector.size() - 1; i >= 0; --i )
    {
    Block & b = *blockp_vector[i];
    if( b.height() > 35 * b.width() || b.width() > 25 * b.height() )
      { delete blockp_vector[i];
      blockp_vector.erase( blockp_vector.begin() + i ); }
    }
  }


void ignore_small_blocks( std::vector< Block * > & blockp_vector ) throw()
  {
  int to = 0, blocks = blockp_vector.size();
  for( int from = 0; from < blocks; ++from )
    {
    Block * p = blockp_vector[from];
    if( p->height() > 4 || p->width() > 4 ||
        ( ( p->height() > 2 || p->width() > 2 ) && p->area() > 5 ) )
      { blockp_vector[from] = blockp_vector[to]; blockp_vector[to] = p; ++to; }
    }
  if( to < blocks )
    {
    for( int i = to; i < blocks; ++i ) delete blockp_vector[i];
    blockp_vector.erase( blockp_vector.begin() + to, blockp_vector.end() );
    }
  }


void ignore_wide_blocks( Zone & zone ) throw()
  {
  const Rectangle & r = zone.rectangle;
  std::vector< std::vector< Block * > > & blockp_matrix = zone.blockp_matrix;

  for( unsigned int cut = 0; cut < blockp_matrix.size(); ++cut )
    {
    bool frame_found = false;
    for( unsigned int i = 0; i < blockp_matrix[cut].size(); )
      {
      std::vector< Block * > & blockp_vector = blockp_matrix[cut];
      const Block & b = *blockp_vector[i];
      if( 2 * b.width() < r.width() ) { ++i; continue; }
      blockp_vector.erase( blockp_vector.begin() + i );
      if( 4 * b.area() > b.size() )		// image, not frame
        {
        if( 5 * b.width() > 4 * r.width() && 5 * b.height() > 4 * r.height() )
          {
          for( unsigned int j = 0; j < blockp_vector.size(); ++j )
            delete blockp_vector[j];
          blockp_vector.clear(); break;
          }
        for( int j = blockp_vector.size() - 1; j >= 0; --j )
          {
          const Block & b2 = *blockp_vector[j];
          if( b.includes( b2 ) )
            { delete &b2; blockp_vector.erase( blockp_vector.begin() + j ); }
          else if( b2.top() < b.top() ) break;
          }
        }
      else frame_found = true;
      }
    if( frame_found )	// Make cuts from blocks inside deleted frame(s)
      {
      int bottom = 0;
      for( unsigned int i = 0; i < blockp_matrix[cut].size(); ++i )
        {
        std::vector< Block * > & blockp_vector = blockp_matrix[cut];
        const Block & b = *blockp_vector[i];
        if( b.bottom() > bottom )
          {
          int old_bottom = bottom; bottom = b.bottom();
          if( b.top() > old_bottom && i > 0 )
            {
            std::vector< Block * > new_blockp_vector( blockp_vector.begin() + i, blockp_vector.end() );
            blockp_vector.erase( blockp_vector.begin() + i, blockp_vector.end() );
            blockp_matrix.insert( blockp_matrix.begin() + cut + 1, new_blockp_vector );
            ++cut; i = 0;
            }
          }
        }
      }
    }
  }


void remove_top_bottom_noise( std::vector< Block * > & blockp_vector ) throw()
  {
  int blocks = blockp_vector.size();
  for( int i = 0; i < blocks; ++i )
    {
    Block & b = *blockp_vector[i];
    if( b.height() < 11 ) continue;

    int c = 0;
    for( int col = b.left(); col <= b.right(); ++col )
      if( b.get_bit( b.top(), col ) && ++c > 1 ) break;
    if( c <= 1 ) b.top( b.top() + 1 );

    c = 0;
    for( int col = b.left(); col <= b.right(); ++col )
      if( b.get_bit( b.bottom(), col ) && ++c > 1 ) break;
    if( c <= 1 ) b.bottom( b.bottom() - 1 );
    }
  }


void remove_left_right_noise( std::vector< Block * > & blockp_vector ) throw()
  {
  int blocks = blockp_vector.size();
  for( int i = 0; i < blocks; ++i )
    {
    Block & b = *blockp_vector[i];
    if( b.width() < 6 ) continue;

    int c = 0;
    for( int row = b.top(); row <= b.bottom(); ++row )
      if( b.get_bit( row, b.left() ) && ++c > 1 ) break;
    if( c <= 1 ) b.left( b.left() + 1 );

    c = 0;
    for( int row = b.top(); row <= b.bottom(); ++row )
      if( b.get_bit( row, b.right() ) && ++c > 1 ) break;
    if( c <= 1 ) b.right( b.right() - 1 );
    }
  }


void find_holes( std::vector< Zone > & zone_vector ) throw()
  {
  for( unsigned int zi = 0; zi < zone_vector.size(); ++zi )
    {
    std::vector< std::vector< Block * > > & blockp_matrix = zone_vector[zi].blockp_matrix;
    for( unsigned int bmi = 0; bmi < blockp_matrix.size(); ++bmi )
      {
      std::vector< Block * > & blockp_vector = blockp_matrix[bmi];
      for( unsigned int bvi = 0; bvi < blockp_vector.size(); ++bvi )
        blockp_vector[bvi]->find_holes();
      }
    }
  }


void scan_page( const Page_image & page_image,
                std::vector< Zone > & zone_vector, const int debug_level ) throw()
  {
  for( int zindex = 0; zindex < page_image.zones(); ++zindex )
    {
    zone_vector.push_back( Zone( page_image.rectangle( zindex ) ) );
    const Rectangle & re = zone_vector[zindex].rectangle;
    const int zthreshold = page_image.threshold( zindex );
    std::vector< std::vector< Block * > > & blockp_matrix = zone_vector[zindex].blockp_matrix;
    std::vector< Block * > blockp_vector;
    std::vector< Block * > old_data( re.width(), (Block *) 0 );
    std::vector< Block * > new_data( re.width(), (Block *) 0 );

    for( int row = re.top() + 1; row <= re.bottom(); ++row )
      {
      bool blank_row = true;
      old_data.swap( new_data );
      for( int col = re.left() + 1; col < re.right(); ++col )
        {
        const int dcol = col - re.left();
        if( !page_image.get_bit( row, col, zthreshold ) )
          new_data[dcol] = 0;			// white point
        else					// black point
          {
          blank_row = false;
          Block *p;
          Block *lp  = new_data[dcol-1];
          Block *ltp = old_data[dcol-1];
          Block *tp  = old_data[dcol];
          Block *rtp = old_data[dcol+1];
          if( lp )       { p = lp;  p->add_point( row, col ); }
          else if( ltp ) { p = ltp; p->add_point( row, col ); }
          else if( tp )  { p = tp;  p->add_point( row, col ); }
          else if( rtp ) { p = rtp; p->add_point( row, col ); }
          else
            {
            p = new Block( col, row, col, row );
            p->set_bit( row, col, true );
            blockp_vector.push_back( p );
            }
          new_data[dcol] = p;
          if( rtp && p != rtp )
            join_blocks( blockp_vector, old_data, new_data, p, rtp, dcol );
          }
        }
      if( blank_row && blockp_vector.size() )
        {
        blockp_matrix.push_back( std::vector< Block * >() );
        blockp_vector.swap( blockp_matrix.back() );
        }
      }
    if( blockp_vector.size() )
      {
      blockp_matrix.push_back( std::vector< Block * >() );
      blockp_vector.swap( blockp_matrix.back() );
      }
    }

  if( debug_level <= 99 )
    for( unsigned int zindex = 0; zindex < zone_vector.size(); ++zindex )
      {
      ignore_wide_blocks( zone_vector[zindex] );
      std::vector< std::vector< Block * > > & blockp_matrix = zone_vector[zindex].blockp_matrix;
      for( unsigned int cut = 0; cut < blockp_matrix.size(); ++cut )
        {
        std::vector< Block * > & blockp_vector = blockp_matrix[cut];
        ignore_small_blocks( blockp_vector );
        ignore_abnormal_blocks( blockp_vector );
        remove_top_bottom_noise( blockp_vector );
        remove_left_right_noise( blockp_vector );
        }
      }

  for( int zindex = zone_vector.size() - 1; zindex >= 0; --zindex )
    {
    std::vector< std::vector< Block * > > & blockp_matrix = zone_vector[zindex].blockp_matrix;
    for( int cut = blockp_matrix.size() - 1; cut >= 0; --cut )
      if( !blockp_matrix[cut].size() )
        blockp_matrix.erase( blockp_matrix.begin() + cut );
    if( !blockp_matrix.size() )
      zone_vector.erase( zone_vector.begin() + zindex );
    }

  find_holes( zone_vector );
  }

} // end namespace


Textpage::Textpage( const Page_image & page_image,
                    const char * filename, const Control & control ) throw()
  : Rectangle( page_image ), name( filename )
  {
  const int debug_level = control.debug_level;
  if( debug_level < 0 || debug_level > 100 ) return;

  std::vector< Zone > zone_vector;
  scan_page( page_image, zone_vector, debug_level );

  if( debug_level >= 98 )
    {
    if( control.outfile ) bprint( zone_vector, control.outfile );
    return;
    }
  if( debug_level > 95 || ( debug_level > 89 && debug_level < 94 ) ) return;

  // build a Textblock for every zone with text
  for( unsigned int i = 0; i < zone_vector.size(); ++i )
    {
    Textblock * tbp = new Textblock( zone_vector[i].rectangle,
                                     zone_vector[i].blockp_matrix );
    if( debug_level < 90 )
      tbp->recognize( control.charset, control.filter );
    if( tbp->textlines() ) tbpv.push_back( tbp );
    else delete tbp;
    }
  if( debug_level == 0 )
    {
    if( control.outfile ) print( control );
    if( control.exportfile ) xprint( control );
    return;
    }
  if( !control.outfile ) return;
  if( debug_level >= 86 )
    {
    bool graph = ( debug_level >= 88 );
    bool recursive = ( debug_level & 1 );
    for( int i = 0; i < textblocks(); ++i )
      tbpv[i]->dprint( control, graph, recursive );
    return;
    }
  if( debug_level > 77 ) return;
  if( debug_level >= 70 )
    {
    Page_image tmp( page_image );
    if( ( debug_level - 70 ) & 1 )	// mark zones
      {
      for( unsigned int i = 0; i < zone_vector.size(); ++i )
        tmp.draw_rectangle( zone_vector[i].rectangle );
      }
    if( ( debug_level - 70 ) & 2 )	// mark lines
      {
      for( int i = 0; i < textblocks(); ++i ) tbpv[i]->lmark( tmp );
      }
    if( ( debug_level - 70 ) & 4 )	// mark characters
      {
      for( int i = 0; i < textblocks(); ++i ) tbpv[i]->cmark( tmp );
      }
    tmp.save( control.outfile, control.filetype );
    return;
    }
  }


Textpage::~Textpage() throw()
  {
  for( int i = textblocks() - 1; i >= 0; --i ) delete tbpv[i];
  }

/*
const Textblock & Textpage::textblock( int i ) const throw()
  {
  if( i < 0 || i >= textblocks() )
    Ocrad::internal_error( "Textpage::textblock, index out of bounds" );
  return *(tbpv[i]);
  }
*/

void Textpage::print( const Control & control ) const throw()
  {
  if( control.outfile )
    for( int i = 0; i < textblocks(); ++i )
      tbpv[i]->print( control );
  }


void Textpage::xprint( const Control & control ) const throw()
  {
  if( !control.exportfile ) return;

  std::fprintf( control.exportfile, "source file %s\n", name.c_str() );
  std::fprintf( control.exportfile, "total text blocks %d\n", textblocks() );

  for( int i = 0; i < textblocks(); ++i )
    {
    const Textblock & tb = *(tbpv[i]);
    std::fprintf( control.exportfile, "text block %d %d %d %d %d\n", i + 1,
                  tb.left(), tb.top(), tb.width(), tb.height() );
    tb.xprint( control );
    }
  }
