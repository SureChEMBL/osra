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
#include "bitmap.h"
#include "block.h"
#include "character.h"
#include "page_image.h"
#include "textline.h"
#include "textblock.h"


namespace {

void add_line( std::vector< Textline * > & textlinep_vector ) throw()
  { textlinep_vector.push_back( new Textline() ); }


void insert_line( std::vector< Textline * > & textlinep_vector, int i ) throw()
  { textlinep_vector.insert( textlinep_vector.begin() + i, new Textline() ); }


void delete_line( std::vector< Textline * > & textlinep_vector, int i ) throw()
  {
  delete textlinep_vector[i];
  textlinep_vector.erase( textlinep_vector.begin() + i );
  }


// Build the vertical composite characters.
//
void join_characters( std::vector< Textline * > & tlpv ) throw()
  {
  for( unsigned int current_line = 0; current_line < tlpv.size(); ++current_line )
    {
    Textline & line = *tlpv[current_line];
    for( int i = 0 ; i < line.characters() - 1; )
      {
      Character & c1 = line.character( i );
      bool joined = false;
      for( int j = i + 1 ; j < line.characters(); ++j )
        {
        Character & c2 = line.character( j );
        if( !c1.h_overlaps( c2 ) ) continue;
        Character *cup, *cdn;
        if( c1.vcenter() < c2.vcenter() ) cup = &c1, cdn = &c2;
        else cup = &c2, cdn = &c1;
        if( cdn->includes_hcenter( *cup ) || cup->includes_hcenter( *cdn ) ||
            ( cdn->top() > cup->bottom() && cdn->hcenter() < cup->hcenter() ) ||
            ( cdn->blocks() == 2 &&
              2 * cdn->block( 0 ).size() < cdn->block( 1 ).size() &&
              cdn->block( 0 ).includes_vcenter( *cup ) ) )
          {
          int k;
          if( 64 * c1.size() < c2.main_block().size() ) k = i;
          else if( 64 * c2.size() < c1.main_block().size() ) k = j;
          else if( cdn == &c2 ) { c2.join( c1 ); k = i; }
          else { c1.join( c2 ); k = j; }
          line.delete_character( k );
          joined = true; break;
          }
        }
      if( !joined ) ++i;
      }
    }
  }


// Return the character position >= first preceding a big gap or eol.
//
int find_big_gap( const Textline & line, const int first, const int space_width_limit ) throw()
  {
  int i = first;
  while( i + 1 < line.characters() )
    {
    Character & c1 = line.character( i );
    Character & c2 = line.character( i + 1 );
    const int gap = c2.left() - c1.right() - 1;
    if( gap > space_width_limit ) break; else ++i;
    }
  return i;
  }


// Insert spaces between characters.
//
void insert_spaces( std::vector< Textline * > & tlpv ) throw()
  {
  for( unsigned int current_line = 0; current_line < tlpv.size(); ++current_line )
    {
    Textline & line1 = *tlpv[current_line];
    const Rational mw = line1.mean_width();
    if( mw < 2 ) continue;
    const int space_width_limit = ( 3 * mw ).trunc();
    int first = 0;
    while( first + 1 < line1.characters() )
      {
      int last = find_big_gap( line1, first, space_width_limit );
      const Rational mg = line1.mean_gap_width( first, last );
      if( first < last && mg >= 0 )
        {
        int spaces = 0, nospaces = 0, spsum = 0, nospsum = 0;
        for( int i = first ; i < last; ++i )
          {
          Character & c1 = line1.character( i );
          Character & c2 = line1.character( i + 1 );
          const int gap = c2.left() - c1.right() - 1;
          if( gap >= mw.trunc() || gap > 3 * mg ||
              ( 5 * gap > 2 * mw && gap > 2 * mg ) ||
              ( 3 * c1.width() > 2 * mw && 3 * c2.width() > 2 * mw && 2 * gap > mw && 5 * gap > 8 * mg ) )
            { ++spaces; spsum += gap;
              if( line1.insert_space( i + 1 ) ) { ++i; ++last; } }
          else { ++nospaces; nospsum += gap; }
          }
        if( spaces && nospaces )
          {
          const Rational th = ( Rational( 3 * spsum, spaces ) + Rational( nospsum, nospaces ) ) / 4;
          for( int i = first ; i < last; ++i )
            {
            Character & c1 = line1.character( i );
            Character & c2 = line1.character( i + 1 );
            const int gap = c2.left() - c1.right() - 1;
            if( gap > th && line1.insert_space( i + 1 ) ) { ++i; ++last; }
            }
          }
        }
      if( ++last < line1.characters() && line1.insert_space( last, true ) )
        ++last;
      first = last;
      }
    }
  }

} // end namespace


Textblock::Textblock( const Rectangle & r,
                      std::vector< std::vector< Block * > > & blockp_matrix ) throw()
  : Rectangle( r )
  {
  int cuts = blockp_matrix.size();
  std::vector< int > mean_height( cuts, 0 );
  std::vector< std::vector< Block * > > pending( cuts );
  std::vector< std::vector< Block * > > pending_tall( cuts );
  std::vector< std::vector< Block * > > pending_short( cuts );

  // Classify blocks by height.
  for( int cut = 0; cut < cuts; ++cut )
    {
    const std::vector< Block * > & blockp_vector = blockp_matrix[cut];
    if( !blockp_vector.size() ) continue;
    unsigned int samples = 0;
    std::vector< int > height_distrib;
    for( unsigned int i = 0; i < blockp_vector.size(); ++i )
      {
      unsigned int h = blockp_vector[i]->height();
      unsigned int w = blockp_vector[i]->width();
      if( h < 10 || w >= 3 * h ) continue;
      if( h >= height_distrib.size() ) height_distrib.resize( h + 1 );
      ++height_distrib[h]; ++samples;
      }
    if( !height_distrib.size() )
      for( unsigned int i = 0; i < blockp_vector.size(); ++i )
        {
        unsigned int h = blockp_vector[i]->height();
        if( h >= height_distrib.size() ) height_distrib.resize( h + 1 );
        ++height_distrib[h]; ++samples;
        }
    int valid_samples = 0;
    for( unsigned int i = 0, count = 0; i < height_distrib.size(); ++i )
      {
      int a = height_distrib[i];
      if( 10 * ( count + a ) >= samples && 10 * count < 9 * samples )
        { mean_height[cut] += a * i; valid_samples += a; }
      count += a;
      }
    if( valid_samples ) mean_height[cut] /= valid_samples;

    for( unsigned int i = 0; i < blockp_vector.size(); ++i )
      {
      Block * p = blockp_vector[i];
      if( p->height() >= 2 * mean_height[cut] )
        pending_tall[cut].push_back( p );
      else if( 2 * p->height() <= mean_height[cut] )
        pending_short[cut].push_back( p );
      else pending[cut].push_back( p );
      }
    }

  // Assign normal blocks to characters and create lines.
  int min_line = 0;	// first line of current cut
  add_line( tlpv );
  for( int cut = 0; cut < cuts; ++cut )
    {
    if( pending[cut].size() )
      {
      if( tlpv.back()->characters() ) add_line( tlpv );
      int current_line = min_line = textlines() - 1;
      tlpv[current_line]->shift_characterp( new Character( pending[cut][0] ) );
      for( unsigned int i = 1; i < pending[cut].size(); ++i )
        {
        Block & b = *pending[cut][i];
        current_line = std::max( min_line, current_line - 2 );
        while( true )
          {
          const Character *cl = 0, *cr = 0;
          for( int j = tlpv[current_line]->characters() - 1; j >= 0; --j )
            {
            const Character & cj = tlpv[current_line]->character( j );
            if( !b.includes_hcenter( cj ) && !cj.includes_hcenter( b ) )
              { if( b.h_precedes( cj ) ) cr = &cj; else { cl = &cj; break; } }
            }
          if( ( cl && ( cl->includes_vcenter( b ) || b.includes_vcenter( *cl ) ) ) ||
              ( cr && ( cr->includes_vcenter( b ) || b.includes_vcenter( *cr ) ) ) )
            { tlpv[current_line]->shift_characterp( new Character( &b ) ); break; }
          else if( ( cl && cl->bottom() < b.top() ) || ( cr && cr->bottom() < b.top() ) )
            {
            if( ++current_line >= textlines() )
              { add_line( tlpv ); current_line = textlines() - 1;
              tlpv[current_line]->shift_characterp( new Character( &b ) ); break; }
            }
          else if( ( cl && cl->top() > b.bottom() ) || ( cr && cr->top() > b.bottom() ) )
            {
            insert_line( tlpv, current_line );
            tlpv[current_line]->shift_characterp( new Character( &b ) ); break;
            }
          else if( ( cl && cl->v_overlaps( b ) ) || ( cr && cr->v_overlaps( b ) ) )
            { tlpv[current_line]->shift_characterp( new Character( &b ) ); break; }
          else { delete &b; break; }
          }
        }
      }
    }

  join_characters( tlpv );

  // Create tracks of lines.
  for( int i = 0; i < textlines(); ++i ) tlpv[i]->set_track();

  // Insert tall blocks.
  // Seek up, then seek down, needed for slanted or curved lines.
  for( int current_line = 0, cut = 0; cut < cuts; ++cut )
    {
    for( unsigned int i = 0; i < pending_tall[cut].size(); ++i )
      {
      Block & b = *pending_tall[cut][i];
      while( current_line > 0 &&
             b.bottom() < tlpv[current_line]->vcenter( b.hcenter() ) )
        --current_line;
      while( current_line < textlines() &&
             b.top() > tlpv[current_line]->vcenter( b.hcenter() ) )
        ++current_line;
      if( current_line >= textlines() )
        { --current_line; delete &b; continue; }
      Textline & l = *tlpv[current_line];
      if( b.height() <= 3 * l.mean_height() &&
          ( b.height() <= 2 * l.mean_height() ||
            ( l.big_initial() && l.big_initial()->left() < b.left() ) ||
            ( !l.big_initial() && l.character(0).left() < b.left() ) ) )
        l.shift_characterp( new Character( &b ) );
      else		// possible big initial
        {
        if( l.big_initial() && l.big_initial()->height() <= 2 * l.mean_height() )
          {
          l.shift_characterp( new Character( *l.big_initial() ) );
          l.big_initial( 0 );
          }
        if( ( !l.big_initial() || l.big_initial()->left() > b.left() ) &&
            l.character(0).left() > b.hcenter() )
          l.big_initial( new Character( &b ) );
        }
      }
    }

  // Insert short blocks.
  // Seek up, then seek down, needed for slanted or curved lines.
  for( int current_line = 0, cut = 0; cut < cuts; ++cut )
    {
    for( unsigned int i = 0; i < pending_short[cut].size(); ++i )
      {
      Block & b = *pending_short[cut][i];
      while( current_line > 0 &&
             b.bottom() < tlpv[current_line]->top( b.hcenter() ) )
        --current_line;
      int temp = current_line;
      while( current_line < textlines() &&
             b.top() > tlpv[current_line]->bottom( b.hcenter() ) )
        ++current_line;
      if( current_line >= textlines() )
        {
        const Textline & l = *tlpv[--current_line];
        if( b.top() > l.bottom( b.hcenter() ) + l.height() ) continue;
        else temp = current_line;
        }
      if( current_line - temp > 1 ) temp = current_line - 1;
      if( current_line != temp &&
          2 * ( b.top() - tlpv[temp]->bottom( b.hcenter() ) ) <
          tlpv[current_line]->top( b.hcenter() ) - b.bottom() )
        current_line = temp;
      tlpv[current_line]->shift_characterp( new Character( &b ) );
      }
    }

  // remove clipped lines at top or bottom of text block
  if( textlines() )
    {
    const Textline & l = *tlpv[textlines()-1];
    for( int i = 0, c = 0; i < l.characters(); ++i )
      if( l.character( i ).bottom() >= bottom() && 2 * ++c >= l.characters() )
        { delete_line( tlpv, textlines() - 1 ); break; }
    }
  if( textlines() )
    {
    const Textline & l = *tlpv[0];
    const int t = std::max( top(), 1 );
    for( int i = 0, c = 0; i < l.characters(); ++i )
      if( l.character( i ).top() <= t && 2 * ++c >= l.characters() )
        { delete_line( tlpv, 0 ); break; }
    }

  // Second pass. Join lines of i-dots and tildes.
  for( int current_line = 0; current_line < textlines() - 1; )
    {
    bool joined = false;
    Textline & line1 = *tlpv[current_line];
    Textline & line2 = *tlpv[current_line+1];
    if( line1.characters() <= 2 * line2.characters() &&
        2 * line1.mean_height() < line2.mean_height() )
      for( int i1 = 0; !joined && i1 < line1.characters(); ++i1 )
        {
        Character & c1 = line1.character( i1 );
        if( 2 * c1.height() >= line2.mean_height() ) continue;
        for( int i2 = 0; !joined && i2 < line2.characters(); ++i2 )
          {
          Character & c2 = line2.character( i2 );
          if( c2.right() < c1.left() ) continue;
          if( c2.left() > c1.right() ) break;
          if( ( c2.includes_hcenter( c1 ) || c1.includes_hcenter( c2 ) )
              && c2.top() - c1.bottom() < line2.mean_height() )
            {
            joined = true; line2.join( line1 );
            delete_line( tlpv, current_line );
            }
          }
        }
    if( !joined ) ++current_line;
    }

  join_characters( tlpv );

  // Fourth pass. Remove noise lines.
  if( textlines() >= 3 )
    {
    for( int i = 0; i + 2 < textlines(); ++i )
      {
      Textline & line1 = *tlpv[i];
      Textline & line2 = *tlpv[i+1];
      Textline & line3 = *tlpv[i+2];
      if( line2.characters() > 2 || line1.characters() < 4 ||
          line3.characters() < 4 ) continue;
      if( !Ocrad::similar( line1.height(), line3.height(), 10 ) ) continue;
      if( 8 * line2.height() > line1.height() + line3.height() ) continue;
      delete_line( tlpv, i + 1 );
      }
    }

  // Remove leading and trailing noise characters.
  for( int i = 0; i < textlines(); ++i )
    {
    Textline & l = *tlpv[i];
    if( !l.big_initial() && l.characters() > 2 )
      {
      const Character & c0 = l.character( 0 );
      const Character & c1 = l.character( 1 );
      const Character & c2 = l.character( 2 );
      if( c0.blocks() == 1 &&
          4 * c0.size() < c1.size() && c1.left() - c0.right() > 2 * l.height() &&
          4 * c0.size() < c2.size() && c2.left() - c1.right() < l.height() )
        l.delete_character( 0 );
      }
    if( l.characters() > 2 )
      {
      const Character & c0 = l.character( l.characters() - 1 );
      const Character & c1 = l.character( l.characters() - 2 );
      const Character & c2 = l.character( l.characters() - 3 );
      if( c0.blocks() == 1 &&
          4 * c0.size() < c1.size() && c0.left() - c1.right() > 2 * l.height() &&
          4 * c0.size() < c2.size() && c1.left() - c2.right() < l.height() )
        l.delete_character( l.characters() - 1 );
      }
    }
  insert_spaces( tlpv );
  }


Textblock::~Textblock() throw()
  {
  for( int i = textlines() - 1; i >= 0; --i ) delete tlpv[i];
  }


void Textblock::recognize( const Charset & charset, const Filter & filter ) throw()
  {
  // Recognize characters.
  for( int i = 0; i < textlines(); ++i )
    {
    // First pass. Recognize the easy characters.
    tlpv[i]->recognize1( charset );
    // Second pass. Use context to clear up ambiguities.
    tlpv[i]->recognize2( charset );
    }

  if( filter.type() != Filter::none )
    for( int i = 0; i < textlines(); ++i )
      tlpv[i]->apply_filter( filter );

  // Remove unrecognized lines.
  for( int i = textlines() - 1; i >= 0; --i )
    {
    Textline & line1 = *tlpv[i];
    bool flag = false;
    if( line1.big_initial() && line1.big_initial()->guesses() ) flag = true;
    else
      for( int j = 0 ; j < line1.characters(); ++j )
        { if( line1.character( j ).guesses() ) { flag = true; break; } }
    if( !flag ) delete_line( tlpv, i );
    }

  // Add blank lines.
  if( textlines() >= 3 )
    {
    int min_vdistance = ( tlpv.back()->mean_vcenter() - tlpv.front()->mean_vcenter() ) / ( textlines() - 1 );
    for( int i = 0; i + 1 < textlines(); ++i )
      {
      const Textline & line1 = *tlpv[i];
      const Textline & line2 = *tlpv[i+1];
      if( !Ocrad::similar( line1.characters(), line2.characters(), 50 ) ||
          !Ocrad::similar( line1.width(), line2.width(), 30 ) ) continue;
      const int vdistance = line2.mean_vcenter() - line1.mean_vcenter();
      if( vdistance >= min_vdistance ) continue;
      const int mh1 = line1.mean_height(), mh2 = line2.mean_height();
      if( mh1 < 10 || mh2 < 10 ) continue;
      if( Ocrad::similar( mh1, mh2, 20 ) && 2 * vdistance > mh1 + mh2 )
        min_vdistance = vdistance;
      }
    if( min_vdistance > 0 )
      for( int i = 0; i + 1 < textlines(); ++i )
        {
        const Textline & line1 = *tlpv[i];
        const Textline & line2 = *tlpv[i+1];
        int vdistance = line2.mean_vcenter() - line1.mean_vcenter() - min_vdistance;
        while( 2 * vdistance > min_vdistance )
          { insert_line( tlpv, ++i ); vdistance -= min_vdistance; }
        }
    }
  }


int Textblock::characters() const throw()
  {
  int total = 0;
  for( int i = 0; i < textlines(); ++i )
    total += tlpv[i]->characters();
  return total;
  }

/*
Textline & Textblock::textline( int i ) const throw()
  {
  if( i < 0 || i >= textlines() )
    Ocrad::internal_error( "line, index out of bounds" );
  return *tlpv[i];
  }
*/

void Textblock::print( const Control & control ) const throw()
  {
  for( int i = 0; i < textlines(); ++i )
    tlpv[i]->print( control );
  std::fputs( "\n", control.outfile );
  }


void Textblock::dprint( const Control & control, bool graph, bool recursive )
								const throw()
  {
  std::fprintf( control.outfile, "%d lines\n\n", textlines() );

  for( int i = 0; i < textlines(); ++i )
    {
    std::fprintf( control.outfile, "%d characters in line %d\n",
                  tlpv[i]->characters(), i + 1 );
    tlpv[i]->dprint( control, graph, recursive );
    }
  std::fputs( "\n", control.outfile );
  }


void Textblock::xprint( const Control & control ) const throw()
  {
  std::fprintf( control.exportfile, "lines %d\n", textlines() );

  for( int i = 0; i < textlines(); ++i )
    {
    std::fprintf( control.exportfile, "line %d chars %d height %d\n", i + 1,
                  tlpv[i]->characters(), tlpv[i]->mean_height() );
    tlpv[i]->xprint( control );
    }
  }


void Textblock::cmark( Page_image & page_image ) const throw()
  {
  for( int i = 0; i < textlines(); ++i ) tlpv[i]->cmark( page_image );
  }


void Textblock::lmark( Page_image & page_image ) const throw()
  {
  for( int i = 0; i < textlines(); ++i ) page_image.draw_track( *tlpv[i] );
  }
