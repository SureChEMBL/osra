/* Added Character::result( const Control & control )
   Igor Filippov, 2007 */

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
#include "ucs.h"
#include "bitmap.h"
#include "block.h"
#include "character.h"


Character::Character( const Character & c ) throw()
  : Rectangle( c ), bpv( c.bpv ), gv( c.gv )
  {
  for( unsigned int i = 0; i < bpv.size(); ++i )
    bpv[i] = new Block( *c.bpv[i] );
  }


Character & Character::operator=( const Character & c ) throw()
  {
  if( this != &c )
    {
    Rectangle::operator=( c );
    for( unsigned int i = 0; i < bpv.size(); ++i ) delete bpv[i];
    bpv = c.bpv;
    for( unsigned int i = 0; i < bpv.size(); ++i )
      bpv[i] = new Block( *c.bpv[i] );
    gv = c.gv;
    }
  return *this;
  }


Character::~Character() throw()
  {
  for( unsigned int i = 0; i < bpv.size(); ++i ) delete bpv[i];
  }


// Return the filled area of the main blocks only (no recursive)
//
int Character::area() const throw()
  {
  int a = 0;
  for( int i = 0; i < blocks(); ++i ) a += bpv[i]->area();
  return a;
  }


Block & Character::block( int i ) throw()
  {
  if( i < 0 || i >= blocks() )
    Ocrad::internal_error( "block, index out of bounds" );
  return *bpv[i];
  }


Block & Character::main_block() throw()
  {
  int imax = 0;
  for( int i = 1; i < blocks(); ++i )
    if( bpv[i]->size() > bpv[imax]->size() )
      imax = i;
  return *bpv[imax];
  }


void Character::shift_blockp( Block * p ) throw()
  {
  add_rectangle( *p );
  int i = blocks() - 1;
  for( ; i >= 0; --i )
    {
    Block & bi = *bpv[i];
    if( p->vcenter() > bi.vcenter() ) break;
    if( p->vcenter() == bi.vcenter() && p->hcenter() >= bi.hcenter() ) break;
    }
  bpv.insert( bpv.begin() + i + 1, p );
  }


void Character::add_guess( int code, int value ) throw()
  {
  gv.push_back( Guess( code, value ) );
  }


void Character::insert_guess( int i, int code, int value ) throw()
  {
  if( i < 0 || i > guesses() )
    Ocrad::internal_error( "insert_guess, index out of bounds" );
  gv.insert( gv.begin() + i, Guess( code, value ) );
  }


void Character::delete_guess( int i ) throw()
  {
  if( i < 0 || i >= guesses() )
    Ocrad::internal_error( "delete_guess, index out of bounds" );
  gv.erase( gv.begin() + i );
  }


void Character::only_guess( int code, int value ) throw()
  { gv.clear(); add_guess( code, value ); }


void Character::swap_guesses( int i, int j ) throw()
  {
  if( i < 0 || i >= guesses() || j < 0 || j >= guesses() )
    Ocrad::internal_error( "swap_guesses, index out of bounds" );
  int code = gv[i].code;
  gv[i].code = gv[j].code; gv[j].code = code;
  }


const Character::Guess & Character::guess( int i ) const throw()
  {
  if( i < 0 || i >= guesses() )
    Ocrad::internal_error( "guess, index out of bounds" );
  return gv[i];
  }


bool Character::maybe( int code ) const throw()
  {
  for( int i = 0; i < guesses(); ++i )
    if( code == gv[i].code ) return true;
  return false;
  }

/*
bool Character::maybe_digit() const throw()
  {
  for( int i = 0; i < guesses(); ++i )
    if( UCS::isdigit( gv[i].code ) ) return true;
  return false;
  }


bool Character::maybe_letter() const throw()
  {
  for( int i = 0; i < guesses(); ++i )
    if( UCS::isalpha( gv[i].code ) ) return true;
  return false;
  }
*/

void Character::join( Character & c ) throw()
  {
  for( int i = 0; i < c.blocks(); ++i ) shift_blockp( c.bpv[i] );
  c.bpv.clear();
  }


void Character::print( const Control & control ) const throw()
  {
  if( guesses() )
    {
//    if( maybe( '\t' ) ) std::putc( '\t', control.outfile ); else
    switch( control.format )
      {
      case Control::byte:
        { char ch = UCS::map_to_byte( gv[0].code );
        if( ch ) std::putc( ch, control.outfile ); }
        break;
      case Control::utf8:
        if( gv[0].code )
          std::fputs( UCS::ucs_to_utf8( gv[0].code ), control.outfile );
      }
    }
  else std::putc( '_', control.outfile );
  }

char Character::result( const Control & control ) const throw()
  {
  if( guesses() )
    {
    switch( control.format )
      {
      case Control::byte:
        { char ch = UCS::map_to_byte( gv[0].code );
        if( ch ) return(ch); }
        break;
      case Control::utf8:
        if( gv[0].code )
	  {
	    char *ch= UCS::ucs_to_utf8( gv[0].code );
	    return( ch[0]);
	  }
	break;
      default:
        { char ch = UCS::map_to_byte( gv[0].code );
	  if( ch ) return(ch); }
        break;
      }
    }
  return('_');
  }


void Character::dprint( const Control & control, const Rectangle & charbox,
                        bool graph, bool recursive ) const throw()
  {
  if( graph || recursive )
    std::fprintf( control.outfile, "%d guesses    ", guesses() );
  for( int i = 0; i < guesses(); ++i )
    {
    switch( control.format )
      {
      case Control::byte:
        { char ch = UCS::map_to_byte( gv[i].code );
        if( ch ) std::fprintf( control.outfile, "guess '%c', confidence %d    ",
                               ch, gv[i].value ); }
        break;
      case Control::utf8:
        std::fprintf( control.outfile, "guess '%s', confidence %d    ",
                      UCS::ucs_to_utf8( gv[i].code ), gv[i].value );
      }
    if( !graph && !recursive ) break;
    }
  std::fputs( "\n", control.outfile );
  if( graph )
    {
    std::fprintf( control.outfile,
                  "left = %d, top = %d, right = %d, bottom = %d\n",
                  left(), top(), right(), bottom() );
    std::fprintf( control.outfile,
                  "width = %d, height = %d, hcenter = %d, vcenter = %d, black area = %d%%\n\n",
                  width(), height(), hcenter(), vcenter(), ( area() * 100 ) / size() );

    const int minrow = std::min( top(), charbox.top() );
    const int maxrow = std::max( bottom(), charbox.bottom() );
    for( int row = minrow; row <= maxrow; ++row )
      {
      bool istop = ( row == top() );
      bool isvc = ( row == vcenter() );
      bool isbot = ( row == bottom() );
      bool iscbtop = ( row == charbox.top() );
      bool iscbvc = ( row == charbox.vcenter() );
      bool iscbbot = ( row == charbox.bottom() );

      bool ish1top = false, ish1bot = false, ish2top = false, ish2bot = false;
      if( blocks() == 1 && bpv[0]->holes() )
        {
        const Block & b = *bpv[0];
        ish1top = ( row == b.hole(0).top() );
        ish1bot = ( row == b.hole(0).bottom() );
        if( b.holes() > 1 )
          {
          ish2top = ( row == b.hole(1).top() );
          ish2bot = ( row == b.hole(1).bottom() );
          }
        }

      for( int col = left(); col <= right(); ++col )
        {
        char ch = ( isvc && col == hcenter() ) ? '+' : '·';
        for( int i = 0; i < blocks(); ++i )
          {
          int id = bpv[i]->id( row, col );
          if( id != 0 )
            {
            if( id > 0 ) ch = (ch == '+') ? 'C' : 'O';
            else ch = (ch == '+') ? '=' : '-';
            break;
            }
          }
        std::fprintf( control.outfile, " %c", ch );
        }
      if( istop ) std::fprintf( control.outfile, "  top(%d)", row );
      if( isvc ) std::fprintf( control.outfile, "  vcenter(%d)", row );
      if( isbot ) std::fprintf( control.outfile, "  bottom(%d)", row );

      if( iscbtop ) std::fprintf( control.outfile, "  box.top(%d)", row );
      if( iscbvc ) std::fprintf( control.outfile, "  box.vcenter(%d)", row );
      if( iscbbot ) std::fprintf( control.outfile, "  box.bottom(%d)", row );

      if( ish1top ) std::fprintf( control.outfile, "  h1.top(%d)", row );
      if( ish1bot ) std::fprintf( control.outfile, "  h1.bottom(%d)", row );
      if( ish2top ) std::fprintf( control.outfile, "  h2.top(%d)", row );
      if( ish2bot ) std::fprintf( control.outfile, "  h2.bottom(%d)", row );

      std::fputs( "\n", control.outfile );
      }
    std::fputs( "\n\n", control.outfile );
    }
  }


void Character::xprint( const Control & control ) const throw()
  {
  std::fprintf( control.exportfile, "%3d %3d %2d %2d; %d",
                left(), top(), width(), height(), guesses() );

  for( int i = 0; i < guesses(); ++i )
    switch( control.format )
      {
      case Control::byte:
        { char ch = UCS::map_to_byte( gv[i].code );
        if( !ch ) ch = '_';
        std::fprintf( control.exportfile, ", '%c'%d", ch, gv[i].value ); }
        break;
      case Control::utf8:
        std::fprintf( control.exportfile, ", '%s'%d",
                      UCS::ucs_to_utf8( gv[i].code ), gv[i].value );
      }
  std::fputs( "\n", control.exportfile );
  }


void Character::apply_filter( const Filter & filter ) throw()
  {
  if( filter.type() == Filter::none ) return;
  const int code = guesses() ? gv[0].code : 0;
  bool remove = false;

  switch( filter.type() )
    {
    case Filter::none:			// only for completeness
      break;
    case Filter::letters_only:
      remove = true;
    case Filter::letters:
      if( !UCS::isalpha( code ) && !UCS::isspace( code ) )
        {
        for( int i = 1; i < guesses(); ++i )
          if( UCS::isalpha( gv[i].code ) ) { swap_guesses( 0, i ); break; }
        if( guesses() && !UCS::isalpha( gv[0].code ) )
          gv[0].code = UCS::to_nearest_letter( gv[0].code );
        if( remove && ( !guesses() || !UCS::isalpha( gv[0].code ) ) )
          only_guess( 0, 0 );
        }
      break;
    case Filter::numbers_only:
      remove = true;
    case Filter::numbers:
      if( !UCS::isdigit( code ) && !UCS::isspace( code ) )
        {
        for( int i = 1; i < guesses(); ++i )
          if( UCS::isdigit( gv[i].code ) ) { swap_guesses( 0, i ); break; }
        if( guesses() && !UCS::isdigit( gv[0].code ) )
          gv[0].code = UCS::to_nearest_digit( gv[0].code );
        if( remove && ( !guesses() || !UCS::isdigit( gv[0].code ) ) )
          only_guess( 0, 0 );
        }
      break;
    }
  }
