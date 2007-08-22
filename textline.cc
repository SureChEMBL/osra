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
#include "track.h"
#include "ucs.h"
#include "bitmap.h"
#include "block.h"
#include "character.h"
#include "page_image.h"
#include "textline.h"


Textline::Textline( const Textline & tl ) throw()
  : Track( tl ), _big_initial( 0 )
  {
  if( tl._big_initial ) big_initial( new Character( *tl._big_initial ) );
  for( unsigned int i = 0; i < tl.cpv.size(); ++i )
    shift_characterp( new Character( *tl.cpv[i] ) );
  }


Textline & Textline::operator=( const Textline & tl ) throw()
  {
  if( this != &tl )
    {
    Track::operator=( tl );
    if( tl._big_initial ) big_initial( new Character( *tl._big_initial ) );
    else big_initial( 0 );
    for( unsigned int i = 0; i < cpv.size(); ++i ) delete cpv[i];
    cpv.clear();
    for( unsigned int i = 0; i < tl.cpv.size(); ++i )
      shift_characterp( new Character( *tl.cpv[i] ) );
    }
  return *this;
  }


Textline::~Textline() throw()
  {
  big_initial( 0 );
  for( unsigned int i = 0; i < cpv.size(); ++i ) delete cpv[i];
  }


void Textline::set_track() throw()
  {
  std::vector< Rectangle > rv;
  for( unsigned int i = 0; i < cpv.size(); ++i )
    if( !cpv[i]->maybe(' ') ) rv.push_back( *cpv[i] );
  Track::set_track( rv );
  }


void Textline::big_initial( Character * p ) throw()
  {
  if( _big_initial ) delete _big_initial;
  _big_initial = p;
  }

/*
void Textline::delete_big_initial() throw()
  {
  if( _big_initial ) { delete _big_initial; _big_initial = 0; }
  }
*/

void Textline::delete_character( int i ) throw()
  {
  if( i < 0 || i >= characters() )
    Ocrad::internal_error( "delete_character, index out of bounds" );
  delete cpv[i]; cpv.erase( cpv.begin() + i );
  }


int Textline::shift_characterp( Character * p ) throw()
  {
  int i = characters();

  while( i > 0 && p->h_precedes( *cpv[i-1] ) ) --i;
  cpv.insert( cpv.begin() + i, p );
  return i;
  }


bool Textline::insert_space( int i, bool tab ) throw()
  {
  if( i <= 0 || i >= characters() )
    Ocrad::internal_error( "insert_space, index out of bounds" );
  if( height() == 0 )
    Ocrad::internal_error( "insert_space, track not set yet" );
  Character & c1 = *cpv[i-1];
  Character & c2 = *cpv[i];
  int l = c1.right() + 1;
  int r = c2.left() - 1;
  if( l > r ) return false;
  int t = top( ( l + r ) / 2 );
  int b = bottom( ( l + r ) / 2 );
  Rectangle re( l, t, r, b );
  Character * p = new Character( re, ' ', tab ? 1 : 0 );
  if( tab ) p->add_guess( '\t', 0 );
  cpv.insert( cpv.begin() + i, p );
  return true;
  }


void Textline::join( Textline & tl ) throw()
  {
  if( tl._big_initial &&
      ( !_big_initial || tl._big_initial->left() < _big_initial->left() ) )
    { big_initial( tl._big_initial ); tl._big_initial = 0; }

  for( int i = 0; i < tl.characters(); ++i ) shift_characterp( tl.cpv[i] );

  tl.big_initial( 0 ); tl.cpv.clear();
  }


Character & Textline::character( int i ) const throw()
  {
  if( i < 0 || i >= characters() )
    Ocrad::internal_error( "character, index out of bounds" );
  return *cpv[i];
  }


Rectangle Textline::charbox( const Character & c ) const throw()
  {
  return Rectangle( c.left(), top( c.hcenter() ), c.right(), bottom( c.hcenter() ) );
  }


int Textline::mean_height() const throw()
  {
  int c = 0, sum = 0;

  for( int i = 0; i < characters(); ++i )
    if( !cpv[i]->maybe(' ') ) { ++c; sum += cpv[i]->height(); }
  if( c ) sum /= c;
  return sum;
  }


Rational Textline::mean_width() const throw()
  {
  int c = 0, sum = 0;

  for( int i = 0; i < characters(); ++i )
    if( !cpv[i]->maybe(' ') ) { ++c; sum += cpv[i]->width(); }

  if( c ) return Rational( sum, c );
  return Rational( 0 );
  }


Rational Textline::mean_gap_width( const int first, int last ) const throw()
  {
  if( last < 0 ) last = characters() - 1;
  int sum = 0;

  for( int i = first; i < last; ++i )
    sum += std::max( 0, cpv[i+1]->left() - cpv[i]->right() - 1 );

  if( last > first ) return Rational( sum, last - first );
  return Rational( 0 );
  }


int Textline::mean_hcenter() const throw()
  {
  int sum = 0;

  for( int i = 0; i < characters(); ++i ) sum += cpv[i]->hcenter();
  if( characters() ) sum /= characters();
  return sum;
  }


int Textline::mean_vcenter() const throw()
  {
  int sum = 0;

  for( int i = 0; i < characters(); ++i ) sum += cpv[i]->vcenter();
  if( characters() ) sum /= characters();
  return sum;
  }


void Textline::print( const Control & control ) const throw()
  {
  if( _big_initial ) _big_initial->print( control );
  for( int i = 0; i < characters(); ++i ) character( i ).print( control );
  std::fputs( "\n", control.outfile );
  }


void Textline::dprint( const Control & control, bool graph, bool recursive )
								const throw()
  {
  if( graph || recursive )
    std::fprintf( control.outfile, "mean height = %d, track segments = %d\n",
                  mean_height(), segments() );

  if( _big_initial )
    {
    Character & c = *_big_initial;
    c.dprint( control, c, graph, recursive );
    }
  for( int i = 0; i < characters(); ++i )
    {
    Character & c = character( i );
    c.dprint( control, charbox( c ), graph, recursive );
    }
  std::fputs( "\n", control.outfile );
  }


void Textline::xprint( const Control & control ) const throw()
  {
  if( _big_initial ) _big_initial->xprint( control );
  for( int i = 0; i < characters(); ++i )
    character( i ).xprint( control );
  }


void Textline::cmark( Page_image & page_image ) const throw()
  {
  if( _big_initial ) page_image.draw_rectangle( *_big_initial );
  for( int i = 0; i < characters(); ++i )
    page_image.draw_rectangle( character(i) );
  }


void Textline::recognize1( const Charset & charset ) const throw()
  {
  if( _big_initial )
    {
    Character & c = *_big_initial;
    c.recognize1( charset, c );
    if( c.guesses() == 1 )
      {
      int code = c.guess( 0 ).code;
      if( UCS::islower( code ) ) c.only_guess( UCS::toupper( code ), 0 );
      }
    else if( c.guesses() > 1 ) c.clear_guesses();
    }

  for( int i = 0; i < characters(); ++i )
    {
    Character & c = character( i );
    c.recognize1( charset, charbox( c ) );
    }
  }


void Textline::apply_filter( const Filter & filter ) throw()
  {
  if( _big_initial )
    {
    Character & c = *_big_initial;
    const int guesses = c.guesses();
    c.apply_filter( filter );
    if( guesses && !c.guesses() ) big_initial( 0 );
    }

  bool flag = false;
  for( int i = 0; i < characters(); )
    {
    Character & c = character( i );
    const int guesses = c.guesses();
    c.apply_filter( filter );
    if( guesses && !c.guesses() ) { delete_character( i ); flag = true; }
    else ++i;
    }
  if( flag )
    for( int i = characters() - 1; i >= 0; --i )
      if( character(i).maybe(' ') &&
          ( i >= characters() - 1 || ( i > 0 && character(i-1).maybe(' ') ) ) )
        delete_character( i );
  }
