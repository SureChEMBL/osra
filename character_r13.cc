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
#include "profile.h"
#include "feats.h"


// Recognizes 3 block characters.
// %ÄËÏÖÜäëïöüÿ÷
//
void Character::recognize13( const Charset & charset, const Rectangle & charbox ) throw()
  {
  const Block & b1 = block( 0 );
  const Block & b2 = block( 1 );
  const Block & b3 = block( 2 );		// lower block
  Character c( new Block( b3 ) );
  int code = 0;

  c.recognize1( charset, charbox );
  if( c.guesses() )
    {
    if( c.maybe('.') ||
        ( c.height() < 2 * c.width() && c.maybe(',') && 2 * b3.area() >= b3.size() ) )
      {
      if( b1.bottom() <= b2.top() && b2.bottom() <= b3.top() )
        { if( b2.width() >= 2 * b2.height() ) code = UCS::DIV; }
      else code = '%';
      }
    else if( std::max( b1.width(), b2.width() ) < b3.width() &&
             Ocrad::similar( b1.height(), b2.height(), 20, 2 ) &&
             2 * std::max( b1.height(), b2.height() ) < b3.height() )
      code = UCS::compose( c.guess( 0 ).code, ':' );
    else if( c.maybe('o') )
      {
      if( ( b1.hcenter() < b2.hcenter() && b1.holes() == 1 && !b2.holes() ) ||
          ( b2.hcenter() < b1.hcenter() && b2.holes() == 1 && !b1.holes() ) )
        code = '%';
      }
    }
  if( charset.only( Charset::ascii ) )
    {
    if( code == UCS::DIV ) code = '%';
    else code = UCS::base_letter( code );
    }
  if( code ) add_guess( code, 0 );
  }
