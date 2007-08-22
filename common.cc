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
#include <cstdlib>
#include "common.h"


namespace {

const int charsets = 3;

const Charset::Value charset_value[charsets] =
  { Charset::ascii, Charset::iso_8859_9, Charset::iso_8859_15 };

const char * const charset_name[charsets] =
  { "ascii", "iso-8859-9", "iso-8859-15" };

struct F_entry
  {
  const char * name;
  Filter::Type type;
  };

const F_entry F_table[] =
  {
  { "none",         Filter::none },
  { "letters",      Filter::letters },
  { "letters_only", Filter::letters_only },
  { "numbers",      Filter::numbers },
  { "numbers_only", Filter::numbers_only },
  { 0, Filter::none }
  };

struct T_entry
  {
  const char * name;
  Transformation::Type type;
  };

const T_entry T_table[] =
  {
  { "none",      Transformation::none },
  { "rotate90",  Transformation::rotate90 },
  { "rotate180", Transformation::rotate180 },
  { "rotate270", Transformation::rotate270 },
  { "mirror_lr", Transformation::mirror_lr },
  { "mirror_tb", Transformation::mirror_tb },
  { "mirror_d1", Transformation::mirror_d1 },
  { "mirror_d2", Transformation::mirror_d2 },
  { 0, Transformation::none }
  };

} // end namespace


bool Ocrad::similar( int a, int b, int percent_dif, int abs_dif ) throw()
  {
  int difference = std::abs( a - b );
  if( percent_dif > 0 && difference <= abs_dif ) return true;
  int max_abs = std::max( std::abs(a), std::abs(b) );
  return ( difference * 100 <= max_abs * percent_dif );
  }


bool Charset::enable( const char * name ) throw()
  {
  for( int i = 0; i < charsets; ++i )
    if( std::strcmp( name, charset_name[i] ) == 0 )
      { _charset |= charset_value[i]; return true; }
  return false;
  }


bool Charset::enabled( Value cset ) const throw()
  {
  if( !_charset ) return cset == iso_8859_15;		// default charset
  return _charset & cset;
  }


bool Charset::only( Value cset ) const throw()
  {
  if( !_charset ) return cset == iso_8859_15;		// default charset
  return _charset == cset;
  }


void Charset::show_error( const char * program_name, const char * arg ) const throw()
  {
  if( arg && std::strcmp( arg, "help" ) )
    std::fprintf( stderr,"%s: bad charset `%s'\n", program_name, arg );
  std::fputs( "Valid charset names are:", stderr );
  for( int i = 0; i < charsets; ++i )
    std::fprintf( stderr, "  %s", charset_name[i] );
  std::fputs( "\n", stderr );
  }


bool Filter::set( const char * name ) throw()
  {
  for( int i = 0; F_table[i].name != 0; ++i )
    if( std::strcmp( name, F_table[i].name ) == 0 )
      { _type = F_table[i].type; return true; }
  return false;
  }


void Filter::show_error( const char * program_name, const char * arg ) const throw()
  {
  if( arg && std::strcmp( arg, "help" ) )
    std::fprintf( stderr,"%s: bad filter `%s'\n", program_name, arg );
  std::fputs( "Valid filter names are:", stderr );
  for( int i = 0; F_table[i].name != 0; ++i )
    std::fprintf( stderr, "  %s", F_table[i].name );
  std::fputs( "\n", stderr );
  }


bool Transformation::set( const char * name ) throw()
  {
  for( int i = 0; T_table[i].name != 0; ++i )
    if( std::strcmp( name, T_table[i].name ) == 0 )
      { _type = T_table[i].type; return true; }
  return false;
  }


void Transformation::show_error( const char * program_name, const char * arg ) const throw()
  {
  if( arg && std::strcmp( arg, "help" ) )
    std::fprintf( stderr,"%s: bad bitmap trasformation `%s'\n", program_name, arg );
  std::fputs( "Valid transformation names are:", stderr );
  for( int i = 0; T_table[i].name != 0; ++i )
    std::fprintf( stderr, "  %s", T_table[i].name );
  std::fputs( "\nRotations are made counter-clockwise.\n", stderr );
  }


bool Control::set_format( const char * name ) throw()
  {
  if( std::strcmp( name, "byte" ) == 0 ) { format = byte; return true; }
  if( std::strcmp( name, "utf8" ) == 0 ) { format = utf8; return true; }
  return false;
  }
