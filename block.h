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

class Block : public Bitmap
  {
  std::vector< Bitmap * > bpv;		// vector of holes

public:
  Block( const int l, const int t, const int r, const int b ) throw()
    : Bitmap( l, t, r, b ) {}

  Block( const Bitmap & source, const Rectangle & re ) throw()
    : Bitmap( source, re ) {}

  Block( const Block & b ) throw();
  Block & operator=( const Block & b ) throw();

  ~Block() throw();

  using Bitmap::left;
  using Bitmap::top;
  using Bitmap::right;
  using Bitmap::bottom;
  using Bitmap::height;
  using Bitmap::width;
  void left  ( const int l ) throw();
  void top   ( const int t ) throw();
  void right ( const int r ) throw();
  void bottom( const int b ) throw();
  void height( const int h ) throw() { bottom( top() + h - 1 ); }
  void width ( const int w ) throw() { right( left() + w - 1 ); }

  const Bitmap & hole( const int i ) const throw();
  int holes() const throw() { return bpv.size(); }
  //  id = 1 for block dots, negative for hole dots, 0 otherwise
  int id( const int row, const int col ) const throw();

  void print( FILE * outfile ) const throw();

  void fill_hole( const int i ) throw();
  void find_holes() throw();
  };
