/* Added result( const Control & control )
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

class Character : public Rectangle
  {
  struct Guess
    {
    int code;
    int value;
    Guess( int c, int v ) throw() : code( c ), value( v ) {}
    };

  std::vector< Block * > bpv;		// the blocks forming this Character
  std::vector< Guess > gv;		// vector of possible char codes
					// and their associated values

  void recognize11( const Charset & charset, const Rectangle & charbox ) throw();
    void recognize110( const Charset & charset, const Rectangle & charbox ) throw();
    void recognize111( const Charset & charset, const Rectangle & charbox ) throw();
    void recognize112( const Rectangle & charbox ) throw();
  void recognize12( const Charset & charset, const Rectangle & charbox ) throw();
  void recognize13( const Charset & charset, const Rectangle & charbox ) throw();

public:
  Character( Block * p ) throw()
    : Rectangle( *p ), bpv( 1, p ) {}

  Character( const Rectangle & re, int code, int value ) throw()
    : Rectangle( re ), gv( 1, Guess( code, value ) ) {}

  Character( const Character & c ) throw();
  Character & operator=( const Character & c ) throw();

  ~Character() throw();

  int area() const throw();
  Block & block( int i ) throw();
  int blocks() const throw() { return bpv.size(); }
  Block & main_block() throw();

  void shift_blockp( Block * p ) throw();

  void add_guess( int code, int value ) throw();
  void clear_guesses() throw() { gv.clear(); }
  void insert_guess( int i, int code, int value ) throw();
  void delete_guess( int i ) throw();
  void only_guess( int code, int value ) throw();
  void swap_guesses( int i, int j ) throw();
  const Guess & guess( int i ) const throw();
  int guesses() const throw() { return gv.size(); }
  bool maybe( int code ) const throw();
//  bool maybe_digit() const throw();
//  bool maybe_letter() const throw();

  void join( Character & c ) throw();
  void print( const Control & control ) const throw();
  char result( const Control & control ) const throw();

  void dprint( const Control & control, const Rectangle & charbox,
               bool graph, bool recursive ) const throw();
  void xprint( const Control & control ) const throw();

  void recognize1( const Charset & charset, const Rectangle & charbox ) throw();
  void apply_filter( const Filter & filter ) throw();
  };
