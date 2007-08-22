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

class Page_image;
class Rational;

class Textline : public Track
  {
  Character * _big_initial;
  mutable std::vector< Character * > cpv;

public:
  Textline() throw() : _big_initial( 0 ) {}
  Textline( const Textline & tl ) throw();
  Textline & operator=( const Textline & tl ) throw();
  ~Textline() throw();
  void set_track() throw();

  const Character * big_initial() const throw() { return _big_initial; }
  void big_initial( Character * p ) throw();
//  void delete_big_initial() throw();

  void delete_character( int i ) throw();
  int  shift_characterp( Character * p ) throw();
  bool insert_space( int i, bool tab = false ) throw();
  void join( Textline & tl ) throw();

  Character & character( int i ) const throw();
  int characters() const throw() { return cpv.size(); }
  Rectangle charbox( const Character & c ) const throw();
  int width() const throw()
    { return cpv.size() ? cpv.back()->right() - cpv[0]->left() : 0; }

  int mean_height() const throw();
  Rational mean_width() const throw();
  Rational mean_gap_width( const int first = 0, int last = -1 ) const throw();
  int mean_hcenter() const throw();
  int mean_vcenter() const throw();

  void print( const Control & control ) const throw();
  void dprint( const Control & control, bool graph, bool recursive ) const throw();
  void xprint( const Control & control ) const throw();
  void cmark( Page_image & page_image ) const throw();

  void recognize1( const Charset & charset ) const throw();
  void recognize2( const Charset & charset ) throw();
  void apply_filter( const Filter & filter ) throw();
  };
