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

class Textblock : public Rectangle
  {
  mutable std::vector< Textline * > tlpv;

public:
  Textblock( const Rectangle & r,
             std::vector< std::vector< Block * > > & blockp_matrix ) throw();
  ~Textblock() throw();
  void recognize( const Charset & charset, const Filter & filter ) throw();

//  Textline & textline( int i ) const throw();
  int textlines() const throw() { return tlpv.size(); }
  int characters() const throw();

  void print( const Control & control ) const throw();
  void dprint( const Control & control, bool graph, bool recursive ) const throw();
  void xprint( const Control & control ) const throw();
  void cmark( Page_image & page_image ) const throw();
  void lmark( Page_image & page_image ) const throw();
  };
