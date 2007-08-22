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

class Rectangle
  {
  int _left, _top, _right, _bottom;

public:
  Rectangle( const int l, const int t, const int r, const int b ) throw();

  void left  ( const int l ) throw();
  void top   ( const int t ) throw();
  void right ( const int r ) throw();
  void bottom( const int b ) throw();
  void height( const int h ) throw();
  void width ( const int w ) throw();
  void add_point( const int row, const int col ) throw();
  void add_rectangle( const Rectangle & re ) throw();
  void enlarge( const int scale ) throw();
  void move( const int row, const int col ) throw();

  int left()    const throw() { return _left;   }
  int top()     const throw() { return _top;    }
  int right()   const throw() { return _right;  }
  int bottom()  const throw() { return _bottom; }
  int height()  const throw() { return _bottom - _top + 1; }
  int width()   const throw() { return _right - _left + 1; }
  int size()    const throw() { return height() * width(); }
  int hcenter() const throw() { return ( _left + _right ) / 2; }
  int vcenter() const throw() { return ( _top + _bottom ) / 2; }
  int hpos( const int p ) const throw()
    { return _left + ( ( ( _right - _left ) * p ) / 100 ); }
  int vpos( const int p ) const throw()
    { return _top + ( ( ( _bottom - _top ) * p ) / 100 ); }

  bool operator==( const Rectangle & re ) const throw()
    { return ( _left == re._left && _top == re._top && _right == re._right && _bottom == re._bottom ); }
  bool operator!=( const Rectangle & re ) const throw() { return !( *this == re ); }

  bool includes( const Rectangle & re ) const throw();
  bool includes( const int row, const int col ) const throw();
  bool strictly_includes( const Rectangle & re ) const throw();
  bool strictly_includes( const int row, const int col ) const throw();
  bool includes_hcenter( const Rectangle & re ) const throw();
  bool includes_vcenter( const Rectangle & re ) const throw();
  bool h_includes( const Rectangle & re ) const throw();
  bool h_includes( const int col ) const throw();
  bool v_includes( const Rectangle & re ) const throw();
  bool v_includes( const int row ) const throw();
  bool h_overlaps( const Rectangle & re ) const throw();
  bool v_overlaps( const Rectangle & re ) const throw();
  bool is_hcentred_in( const Rectangle & re ) const throw();
  bool is_vcentred_in( const Rectangle & re ) const throw();
  bool h_precedes( const Rectangle & re ) const throw();
  bool v_precedes( const Rectangle & re ) const throw();

  int distance( const Rectangle & re ) const throw();
  int distance( const int row, const int col ) const throw();
  int h_distance( const Rectangle & re ) const throw();
  int h_distance( const int col ) const throw();
  int v_distance( const Rectangle & re ) const throw();
  int v_distance( const int row ) const throw();
  };
