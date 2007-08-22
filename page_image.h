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

class Rational;
class Track;

class Page_image : public Rectangle		// left,top is always 0,0
  {
public:
  struct Error
    {
    const char * s;
    Error( const char * p ) { s = p; }
    };

private:
  std::vector< std::vector< unsigned char > > data;	// 256 level greymap
  std::vector< Rectangle > rv;			// layout zones
  std::vector< unsigned char > tv;		// threshold for every zone
  unsigned char _maxval, _threshold;		// x > threshold == white

  void read_p1( FILE * f, const bool invert ) throw( Error );
  void read_p4( FILE * f, const bool invert ) throw( Error );
  void read_p2( FILE * f, const Rational & th, const bool invert ) throw( Error );
  void read_p5( FILE * f, const Rational & th, const bool invert ) throw( Error );
  void read_p3( FILE * f, const Rational & th, const bool invert ) throw( Error );
  void read_p6( FILE * f, const Rational & th, const bool invert ) throw( Error );
  void find_columns( const Rectangle & rin, bool recursive ) throw();
  void find_rows( const Rectangle & rin, bool recursive ) throw();

public:
  // Creates a Page_image from a pbm, pgm or ppm file
  Page_image( FILE * f, const Rational & th, const bool invert ) throw( Error );

  // Creates a reduced Page_image
  Page_image( const Page_image & source, const int scale ) throw();

  // Creates a reduced, b/w Page_image
  Page_image( const Page_image & source, const int scale, const Rational & th ) throw();

  using Rectangle::left;
  using Rectangle::top;
  using Rectangle::right;
  using Rectangle::bottom;
  using Rectangle::height;
  using Rectangle::width;
  void left  ( int ) throw( Error ) { throw Error( "Page_image resize not allowed." ); }
  void top   ( int ) throw( Error ) { left( 0 ); }
  void right ( int ) throw( Error ) { left( 0 ); }
  void bottom( int ) throw( Error ) { left( 0 ); }
  void height( int ) throw( Error ) { left( 0 ); }
  void width ( int ) throw( Error ) { left( 0 ); }

  bool get_bit( const int row, const int col ) const throw()
    { return data[row-top()][col-left()] <= _threshold; }
  bool get_bit( const int row, const int col, const unsigned char th ) const throw()
    { return data[row-top()][col-left()] <= th; }
  void set_bit( const int row, const int col, const bool bit ) throw()
    { data[row-top()][col-left()] = ( bit ? 0 : _maxval ); }

  unsigned char maxval() const throw() { return _maxval; }
  unsigned char threshold() const throw() { return _threshold; }
  unsigned char threshold( const int i ) const throw() { return tv[i]; }
  const Rectangle & rectangle( const int i ) const throw() { return rv[i]; }
  int zones() const throw() { return rv.size(); }

  void adapt_thresholds() throw();
  int  analyse_layout( const int layout_level ) throw();
  bool crop( const Rational ltrb[4], const Rational & th ) throw();
  void draw_rectangle( const Rectangle & re ) throw();
  void draw_track( const Track & tr ) throw();
  void histogramize( const bool vertical = true ) throw();
  bool save( FILE * f, const char filetype, const int i = -1 ) const throw();
  bool scale( int n ) throw( Error );
  void transform( const Transformation & t ) throw();
  };
