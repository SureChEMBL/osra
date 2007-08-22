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

class Profile
  {
public:
  enum Type { left, top, right, bottom, height, width };

private:
  const Bitmap * _bitmap;	// Bitmap to witch this profile belongs
  Type type;
  int _limit, _max, _min, _mean;
  signed char _isconcave, _isconvex, _isflat, _isflats,
              _ispit, _istpit, _isupit, _isvpit, _istip;
  std::vector< int > data;
  void initialize() throw();
  int mean() throw();

public:
  Profile( const Bitmap & b, Type t ) throw();

//  const Bitmap * bitmap() const throw() { return _bitmap; }

  int limit() throw() { if( _limit < 0 ) initialize(); return _limit; }
  int max() throw();
  int max( int l, int r = -1 ) throw();
  int min() throw();
  int min( int l, int r = -1 ) throw();
  int operator[]( int i ) throw();
  int pos( int p ) throw() { return ( ( samples() - 1 ) * p ) / 100; }
  int range() throw() { return max() - min(); }
  int samples() throw() { if( _limit < 0 ) initialize(); return data.size(); }

  int  area( int l = 0, int r = -1 ) throw();
  bool increasing( int i = 1 ) throw();
  bool decreasing( int i = 1 ) throw();
  bool isconcave() throw();
  bool isconvex() throw();
  bool isflat() throw();
  bool isflats() throw();
  bool ispit() throw();
  bool iscpit( const int cpos = 50 ) throw();
  bool islpit() throw();
  bool istpit() throw();
  bool isupit() throw();
  bool isvpit() throw();
  bool istip() throw();
  bool isctip( const int cpos = 50 ) throw();
  int  imaximum() throw();
  int  iminimum( int m = 0, int th = -1 ) throw();
  int  minima( int th = -1 ) throw();
  bool straight( int * dy ) throw();
  };
