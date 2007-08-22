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

// Rationals are kept normalized at all times.
// Invariant = ( gcd( num, den ) == 1 && den > 0 ).
//
class Rational
  {
  int num, den;

  void normalize() throw();

public:
  explicit Rational( const int n = 0 ) throw() : num( n ), den( 1 ) {}
  Rational( const int n, const int d ) throw() : num( n ), den( d ) { normalize(); }

  Rational & assign( const int n, const int d ) throw()
    { num = n; den = d; normalize(); return *this; }
  Rational & operator=( const int n ) throw()
    { num = n; den = 1; return *this; }

  int numerator() const throw() { return num; }
  int denominator() const throw() { return den; }

  Rational & operator+() throw() { return *this; }		// Unary plus
  Rational operator-() const throw()				// Unary minus
    { Rational tmp = *this; tmp.num = -tmp.num; return tmp; }

  Rational & operator+=( const Rational & r ) throw();
  Rational & operator-=( const Rational & r ) throw() { return operator+=( -r ); }
  Rational & operator*=( const Rational & r ) throw();
  Rational & operator/=( const Rational & r ) throw();

  Rational & operator+=( const int n ) throw() { num += n * den; return *this; }
  Rational & operator-=( const int n ) throw() { num -= n * den; return *this; }
  Rational & operator*=( const int n ) throw() { return operator*=( Rational( n ) ); }
  Rational & operator/=( const int n ) throw() { return operator/=( Rational( n ) ); }

  Rational operator+( const Rational & r ) const throw()
    { Rational tmp = *this; tmp += r; return tmp; }
  Rational operator-( const Rational & r ) const throw()
    { Rational tmp = *this; tmp -= r; return tmp; }
  Rational operator*( const Rational & r ) const throw()
    { Rational tmp = *this; tmp *= r; return tmp; }
  Rational operator/( const Rational & r ) const throw()
    { Rational tmp = *this; tmp /= r; return tmp; }

  Rational operator+( const int n ) const throw()
    { Rational tmp = *this; tmp += n; return tmp; }
  Rational operator-( const int n ) const throw()
    { Rational tmp = *this; tmp -= n; return tmp; }
  Rational operator*( const int n ) const throw()
    { Rational tmp = *this; tmp *= n; return tmp; }
  Rational operator/( const int n ) const throw()
    { Rational tmp = *this; tmp /= n; return tmp; }

  Rational & operator++() throw() { num += den; return *this; }	// prefix
  Rational operator++( int ) throw()				// suffix
    { Rational tmp = *this; num += den; return tmp; }
  Rational & operator--() throw() { num -= den; return *this; }	// prefix
  Rational operator--( int ) throw()				// suffix
    { Rational tmp = *this; num -= den; return tmp; }

  bool operator!() const throw() { return !num; }

  bool operator==( const Rational & r ) const throw()
    { return ( num == r.num && den == r.den ); }
  bool operator==( const int n ) const throw()
    { return ( num == n && den == 1 ); }
  bool operator!=( const Rational & r ) const throw() { return !( *this == r ); }
  bool operator!=( const int n ) const throw()        { return !( *this == n ); }

  bool operator< ( const Rational & r ) const throw();
  bool operator<=( const Rational & r ) const throw() { return ( *this < r || *this == r ); }
  bool operator> ( const Rational & r ) const throw() { return !( *this <= r ); }
  bool operator>=( const Rational & r ) const throw() { return !( *this <  r ); }

  bool operator< ( const int n ) const throw() { return operator< ( Rational( n ) ); }
  bool operator<=( const int n ) const throw() { return operator<=( Rational( n ) ); }
  bool operator> ( const int n ) const throw() { return operator> ( Rational( n ) ); }
  bool operator>=( const int n ) const throw() { return operator>=( Rational( n ) ); }

  Rational abs() const throw()
    { if( num >= 0 ) return *this; else return -*this; }
  Rational inverse() const throw();
  int round() const throw();
  int trunc() const throw() { return ( num / den ); }

  int parse( const char * ptr ) throw();
  const std::string to_decimal( const int iwidth = 1, int prec = -2 ) const throw();
  };


inline Rational operator+( const int n, const Rational & r ) throw() { return r + n; }
inline Rational operator-( const int n, const Rational & r ) throw() { return -r + n; }
inline Rational operator*( const int n, const Rational & r ) throw() { return r * n; }
inline Rational operator/( const int n, const Rational & r ) throw() { return Rational( n ) / r; }

inline bool operator< ( const int n, const Rational & r ) throw() { return r >  n; }
inline bool operator<=( const int n, const Rational & r ) throw() { return r >= n; }
inline bool operator> ( const int n, const Rational & r ) throw() { return r <  n; }
inline bool operator>=( const int n, const Rational & r ) throw() { return r <= n; }
