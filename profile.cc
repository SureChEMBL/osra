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
#include "bitmap.h"
#include "block.h"
#include "profile.h"


Profile::Profile( const Bitmap & b, Type t ) throw()
  : _bitmap( &b ), type( t ),
    _limit( -1 ), _max( -1 ), _min( -1 ), _mean( -1 ),
    _isconcave( -1 ), _isconvex( -1 ), _isflat( -1 ), _isflats( -1 ),
    _ispit( -1 ), _istpit( -1 ), _isupit( -1 ), _isvpit( -1 ), _istip( -1 ) {}


void Profile::initialize() throw()
  {
  const Bitmap & b = *_bitmap;

  switch( type )
    {
    case left :
      data.resize( b.height() ); _limit = b.width();
      for( int row = b.top(); row <= b.bottom(); ++row )
        {
        int j = b.left();
        while( j <= b.right() && !b.get_bit( row, j ) ) ++j;
        data[row-b.top()] = j - b.left();
        } break;
    case top :
      data.resize( b.width() ); _limit = b.height();
      for( int col = b.left(); col <= b.right(); ++col )
        {
        int j = b.top();
        while( j <= b.bottom() && !b.get_bit( j, col ) ) ++j;
        data[col-b.left()] = j - b.top();
        } break;
    case right :
      data.resize( b.height() ); _limit = b.width();
      for( int row = b.top(); row <= b.bottom(); ++row )
        {
        int j = b.right();
        while( j >= b.left() && !b.get_bit( row, j ) ) --j;
        data[row-b.top()] = b.right() - j;
        } break;
    case bottom :
      data.resize( b.width() ); _limit = b.height();
      for( int col = b.left(); col <= b.right(); ++col )
        {
        int j = b.bottom();
        while( j >= b.top() && !b.get_bit( j, col ) ) --j;
        data[col-b.left()] = b.bottom() - j;
        } break;
    case height :
      data.resize( b.width() ); _limit = b.height();
      for( int col = b.left(); col <= b.right(); ++col )
        {
        int u = b.top(), d = b.bottom();
        while( u <= d && !b.get_bit( u, col ) ) ++u;
        while( u <= d && !b.get_bit( d, col ) ) --d;
        data[col-b.left()] = d - u + 1;
        } break;
    case width :
      data.resize( b.height() ); _limit = b.width();
      for( int row = b.top(); row <= b.bottom(); ++row )
        {
        int l = b.left(), r = b.right();
        while( l <= r && !b.get_bit( row, l ) ) ++l;
        while( l <= r && !b.get_bit( row, r ) ) --r;
        data[row-b.top()] = r - l + 1;
        } break;
    }
  }


int Profile::mean() throw()
  {
  if( _mean < 0 )
    {
    if( _limit < 0 ) initialize();
    _mean = 0;
    for( int i = 0; i < samples(); ++i ) _mean += data[i];
    if( samples() > 1 ) _mean /= samples();
    }
  return _mean;
  }


int Profile::max() throw()
  {
  if( _max < 0 )
    {
    if( _limit < 0 ) initialize();
    _max = data[0];
    for( int i = 1; i < samples(); ++i ) if( data[i] > _max ) _max = data[i];
    }
  return _max;
  }


int Profile::max( int l, int r ) throw()
  {
  if( _limit < 0 ) initialize();
  if( r < 0 ) r = samples() - 1;
  int m = 0;
  for( int i = l; i <= r; ++i ) if( data[i] > m ) m = data[i];
  return m;
  }


int Profile::min() throw()
  {
  if( _min < 0 )
    {
    if( _limit < 0 ) initialize();
    _min = data[0];
    for( int i = 1; i < samples(); ++i ) if( data[i] < _min ) _min = data[i];
    }
  return _min;
  }


int Profile::min( int l, int r ) throw()
  {
  if( _limit < 0 ) initialize();
  if( r < 0 ) r = samples() - 1;
  int m = _limit;
  for( int i = l; i <= r; ++i ) if( data[i] < m ) m = data[i];
  return m;
  }


int Profile::operator[]( int i ) throw()
  {
  if( _limit < 0 ) initialize();
  if( i < 0 ) i = 0;
  else if( i >= samples() ) i = samples() - 1;
  return data[i];
  }


int Profile::area( int l, int r ) throw()
  {
  if( _limit < 0 ) initialize();
  if( r < 0 ) r = samples() - 1;
  int area = 0;
  for( int i = l; i <= r; ++i ) area += data[i];
  return area;
  }


bool Profile::increasing( int i ) throw()
  {
  if( _limit < 0 ) initialize();
  if( i <= 0 || i > samples() - 2 || data[samples()-1] - data[i] < 2 )
    return false;
  while( ++i < samples() ) if( data[i] < data[i-1] ) return false;
  return true;
  }


bool Profile::decreasing( int i ) throw()
  {
  if( _limit < 0 ) initialize();
  const int noise = ( std::min( samples(), _limit ) / 20 ) + 1;
  if( i < 0 || samples() - i <= 2 * noise ||
      data[samples()-noise] - data[i] > -( noise + 1 ) )
    return false;
  while( ++i < samples() - noise ) if( data[i] > data[i-1] ) return false;
  return true;
  }


bool Profile::isconcave() throw()
  {
  if( _isconcave < 0 )
    {
    _isconcave = false; if( _limit < 0 ) initialize();
    if( samples() < 5 ) return _isconcave;
    int dmax = -1, l = 0, r = 0;

    for( int i = pos( 10 ); i <= pos( 90 ); ++i )
      {
      if( data[i] > dmax ) { dmax = data[i]; l = r = i; }
      else if( data[i] == dmax ) { r = i; }
      }
    if( l > r || l < pos( 25 ) || r > pos( 75 ) ) return _isconcave;
    if( data[pos(10)] >= dmax || data[pos(90)] >= dmax ) return _isconcave;
    int imax = ( l + r ) / 2;

    for( int i = pos( 10 ); i < imax; ++i )
      if( data[i] > data[i+1] ) return _isconcave;
    for( int i = pos( 90 ); i > imax; --i )
      if( data[i] > data[i-1] ) return _isconcave;
    _isconcave = true;
    }
  return _isconcave;
  }


bool Profile::isconvex() throw()
  {
  if( _isconvex < 0 )
    {
    _isconvex = false; if( _limit < 0 ) initialize();
    if( samples() < 9 || _limit < 5 ) return _isconvex;
    int min = _limit, min_begin = 0, min_end = 0;
    int lmin = _limit, rmax = -_limit, l = 0, r = 0;

    for( int i = 1; i < samples(); ++i )
      {
      int d = data[i] - data[i-1];
      if( d < lmin ) { lmin = d; l = i - 1; }
      if( d >= rmax ) { rmax = d; r = i; }
      if( data[i] <= min )
        { min_end = i; if( data[i] < min ) { min = data[i]; min_begin = i; } }
      }
    if( l >= r || l >= pos( 25 ) || r <= pos( 75 ) ) return _isconvex;
    if( lmin >= 0 || rmax <= 0 || data[l] < 2 || data[r] < 2 ||
        3 * ( data[l] + data[r] ) <= std::min( _limit, samples() ) )
      return _isconvex;
    if( 3 * ( min_end - min_begin + 1 ) > 2 * samples() ) return _isconvex;
    if( 2 * l >= min_begin || 2 * r <= min_end + samples() - 1 ) return _isconvex;
    if( min_begin < pos( 10 ) || min_end > pos( 90 ) ) return _isconvex;

    const int noise = ( std::min( samples(), _limit ) / 30 ) + 1;
    int dmax = -_limit;
    for( int i = l + 1; i <= r; ++i )
      {
      if( i >= min_begin && i <= min_end )
        { if( data[i] <= noise ) continue; else return _isconvex; }
      int d = data[i] - data[i-1];
      if( d == 0 ) continue;
      if( d > dmax ) { if( std::abs( d ) <= noise ) ++dmax; else dmax = d; }
      else if( d < dmax - noise ) return _isconvex;
      }
    if( 2 * ( min_end - min_begin + 1 ) < samples() )
      {
      int varea = ( min_begin - l + 1 ) * data[l] / 2;
      varea += ( r - min_end + 1 ) * data[r] / 2;
      if( this->area( l, min_begin - 1 ) + this->area( min_end + 1, r ) >= varea )
        return _isconvex;
      }
    _isconvex = true;
    }
  return _isconvex;
  }


bool Profile::isflat() throw()
  {
  if( _isflat < 0 )
    {
    _isflat = false; if( _limit < 0 ) initialize();
    if( samples() < 15 ) return _isflat;
    int mn = data[samples()/2], mx = mn;

    for( int i = 1; i < samples() - 1; ++i )
      { int d = data[i]; if( d < mn ) mn = d; else if( d > mx ) mx = d; }
    _isflat = (bool)( mx - mn <= 1 + ( samples() / 30 ) );
    }
  return _isflat;
  }


bool Profile::isflats() throw()
  {
  if( _isflats < 0 )
    {
    _isflats = false; if( _limit < 0 ) initialize();
    if( samples() < 15 ) return _isflats;
    const int s1 = std::max( pos( 15 ), 3 );
    const int s2 = std::min( pos( 85 ), samples() - 4 );
    int mn = -1, mx = 0;
    for( int i = s1 + 2; i < s2; ++i )
      if( data[i-1] == data[i] ) { mn = mx = data[i]; break; }
    if( mn < 0 ) return _isflats;

    for( int i = 1; i <= s1; ++i ) if( data[i] > mx ) mx = data[i];
    for( int i = s1 + 1; i < s2; ++i )
      { int d = data[i]; if( d < mn ) mn = d; else if( d > mx ) mx = d; }
    for( int i = s2; i < samples() - 1; ++i ) if( data[i] > mx ) mx = data[i];
    _isflats = (bool)( mx - mn <= 1 + ( samples() / 30 ) );
    }
  return _isflats;
  }


bool Profile::ispit() throw()
  {
  if( _ispit < 0 )
    {
    _ispit = false; if( _limit < 0 ) initialize();
    if( samples() < 5 ) return _ispit;
    const int noise = ( std::min( samples(), _limit ) / 25 ) + 1;
    for( int i = 0; i < noise; ++i )
      if( data[i] <= noise - i || data[samples()-i-1] <= noise - i )
        return _ispit;

    const int dmin = min(), dmax = _limit / 2;
    int begin = 0, end = 0, i, ref;
    for( i = 0, ref = dmax; i < samples(); ++i )
      {
      int d = data[i];
      if( d == dmin ) { begin = i; break; }
      if( d < ref ) ref = d; else if( d > ref + noise && ref < dmax ) return _ispit;
      }
    if( begin < 2 || begin > samples() - 3 ) return _ispit;

    for( i = samples() - 1, ref = dmax; i >= begin; --i )
      {
      int d = data[i];
      if( d == dmin ) { end = i; break; }
      if( d < ref ) ref = d; else if( d > ref + noise && ref < dmax ) return _ispit;
      }
    if( end < begin || end > samples() - 3 ) return _ispit;

    for( i = begin + 1; i < end; ++i )
      if( data[i] > dmin + noise ) return _ispit;

    _ispit = true;
    }
  return _ispit;
  }


bool Profile::iscpit( const int cpos ) throw()
  {
  if( _limit < 0 ) initialize();
  if( samples() < 5 || cpos < 25 || cpos > 75 ) return false;
  int th = ( mean() < 2 ) ? 2 : mean();
  int imin = -1, mid = ( ( samples() - 1 ) * cpos ) / 100;

  for( int i = 0; i < samples() / 4; ++i )
    {
    if( data[mid+i] < th ) { imin = mid + i; break; }
    if( data[mid-i-1] < th ) { imin = mid - i - 1; break; }
    }
  if( imin < 0 ) return false;

  for( int i = imin + 1; i < samples(); ++i )
    if( data[i] > th )
      {
      for( int j = imin - 1; j >= 0; --j ) if( data[j] > th ) return true;
      break;
      }
  return false;
  }


bool Profile::islpit() throw()
  {
  if( _limit < 0 ) initialize();
  if( samples() < 5 ) return false;
  const int noise = samples() / 30;
  if( data[0] < noise + 2 ) return false;

  const int dmin = min();
  int begin = 0, ref = _limit;
  for( int i = 0; i < samples(); ++i )
    {
    int d = data[i];
    if( d == dmin ) { begin = i; break; }
    if( d < ref ) ref = d; else if( d > ref + 1 ) return false;
    }
  if( begin < 2 || 2 * begin >= samples() ) return false;
  return true;
  }


bool Profile::istpit() throw()
  {
  if( _istpit < 0 )
    {
    if( _limit < 0 ) initialize();
    if( _limit < 5 || samples() < 5 || !ispit() )
      { _istpit = false; return _istpit; }

    const int noise = ( std::min( _limit, samples() ) / 30 ) + 1;
    int l = -1, r = 0;
    for( int i = 0; i < samples(); ++i )
      if( data[i] <= noise ) { r = i; if( l < 0 ) l = i; }
    _istpit = (bool)( l > 0 && 4 * ( r - l + 1 ) < samples() );
    }
  return _istpit;
  }


bool Profile::isupit() throw()
  {
  if( _isupit < 0 )
    {
    _isupit = false; if( _limit < 0 ) initialize();
    if( samples() < 5 ) return _isupit;

    int th = ( mean() < 2 && range() > 2 ) ? 2 : mean();
    int status = 0, ucount = 0, lcount = 0, umean =0, lmean = 0;
    for( int i = 0; i < samples(); ++i )
      {
      int d = data[i];
      switch( status )
        {
        case 0: if( d < th )
                  { if( i < pos( 25 ) || i > pos( 70 ) ) return _isupit;
                  status = 1; break; }
                if( d > th ) { ++ucount; umean += d; }
          break;
        case 1: if( d > th )
                  { if( i < pos( 30 ) || i > pos( 75 ) ) return _isupit;
                  status = 2; break; }
                if( d < th ) { ++lcount; lmean += d; }
          break;
        case 2: if( d < th ) return _isupit;
                if( d > th ) { ++ucount; umean += d; }
          break;
        }
      }
    if( ucount > 1 ) umean /= ucount;
    if( lcount > 1 ) lmean /= lcount;
    _isupit = (bool)( status == 2 && umean - lmean > range() / 2 );
    }
  return _isupit;
  }


bool Profile::isvpit() throw()
  {
  if( _isvpit < 0 )
    {
    if( _limit < 0 ) initialize();
    if( _limit < 5 || samples() < 5 || !ispit() )
      { _isvpit = false; return _isvpit; }

    const int noise = _limit / 20;
    const int level = ( _limit / 10 ) + 2;
    int ll = -1, ln = -1, rl = -1, rn = -1;
    for( int i = 0; i < samples(); ++i )
      if( data[i] <= level )
        {
        rl = i; if( ll < 0 ) ll = i;
        if( data[i] <= noise ) { rn = i; if( ln < 0 ) ln = i; }
        }
    const int wl = rl - ll + 1, wn = rn - ln + 1;
    _isvpit = (bool)( ln > 0 && 2 * wl <= samples() + 1 && wl - wn <= 2 * ( level - noise ) );
    }
  return _isvpit;
  }


bool Profile::istip() throw()
  {
  if( _istip < 0 )
    {
    _istip = false; if( _limit < 0 ) initialize();
    if( samples() < 5 ) return _istip;

    int th = ( mean() < 2 && range() > 2 ) ? 2 : mean(); if( th < 2 ) ++th;
    int lth = data[0], rth = data[samples()-1];
    int begin = 0, end = samples() - 1;
    for( int i = 1, j = std::max( 2, samples() / 10 ); i < j; ++i )
      {
      if( data[i] < lth ) { lth = data[i]; begin = i; }
      if( data[samples()-1-i] < rth )
        { rth = data[samples()-1-i]; end = samples() - 1 - i; }
      }
    if( lth >= th || rth >= th ) return _istip;
    if( 3 * lth >= 2 * range() || 3 * rth >= 2 * range() ) return _istip;
    th = std::max( lth, rth );
    int status = 0;
    for( int i = begin + 1; i < end; ++i )
      switch( status )
        {
        case 0: if( data[i] > th + 1 ) status = 1; break;
        case 1: if( data[i] > th + 1 ) status = 2; else status = 0; break;
        case 2: if( data[i] <= th ) status = 3; break;
        case 3: if( data[i] > th + 1 ) return _istip;
        }
    _istip = (bool)( status >= 2 );
    }
  return _istip;
  }


bool Profile::isctip( const int cpos ) throw()
  {
  if( _limit < 0 ) initialize();
  if( samples() < 5 || cpos < 25 || cpos > 75 ) return false;
  const int mid = ( ( samples() - 1 ) * cpos ) / 100;
  int th = std::max( 2, std::min( mean(), _limit / 3 ) );
  int imax = -1;

  for( int i = 0; i < samples() / 4; ++i )
    {
    if( data[mid+i] > th ) { imax = mid + i; break; }
    if( data[mid-i-1] > th ) { imax = mid - i - 1; break; }
    }
  if( imax < 0 && mean() == 0 )
    {
    --th;
    for( int i = 0; i < samples() / 4; ++i )
      {
      if( data[mid+i] > th ) { imax = mid + i; break; }
      if( data[mid-i-1] > th ) { imax = mid - i - 1; break; }
      }
    }
  if( imax < 0 ) return false;

  for( int i = imax + 1; i < samples(); ++i )
    if( data[i] < th )
      {
      for( int j = imax - 1; j >= 0; --j ) if( data[j] < th ) return true;
      break;
      }
  return false;
  }


int Profile::imaximum() throw()
  {
  if( _limit < 0 ) initialize();
  const int margin = ( samples() / 30 ) + 1;
  int mbegin = 0, mend, mvalue = 0;

  for( int i = margin; i < samples() - margin; ++i )
    if( data[i] > mvalue ) { mvalue = data[i]; mbegin = i; }

  for( mend = mbegin + 1; mend < samples(); ++mend )
    if( data[mend] < mvalue ) break;

  return ( mbegin + mend - 1 ) / 2;
  }


int Profile::iminimum( int m, int th ) throw()
  {
  if( _limit < 0 ) initialize();
  const int margin = ( samples() / 30 ) + 1;
  if( samples() < 2 * margin ) return 0;
  if( th < 2 ) th = ( mean() < 2 ) ? 2 : mean();
  int minima = 0, status = 0;
  int begin = 0, end, value = _limit + 1;

  for( end = margin; end < samples() - margin; ++end )
    {
    if( status == 0 )
      { if( data[end] < th ) { status = 1; ++minima; begin = end; } }
    else if( data[end] > th )
      { if( minima == m + 1 ) { --end; break; } else status = 0; }
    }
  if( end >= samples() ) --end;
  if( minima != m + 1 ) return 0;

  for( int i = begin; i <= end; ++i )
    if( data[i] < value ) { value = data[i]; begin = i; }
  for( ; end >= begin; --end ) if( data[end] == value ) break;
  return ( begin + end ) / 2;
  }


int Profile::minima( int th ) throw()
  {
  if( _limit < 0 ) initialize();
  if( !samples() ) return 0;
  if( th < 1 ) th = ( mean() < 2 ) ? 2 : mean();
  const int noise = _limit / 40;
  const int dth = th - ( ( noise + 1 ) / 2 ), uth = th + ( noise / 2 );
  if( dth < 1 ) return 1;
  int minima = ( data[0] < dth ) ? 1 : 0;
  int status = ( minima ) ? 1 : 0;

  for( int i = 1; i < samples(); ++i )
    switch( status )
      {
      case 0: if( data[i] < dth ) { status = 1; ++minima; } break;
      case 1: if( data[i] > uth ) status = 0; break;
      }
  return minima;
  }


bool Profile::straight( int * _dy ) throw()
  {
  if( _limit < 0 ) initialize();
  if( samples() < 5 ) return false;

  const int xl = ( samples() / 30 ) + 1, yl = ( data[xl] + data[xl+1] ) / 2 ;
  const int xr = samples() - xl - 1,     yr = ( data[xr-1] + data[xr] ) / 2 ;
  const int dx = xr - xl, dy = yr - yl;
  if( dx <= 0 ) return false;
  const int dmax = dx * ( ( samples() / 20 ) + 2 );
  int faults = samples() / 10;
  for( int i = 0; i < samples(); ++i )
    {
    int y = ( dx * yl ) + ( ( i - xl ) * dy );
    int d = std::abs( ( dx * data[i] ) - y );
    if( d >= dmax && ( ( dx * data[i] ) < y || ( i >= xl && i <= xr ) ) )
      if( d > dmax || ( d == dmax && --faults < 0 ) ) return false;
    }
  if( _dy ) *_dy = dy;
  return true;
  }
