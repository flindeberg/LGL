/*
 *  Copyright (c) 2026 Fredrik Lindeberg
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 */
#ifndef CELLLIST_HPP_
#define CELLLIST_HPP_

//---------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include "configs.h"

#if defined(__aarch64__)
#include <arm_neon.h>
#define CELLLIST_NEON 1
#endif

//---------------------------------------------------------------
// A flat copy, made every iteration, of the particles that sit in the
// voxels of the active box: one contiguous block per voxel holding the
// coordinates (one array per dimension), the radius, an "anchor factor"
// (2 for an anchor, which pushes twice as hard, 1 otherwise) and the
// particle itself. Each block is padded to a multiple of LANES with
// far-away sentinels, so the pair kernel can take whole lanes without a
// tail loop. The repulsion phase reads these arrays sequentially instead
// of chasing pointers to 160-byte particles, and can use SIMD.
//
// The work is split into units of at most UNIT particles of one voxel,
// so a hub whose thousands of leaves share a voxel is spread over the
// threads instead of being one thread's job while the others wait.
//---------------------------------------------------------------

struct CellList {
  enum { LANES = 4 , UNIT = 256 };
  static constexpr float SENTINEL = 1e18f; // never in range, no overflow

  struct Unit { long box; unsigned begin , end; };

  std::vector< float > x[DIMENSION];
  std::vector< float > r;
  std::vector< float > af;
  std::vector< Node* > p;
  // Per voxel of the box (row-major): start of its block and occupants
  std::vector< unsigned > offset;
  std::vector< unsigned > count;
  std::vector< Unit > units;
  FixedVec_l lo;            // corner of the box, in voxel coordinates
  long extent[DIMENSION];
  long boxSize;

  static unsigned padded( unsigned n ) { return ( n + LANES - 1 ) / LANES * LANES; }

  // Box index to voxel coordinates
  void coords( long b , long c[DIMENSION] ) const {
    for ( unsigned int d=0; d<DIMENSION; ++d ) {
      c[d] = lo[d] + b % extent[d];
      b /= extent[d];
    }
  }

  // Voxel coordinates to box index, or -1 when outside the box
  long boxIndex( const long c[DIMENSION] ) const {
    long b = 0 , stride = 1;
    for ( unsigned int d=0; d<DIMENSION; ++d ) {
      const long cc = c[d] - lo[d];
      if ( cc < 0 || cc >= extent[d] ) { return -1; }
      b += cc * stride;
      stride *= extent[d];
    }
    return b;
  }

  static long gridIndex( const Grid_t& g , const long c[DIMENSION] ) {
    long i = 0;
    for ( unsigned int d=0; d<DIMENSION; ++d ) { i += c[d] * g.voxelsPerDim( d ); }
    return i;
  }

  // Serial part: the box, the block offsets and the work units. Reads
  // the occupancy of every voxel in the box once.
  void prepare( const Grid_t& grid , const FixedVec_l& boxLo , const FixedVec_l& boxHi ) {
    lo = boxLo;
    boxSize = 1;
    for ( unsigned int d=0; d<DIMENSION; ++d ) {
      extent[d] = boxHi[d] - boxLo[d] + 1;
      if ( extent[d] < 0 ) { extent[d] = 0; }
      boxSize *= extent[d];
    }
    offset.resize( boxSize + 1 );
    count.resize( boxSize );
    units.clear();
    unsigned total = 0;
    long c[DIMENSION];
    for ( long b=0; b<boxSize; ++b ) {
      coords( b , c );
      const unsigned n = grid[ gridIndex( grid , c ) ].occupancy();
      count[b] = n;
      offset[b] = total;
      total += padded( n );
      for ( unsigned i=0; i<n; i+=UNIT ) {
	Unit u = { b , i , std::min( i + (unsigned) UNIT , n ) };
	units.push_back( u );
      }
    }
    offset[boxSize] = total;
    for ( unsigned int d=0; d<DIMENSION; ++d ) { x[d].resize( total ); }
    r.resize( total );
    af.resize( total );
    p.resize( total );
  }

  // Parallel part: copy the particles of the box voxels [b0,b1) into
  // their blocks. Disjoint ranges may be filled concurrently.
  void fill( const Grid_t& grid , long b0 , long b1 ) {
    long c[DIMENSION];
    for ( long b=b0; b<b1; ++b ) {
      const unsigned n = count[b];
      if ( n == 0 ) { continue; }
      coords( b , c );
      const Voxel_t& v = grid[ gridIndex( grid , c ) ];
      const unsigned o = offset[b];
      unsigned k = 0;
      for ( Voxel_t::const_occupant_iterator it = v.begin(); it != v.end() && k < n; ++it, ++k ) {
	Node& q = **it;
	for ( unsigned int d=0; d<DIMENSION; ++d ) { x[d][o+k] = q.X()[d]; }
	r[o+k] = q.radius();
	af[o+k] = q.isAnchor() ? 2.f : 1.f;
	p[o+k] = &q;
      }
      for ( ; k < padded( n ); ++k ) {
	for ( unsigned int d=0; d<DIMENSION; ++d ) { x[d][o+k] = SENTINEL; }
	r[o+k] = 0;
	af[o+k] = 0;
	p[o+k] = 0;
      }
    }
  }
};

//---------------------------------------------------------------
// The repulsive spring between particles closer than the cutoff: with
// d the distance and eq the equilibrium (both the neighbourhood radius),
// the force on particle 1 is -k (d - eq) / d * (x1 - x2), times 2 when
// particle 2 is an anchor. Ellipse factors scale the differences before
// the force is computed (not in the cutoff or collision tests), as
// ParticleInteractionHandler::springRepulsiveInteraction does.
//---------------------------------------------------------------

struct PairParams {
  float eq;              // equilibrium distance
  float eq2;             // cutoff distance squared
  float k;               // spring constant
  float ef[DIMENSION];   // ellipse factors
};

// Adds to f the force on a particle at xi with radius ri from the block
// [o, o+n) of the cell list (n padded), skipping slot `skip` (its own
// slot, or -1). Returns how many particles in the block collide with it
// (closer than the sum of the radii); those give no force, the caller
// adds noise instead. The generic version works for any dimension; on
// arm64 the 2D layout uses a NEON version, four pairs at a time.
inline unsigned repulsionFromBlockScalar( const CellList& cl , unsigned o , unsigned n , long skip ,
					  const float xi[DIMENSION] , float ri ,
					  const PairParams& pp , float f[DIMENSION] )
{
  unsigned collisions = 0;
  for ( unsigned j=0; j<n; ++j ) {
    if ( (long) j == skip ) { continue; }
    float dx[DIMENSION];
    float d2 = 0;
    for ( unsigned int d=0; d<DIMENSION; ++d ) {
      dx[d] = xi[d] - cl.x[d][o+j];
      d2 += dx[d] * dx[d];
    }
    if ( ! ( d2 < pp.eq2 ) ) { continue; }
    const float rr = ri + cl.r[o+j];
    if ( d2 <= rr * rr ) { ++collisions; continue; }
    float ds[DIMENSION];
    float m2 = 0;
    for ( unsigned int d=0; d<DIMENSION; ++d ) {
      ds[d] = dx[d] * pp.ef[d];
      m2 += ds[d] * ds[d];
    }
    const float mag = std::sqrt( m2 );
    const float sc = -pp.k * ( mag - pp.eq ) / mag * cl.af[o+j];
    for ( unsigned int d=0; d<DIMENSION; ++d ) { f[d] += ds[d] * sc; }
  }
  return collisions;
}

#ifdef CELLLIST_NEON

// Only valid when DIMENSION == 2 (the other coordinates are ignored)
inline unsigned repulsionFromBlockNeon2D( const CellList& cl , unsigned o , unsigned n , long skip ,
					  const float xi[DIMENSION] , float ri ,
					  const PairParams& pp , float f[DIMENSION] )
{
  const float32x4_t vxi = vdupq_n_f32( xi[0] ), vyi = vdupq_n_f32( xi[1] );
  const float32x4_t vri = vdupq_n_f32( ri );
  const float32x4_t veq2 = vdupq_n_f32( pp.eq2 ), veq = vdupq_n_f32( pp.eq );
  const float32x4_t vmk = vdupq_n_f32( -pp.k );
  const float32x4_t vex = vdupq_n_f32( pp.ef[0] ), vey = vdupq_n_f32( pp.ef[1] );
  const float32x4_t vzero = vdupq_n_f32( 0 );
  const uint32_t idx0[4] = { 0 , 1 , 2 , 3 };
  uint32x4_t vidx = vld1q_u32( idx0 );
  const uint32x4_t vskip = vdupq_n_u32( skip < 0 ? 0xffffffffu : (uint32_t) skip );
  const uint32x4_t vfour = vdupq_n_u32( 4 );
  float32x4_t vfx = vzero , vfy = vzero;
  uint32x4_t vcol = vdupq_n_u32( 0 ); // collisions per lane (a mask lane is -1)
  const float * xs = cl.x[0].data() + o;
  const float * ys = cl.x[1].data() + o;
  const float * rs = cl.r.data() + o;
  const float * as = cl.af.data() + o;
  for ( unsigned j=0; j<n; j+=4 ) {
    const float32x4_t dx = vsubq_f32( vxi , vld1q_f32( xs + j ) );
    const float32x4_t dy = vsubq_f32( vyi , vld1q_f32( ys + j ) );
    const float32x4_t d2 = vmlaq_f32( vmulq_f32( dx , dx ) , dy , dy );
    const float32x4_t rr = vaddq_f32( vri , vld1q_f32( rs + j ) );
    uint32x4_t in = vcltq_f32( d2 , veq2 );
    in = vbicq_u32( in , vceqq_u32( vidx , vskip ) );
    const uint32x4_t col = vandq_u32( vcleq_f32( d2 , vmulq_f32( rr , rr ) ) , in );
    vcol = vsubq_u32( vcol , col );
    const float32x4_t dsx = vmulq_f32( dx , vex ), dsy = vmulq_f32( dy , vey );
    const float32x4_t mag = vsqrtq_f32( vmlaq_f32( vmulq_f32( dsx , dsx ) , dsy , dsy ) );
    float32x4_t sc = vdivq_f32( vmulq_f32( vmk , vsubq_f32( mag , veq ) ) , mag );
    sc = vmulq_f32( sc , vld1q_f32( as + j ) );
    // Lanes out of range, colliding or the particle itself get 0 (a bit
    // select, so a NaN from 0/0 in those lanes is dropped)
    sc = vbslq_f32( vbicq_u32( in , col ) , sc , vzero );
    vfx = vmlaq_f32( vfx , dsx , sc );
    vfy = vmlaq_f32( vfy , dsy , sc );
    vidx = vaddq_u32( vidx , vfour );
  }
  f[0] += vaddvq_f32( vfx );
  f[1] += vaddvq_f32( vfy );
  return vaddvq_u32( vcol );
}

#endif

inline unsigned repulsionFromBlock( const CellList& cl , unsigned o , unsigned n , long skip ,
				    const float xi[DIMENSION] , float ri ,
				    const PairParams& pp , float f[DIMENSION] )
{
#ifdef CELLLIST_NEON
  if ( DIMENSION == 2 ) {
    return repulsionFromBlockNeon2D( cl , o , n , skip , xi , ri , pp , f );
  }
#endif
  return repulsionFromBlockScalar( cl , o , n , skip , xi , ri , pp , f );
}


//---------------------------------------------------------------

#endif
