// This file has functions which the threads
// could call.

#ifndef CALCFUNCS_H_
#define CALCFUNCS_H_

#include <atomic>
#include <climits>
#include "configs.h"
#include "fixedVecOperations.hpp"
#include "sphere.hpp"
#include "cellList.hpp"
#include <iostream>
#include <iomanip>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>

//------------------------------------------------

using namespace boost;
typedef adjacency_list< vecS, vecS, directedS > out_graph;
// Edges as pairs of vertex indices. The threads walk these instead of the
// boost graph's edge set: a contiguous block each, no std::advance through
// the whole set per thread.
typedef std::vector< std::pair< unsigned , unsigned > > EdgeList;

//------------------------------------------------

// An empty box (min > max), and growing a box by the voxel with the given
// grid index (index = x + y*nx + z*nx*ny)
inline void emptyBox( FixedVec_l& lo , FixedVec_l& hi ) {
  for ( unsigned int d=0; d<DIMENSION; ++d ) { lo[d] = LONG_MAX; hi[d] = -1; }
}

inline void addVoxelToBox( const Grid_t& grid , long index ,
			   FixedVec_l& lo , FixedVec_l& hi ) {
  for ( unsigned int d=0; d<DIMENSION; ++d ) {
    const long n = grid.voxelsPerEdge( d );
    const long c = index % n;
    index /= n;
    if ( c < lo[d] ) { lo[d] = c; }
    if ( c > hi[d] ) { hi[d] = c; }
  }
}

inline void growBox( const FixedVec_l& lo2 , const FixedVec_l& hi2 ,
		     FixedVec_l& lo , FixedVec_l& hi ) {
  for ( unsigned int d=0; d<DIMENSION; ++d ) {
    if ( lo2[d] < lo[d] ) { lo[d] = lo2[d]; }
    if ( hi2[d] > hi[d] ) { hi[d] = hi2[d]; }
  }
}

struct ThreadArgs
{
  NodeContainer * nodes;
  NodeInteractionHandler * nodeHandler;
  ParticleStats_t * stats;
  // The particles of the active box, packed for the repulsion phase
  // (shared by all threads, rebuilt every iteration)
  CellList * cells;
  // Voxels that can hold particles lie inside the box activeMin..activeMax
  // (shared by all threads). It is the union of staticMin..staticMax, the
  // particles placed before the simulation (root, anchors, -x), and the
  // per-thread boxMin..boxMax of the particles each thread integrated.
  // Threads take chunks of work through nextVoxel, so the work is
  // balanced dynamically.
  FixedVec_l boxMin, boxMax;
  FixedVec_l * activeMin, * activeMax;
  FixedVec_l * staticMin, * staticMax;
  std::atomic<long> * nextVoxel;
  double busySeconds; // time this thread spent working in calcInteractions
  unsigned long long pairTests; // pair tests in calcInteractions (profile)
  // All edges of the layout graph so far, and those touching the current
  // level (the only ones that count for the convergence statistic). Shared
  // by all threads, rebuilt when a level is added.
  const EdgeList * edges;
  const EdgeList * levelEdges;
  Grid_t * grid;
  prec_t eqDistance;
  EllipseFactors ellipseFactors;
  long threadCount;
  long whichThread;
  prec_t nbhdRadius;
  prec_t casualSpringConstant;
  prec_t specialSpringConstant;
  Graph_t * full_graph;
  Graph_t * layout_graph;
  LevelMap * levels;
  ParentMap * parents;
  unsigned int currentLevel;
};

//------------------------------------------------

void * packCells ( void * arg );
void * calcInteractions ( void * arg );
void * integrateParticles ( void * arg );
void * onlyEdgeInteractions( void * arg );
void * collectEdgeStats( void * arg ); 

//------------------------------------------------

int generateMSTFromNodes( NodeContainer& nodes ,  ThreadArgs * args , long p );

prec_t collectOutput ( ThreadArgs * args , PCChaperone& chaperone  );
void beginSimulation( ThreadContainer& threads , 
		      prec_t cutOffPrecision , TimeKeeper_t& timer, 
		      ThreadArgs * threadArgs , PCChaperone& chaperone ,
		      unsigned int totalLevels , bool b , prec_t placementDistance,
		      prec_t placementRadius , bool placeLeafsClose , bool silentOutput);

void adjustWeightsBasedOnChildrenCount( NodeContainer& nodes );
void solidifyLargeEdges( NodeContainer& nodes );

FixedVec_p calcCenterOfMass( Graph_t& g , NodeContainer& nodes ,  LevelMap& levels ,  unsigned int currentLevel);

void initializeCurrentLayer( Graph_t& g , NodeContainer& nodes , LevelMap& levels ,
			     ParentMap& p , Grid_t& grid , unsigned int currentLevel, Graph_t& full ,
			     prec_t placementDistance, prec_t placementRadius , bool placeLeafsClose );


long activeEdgeCount( Node& n );

prec_t  activateNextLayerOfEdges( NodeContainer& n , unsigned int currentLevel , ThreadArgs * threadArgs );
prec_t getBufferDistance( Node& n );
void printOutput(long i,prec_t d, long ll, ostream& o );

void layerNPlacement( NodeContainer& nodes , Grid_t& grid , out_graph& g , FixedVec_p& cm , unsigned int currentLevel , ParentMap& parents , Graph_t& full, LevelMap& lm , prec_t placementDistance,
		      prec_t placementRadius , bool placeLeafsClose );
void generatePlacementVector( FixedVec_p& d , const FixedVec_p& parentNode,
			      const FixedVec_p& parentParentNode , const FixedVec_p& cm );
prec_t placementFormula( prec_t placementDistance , int vertices2place, int dimension );

void gridPrepAndInit( NodeContainer& nc, Grid_t& g , prec_t voxelLength );
bool doesVertexHaveAnyChildren( Graph_t& G, Graph_t::vertex_descriptor v , out_graph& g , ParentMap& parents );

EllipseFactors parseEllipseFactors( const std::string& optionStr );

// Attempts to initialize any particle positions that are not initialized yet, by way of interpolation from its neighbors which have initialized positions already, if any.
// This is done by picking the "center point" of those neighbors' positions to be used as the initial position for a particle that wasn't initialized before.
// The process continues until either all particles have initial positions properly set, or no further progress can be made (due to isolated and totally uninitialized islands).
// Also, the progress of the stages and the final accomplishment is printed to stdout.
void interpolateUninitializedPositions( PCChaperone& chaperone, const Graph_t::boost_graph& g, bool remove_disconnected_nodes );

//------------------------------------------------

#endif

