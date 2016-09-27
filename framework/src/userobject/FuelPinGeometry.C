/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

#include "FuelPinGeometry.h"
#include "MooseMesh.h"

namespace {

// This struct makes using a std::map easy because default constructed value_types automatically
// get the correct initial values for min and max when we add a new boundary ID.
struct MinMax {
  Real min;
  Real max;
  MinMax()
  : min( std::numeric_limits< Real >::infinity() ),
    max( -std::numeric_limits< Real >::infinity() ) {}
};
typedef std::map< BoundaryID, MinMax > BoundaryExtrema;

MinMax findBoundaryExtremum(
  const MooseMesh & mesh,
  const BoundaryName & name,
  const unsigned dimension,
  const std::vector< BoundaryExtrema > & extrema )
{
  const BoundaryID id = mesh.getBoundaryID( name );
  BoundaryExtrema::const_iterator found = extrema[dimension].find( id );
  if ( found == extrema[dimension].end() ) {
    mooseError( "Boundary \"" << name << "\" does not exist for this subdomain." );
  }
  return found->second;
}

} // namespace anon

template<>
InputParameters validParams<FuelPinGeometry>()
{
  InputParameters params = validParams<GeneralUserObject>();

  // These default sideset numbers are from the BISON mesh script.
  params.addParam<BoundaryName>( "clad_inner_wall", "5", "Sideset for inner wall of cladding, not including end caps.");
  params.addParam<BoundaryName>( "clad_outer_wall", "2", "Sideset for outer wall of cladding.");
  params.addParam<BoundaryName>( "clad_top", "3", "Sideset for top of cladding (top of upper end cap).");
  params.addParam<BoundaryName>( "clad_bottom", "1", "Sideset for bottom of cladding (bottom of lower end cap).");
  params.addParam<BoundaryName>( "pellet_exteriors", "8", "Sideset for all pellet exteriors.");

  params.addClassDescription( "Calculates LWR fuel pin geometry by reading the input mesh. This"
    " object can be coupled to Burnup and other functions as an alternative to having the user"
    " supply parameters such as pellet radius and pellet-cladding gap." );

  params.set<bool>( "use_displaced_mesh" ) = false;
  return params;
}

FuelPinGeometry::FuelPinGeometry(const InputParameters & parameters)
: GeneralUserObject(parameters),
  _mesh(_subproblem.mesh())
{
  const BoundaryName clad_inner_wall = getParam<BoundaryName>("clad_inner_wall");
  const BoundaryName clad_outer_wall = getParam<BoundaryName>("clad_outer_wall");
  const BoundaryName clad_top = getParam<BoundaryName>("clad_top");
  const BoundaryName clad_bottom = getParam<BoundaryName>("clad_bottom");
  const BoundaryName pellet_exteriors = getParam<BoundaryName>("pellet_exteriors");

  // For now, make sure everything is RZ.
  std::set< SubdomainID > subdomains = getSubProblem().mesh().meshSubdomains();
  for ( std::set< SubdomainID >::const_iterator it = subdomains.begin(); it != subdomains.end(); ++it )
  {
    Moose::CoordinateSystemType coord_system = _subproblem.getCoordSystem( *it );
    if ( coord_system != Moose::COORD_RZ ) {
      mooseError( "CoordinateSystemType " << coord_system << " is not supported by this object." );
    }
  }
  const unsigned dimensions = 2;

  // libmesh has a bounding box function that would work for the overall pin dimensions and the fuel
  // stack dimensions, but we would still be stuck iterating through nodes to find the inner surface
  // of the cladding. So long as we are iterating through all boundary nodes, we may as well do all
  // of the work here.

  // Find the bounding boxes for all boundaries by just checking the coordinates of the nodes on
  // the boundaries.
  std::vector< BoundaryExtrema > extrema( dimensions );
  ConstBndNodeRange & all_boundary_nodes = *_mesh.getBoundaryNodeRange();
  for ( ConstBndNodeRange::const_iterator it = all_boundary_nodes.begin(); it != all_boundary_nodes.end(); ++it )
  {
    const BoundaryID boundary = (*it)->_bnd_id;
    const Node & node = *(*it)->_node;
    for ( unsigned dim = 0; dim < dimensions; ++dim ) {
      extrema[dim][boundary].min = std::min( node(dim), extrema[dim][boundary].min );
      extrema[dim][boundary].max = std::max( node(dim), extrema[dim][boundary].max );
    }
  }

  // Get the inner and outer radius of the cladding from the sidesets.
  // Sideset 5 is cladding innner wall not including end caps.
  _clad_inner_radius = findBoundaryExtremum( _mesh, clad_inner_wall, 0, extrema ).min;

  // Sideset 2 is cladding exterior wall not including top/bottom.
  _clad_outer_radius = findBoundaryExtremum( _mesh, clad_outer_wall, 0, extrema ).max;

  // Sideset 3 is top of cladding.
  _top_of_rod = findBoundaryExtremum( _mesh, clad_top, 1, extrema ).max;

  // Sideset 1 is bottom of cladding.
  _bottom_of_rod = findBoundaryExtremum( _mesh, clad_bottom, 1, extrema ).min;

  // Sideset 8 includes all of the pellet exteriors.
  _fuel_outer_radius = findBoundaryExtremum( _mesh, pellet_exteriors, 0, extrema ).max;
  _fuel_inner_radius = findBoundaryExtremum( _mesh, pellet_exteriors, 0, extrema ).min;
  _top_of_stack = findBoundaryExtremum( _mesh, pellet_exteriors, 1, extrema ).max;
  _bottom_of_stack = findBoundaryExtremum( _mesh, pellet_exteriors, 1, extrema ).min;

  Moose::out
    << "\nFuel pin geometry read from mesh:"
    << "\n  Fuel rod height:      " << _top_of_rod - _bottom_of_rod
    << "\n  Cladding OD:          " << 2 * _clad_outer_radius
    << "\n  Cladding thickness:   " << _clad_outer_radius - _clad_inner_radius
    << "\n  Pellet-clad gap:      " << _clad_inner_radius - _fuel_outer_radius
    << "\n  Fuel ID:              " << 2 * _fuel_inner_radius
    << "\n  Fuel OD:              " << 2 * _fuel_outer_radius
    << "\n  Top of fuel stack:    " << _top_of_stack
    << "\n  Bottom of fuel stack: " << _bottom_of_stack << "\n" << std::endl;
}

void FuelPinGeometry::initialize()
{
}

void FuelPinGeometry::execute()
{
}

void FuelPinGeometry::finalize()
{
}

Real FuelPinGeometry::pellet_OD() const
{
  return 2 * _fuel_outer_radius;
}

Real FuelPinGeometry::pellet_ID() const
{
  return 2 * _fuel_inner_radius;
}

Real FuelPinGeometry::top_of_stack() const
{
  return _top_of_stack;
}

Real FuelPinGeometry::bottom_of_stack() const
{
  return _bottom_of_stack;
}

Real FuelPinGeometry::clad_OD() const
{
  return 2 * _clad_outer_radius;
}

Real FuelPinGeometry::clad_ID() const
{
  return 2 * _clad_inner_radius;
}

Real FuelPinGeometry::top_of_rod() const
{
  return _top_of_rod;
}

Real FuelPinGeometry::bottom_of_rod() const
{
  return _bottom_of_rod;
}
