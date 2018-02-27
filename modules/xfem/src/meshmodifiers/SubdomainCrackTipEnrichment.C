//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SubdomainCrackTipEnrichment.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<SubdomainCrackTipEnrichment>()
{
  InputParameters params = validParams<MeshModifier>();
  params.addClassDescription(
      "Changes the subdomain ID of elements to inside and outside the enrichment.");
  params.addRequiredParam<SubdomainID>("block_id_inside_enrichment",
                                       "Subdomain id to set for inside enrichment");
  params.addRequiredParam<SubdomainID>("block_id_outside_enrichment",
                                       "Subdomain id to set for outside enrichment");
  params.addRequiredParam<std::vector<Point>>("crack_front_points",
                                              "Set of points to define crack front");
  params.addRequiredParam<Real>("radius", "The radius of the enrichment domain");
  return params;
}

SubdomainCrackTipEnrichment::SubdomainCrackTipEnrichment(const InputParameters & parameters)
  : MeshModifier(parameters),
    _block_id_inside_enrichment(parameters.get<SubdomainID>("block_id_inside_enrichment")),
    _block_id_outside_enrichment(parameters.get<SubdomainID>("block_id_outside_enrichment")),
    _crack_front_points(parameters.get<std::vector<Point>>("crack_front_points")),
    _radius(parameters.get<Real>("radius"))
{
}

void
SubdomainCrackTipEnrichment::modify()
{
  // Check that we have access to the mesh
  if (!_mesh_ptr)
    mooseError(
        "_mesh_ptr must be initialized before calling SubdomainCrackTipEnrichment::modify()");

  // Loop over the elements
  for (const auto & elem : _mesh_ptr->getMesh().active_element_ptr_range())
  {
    Point elem_center = elem->centroid();
    Real min_distance = (_crack_front_points[0] - elem_center).norm();

    for (unsigned int i = 1; i < _crack_front_points.size(); ++i)
    {
      Real distance = (_crack_front_points[i] - elem_center).norm();
      if (distance < min_distance)
        min_distance = distance;
    }

    if (min_distance < _radius)
      elem->subdomain_id() = _block_id_inside_enrichment;
    else
      elem->subdomain_id() = _block_id_outside_enrichment;
  }
}
