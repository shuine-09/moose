//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SUBDOMAINCRACKTIPENRICHMENT_H
#define SUBDOMAINCRACKTIPENRICHMENT_H

// MOOSE includes
#include "MooseEnum.h"
#include "MeshModifier.h"

// Forward declerations
class SubdomainCrackTipEnrichment;

template <>
InputParameters validParams<SubdomainCrackTipEnrichment>();

/**
 * MeshModifier for defining a Subdomain that uses for crack tip enrichment
 */
class SubdomainCrackTipEnrichment : public MeshModifier
{
public:
  /**
   * Class constructor
   * @param parameters The input parameters
   */
  SubdomainCrackTipEnrichment(const InputParameters & parameters);

  virtual void modify() override;

private:
  /// Block ID to assign to the region
  SubdomainID _block_id_inside_enrichment;
  SubdomainID _block_id_outside_enrichment;
  std::vector<Point> _crack_front_points;
  Real _radius;
};

#endif // SUBDOMAINCRACKTIPENRICHMENT_H
