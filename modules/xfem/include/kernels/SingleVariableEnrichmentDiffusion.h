/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef SINGLEVARIABLEENRICHMENTDIFFUSION_H
#define SINGLEVARIABLEENRICHMENTDIFFUSION_H

#include "Kernel.h"
#include "CrackFrontRTheta.h"

class SingleVariableEnrichmentDiffusion;

template<>
InputParameters validParams<SingleVariableEnrichmentDiffusion>();

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class SingleVariableEnrichmentDiffusion : public Kernel
{
public:
  SingleVariableEnrichmentDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
};

#endif /* SINGLEVARIABLEENRICHMENTDIFFUSION_H */
