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

#ifndef RadialDispAux_H
#define RadialDispAux_H

#include "AuxKernel.h"

//Forward Declarations
class RadialDispAux;

template<>
InputParameters validParams<RadialDispAux>();

/**
 * Constant auxiliary value
 */
class RadialDispAux : public AuxKernel
{
public:
  RadialDispAux(const InputParameters & parameters);

  virtual ~RadialDispAux() {}

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue();

  std::vector<Real> _center;
  const VariableValue & _disp_x;
  const VariableValue & _disp_y;



};

#endif //RadialDispAux_H
