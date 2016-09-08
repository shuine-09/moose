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

#include "XFEMLevelSet.h"
#include "Function.h"

template<>
InputParameters validParams<XFEMLevelSet>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

XFEMLevelSet::XFEMLevelSet(const InputParameters & parameters) :
    AuxKernel(parameters)
{
}

Real
XFEMLevelSet::computeValue()
{
  Real x = (*_current_node)(0);
  Real y = (*_current_node)(1);

  if (isNodal())
  {
    //return std::sqrt(x * x + y * y) - 0.5;
    return x - 0.5;
  }
  else
    mooseError("XFEMLevelSet only supports Nodal AuxVariable");
}

