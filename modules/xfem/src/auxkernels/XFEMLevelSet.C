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
  params.addRequiredParam<FunctionName>("function", "The level set function.");
  return params;
}

XFEMLevelSet::XFEMLevelSet(const InputParameters & parameters) :
    AuxKernel(parameters),
    _func(getFunction("function"))
{
}

Real
XFEMLevelSet::computeValue()
{
  if (!isNodal())
    mooseError("XFEMLevelSet only supports Nodal AuxVariable");
  
  return _func.value(_t,*_current_node);
}

