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

#ifndef XFEMEQUALVALUECONSTRAINT_H
#define XFEMEQUALVALUECONSTRAINT_H

// MOOSE includes
#include "XFEMElementConstraint.h"
#include "MooseMesh.h"

// Forward Declarations
class XFEMEqualValueConstraint;

template<>
InputParameters validParams<XFEMEqualValueConstraint>();

class XFEMEqualValueConstraint : public XFEMElementConstraint 
{
public:
  XFEMEqualValueConstraint(const InputParameters & parameters);
  virtual ~XFEMEqualValueConstraint();

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);
};

#endif /* XFEMEQUALVALUECONSTRAINT_H_ */
