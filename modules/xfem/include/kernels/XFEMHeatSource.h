/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMHEATSOURCE_H
#define XFEMHEATSOURCE_H

#include "BodyForce.h"

//Forward Declarations
class XFEMHeatSource;

template<>
InputParameters validParams<XFEMHeatSource>();

class XFEMHeatSource : public BodyForce
{
public:
  XFEMHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

};

#endif
