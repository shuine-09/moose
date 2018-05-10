//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMMOVINGINTERFACEVELOCITYBASE_H
#define XFEMMOVINGINTERFACEVELOCITYBASE_H

#include "DiscreteElementUserObject.h"
#include "PointValueAtXFEMInterface.h"

class XFEMMovingInterfaceVelocityBase;

template <>
InputParameters validParams<XFEMMovingInterfaceVelocityBase>();

class XFEMMovingInterfaceVelocityBase : public DiscreteElementUserObject
{
public:
  XFEMMovingInterfaceVelocityBase(const InputParameters & parameters);
  virtual ~XFEMMovingInterfaceVelocityBase() {}

  virtual void initialize() override;

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const = 0;

  unsigned int numberPoints() const { return _value_at_interface_uo->numberPoints(); }

protected:
  const PointValueAtXFEMInterface * _value_at_interface_uo;
};

#endif // XFEMMOVINGINTERFACEVELOCITYBASE_H
