/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMMATERIALTENSORMARKERUSEROBJECT_H
#define XFEMMATERIALTENSORMARKERUSEROBJECT_H

#include "XFEMMarkerUserObject.h"
#include "RankTwoTensor.h"

class XFEMMaterialTensorMarkerUserObject;

template <>
InputParameters validParams<XFEMMaterialTensorMarkerUserObject>();

class XFEMMaterialTensorMarkerUserObject : public XFEMMarkerUserObject
{
public:
  XFEMMaterialTensorMarkerUserObject(const InputParameters & parameters);
  virtual ~XFEMMaterialTensorMarkerUserObject() {}

protected:
  bool _use_weibull;
  const MaterialProperty<RankTwoTensor> & _tensor;

  MooseEnum _scalar_type;
  const Point _point1;
  const Point _point2;

  Real _threshold;
  bool _average;
  Real _random_range;

  const MaterialProperty<Real> & _weibull;

  virtual bool doesElementCrack(Point & direction);
};

#endif // XFEMMATERIALTENSORMARKERUSEROBJECT_H
