/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMSURFACEPRESSURE_H
#define XFEMSURFACEPRESSURE_H

// Moose Includes
#include "DiracKernel.h"

class XFEM;
class Function;

class XFEMSurfacePressure : public DiracKernel
{
public:
  XFEMSurfacePressure(const InputParameters & parameters);

  virtual void addPoints();
  virtual Real computeQpResidual();

protected:
  void qRuleOnLine(Point & p1, Point & p2, std::vector<Point> & quad_pts, std::vector<Real> & quad_wts);
  const int _component;
  const Real _factor;  
  Function * const _function;

  std::map<const Elem*, std::map<Point, Point> > _elem_to_point_to_normal;
  std::map<const Elem*, std::map<Point, Real> > _elem_to_point_to_quadrature_weights;
  std::map<const Elem*, std::map<Point, Real> > _elem_to_point_to_integration_jacobian;

private:
  XFEM *_xfem;
};

template<>
InputParameters validParams<XFEMSurfacePressure>();

#endif //XFEMSURFACEPRESSURE_H
