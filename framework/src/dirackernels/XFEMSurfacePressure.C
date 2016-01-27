/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "XFEMSurfacePressure.h"
#include "Function.h"
#include "XFEM.h"

template<>
InputParameters validParams<XFEMSurfacePressure>()
{
  InputParameters params = validParams<DiracKernel>();
  params.addRequiredParam<unsigned int>("component", "The component for the pressure");
  params.addParam<Real>("factor", 1.0, "The magnitude to use in computing the pressure"); 
  params.addParam<FunctionName>("function", "The function that describes the pressure");
  return params;
}



XFEMSurfacePressure::XFEMSurfacePressure(const InputParameters & parameters) :
  DiracKernel(parameters),
  _component(getParam<unsigned int>("component")), 
  _factor(getParam<Real>("factor")),
  _function( isParamValid("function") ? &getFunction("function") : NULL )
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMVolFracAux");
  _xfem = fe_problem->get_xfem();
}

void
XFEMSurfacePressure::qRuleOnLine(Point & p1, Point & p2, std::vector<Point> & quad_pts, std::vector<Real> & quad_wts)
{
  //number of quadrature points
  unsigned int num_qpoints = 2;

  //quadrature coordinates
  Real xi0 = -std::sqrt(1.0/3.0);
  Real xi1 =  std::sqrt(1.0/3.0);
  
  quad_wts.resize(num_qpoints);
  quad_pts.resize(num_qpoints);

  quad_wts[0] = 1.0;
  quad_wts[1] = 1.0;

  quad_pts[0] = (1.0 - xi0) / 2.0 * p1 + (1.0 + xi0) / 2.0 * p2;
  quad_pts[1] = (1.0 - xi1) / 2.0 * p1 + (1.0 + xi1) / 2.0 * p2;
}

void
XFEMSurfacePressure::addPoints()
{
  std::vector<const Elem *> cutted_elems;
  _xfem->getCuttedElements(cutted_elems);

  _elem_to_point_to_normal.clear();
  _elem_to_point_to_quadrature_weights.clear();
  _elem_to_point_to_integration_jacobian.clear();

  for (unsigned int i = 0; i < cutted_elems.size(); i++)
  {
    const Elem * elem = cutted_elems[i];

    std::vector<Point> intersectionPoints;
    Point normal(0.0, 0.0, 0.0);
    _xfem->get_intersection_info(elem, 0, normal, intersectionPoints);

    std::vector<Point> quad_pts;
    std::vector<Real> quad_wts;

    qRuleOnLine(intersectionPoints[0], intersectionPoints[1], quad_pts, quad_wts);
    
    Real integ_jacobian = pow((intersectionPoints[1] -  intersectionPoints[0]).size_sq(), 0.5) * 0.5;

    _elem_to_point_to_normal[elem][quad_pts[0]] = -1.0 * normal;
    _elem_to_point_to_normal[elem][quad_pts[1]] = -1.0 * normal;

    _elem_to_point_to_quadrature_weights[elem][quad_pts[0]] = quad_wts[0];
    _elem_to_point_to_quadrature_weights[elem][quad_pts[1]] = quad_wts[1];

    _elem_to_point_to_integration_jacobian[elem][quad_pts[0]] = integ_jacobian;
    _elem_to_point_to_integration_jacobian[elem][quad_pts[1]] = integ_jacobian;

    addPoint(elem, quad_pts[0]);
    addPoint(elem, quad_pts[1]);
  }
}

Real
XFEMSurfacePressure::computeQpResidual()
{
  Real factor = _factor;
  
  if (_function)
    factor *= _function->value(_t, _current_point);
 
  Point normal = _elem_to_point_to_normal[_current_elem][_current_point];
  Real quad_wts = _elem_to_point_to_quadrature_weights[_current_elem][_current_point];
  Real integ_jac = _elem_to_point_to_integration_jacobian[_current_elem][_current_point];

  factor *= quad_wts * integ_jac;

  return -1.0 * factor * (normal(_component) * _test[_i][_qp]);

}
