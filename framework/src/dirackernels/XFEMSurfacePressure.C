/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "XFEMSurfacePressure.h"
#include "Function.h"
#include "XFEM.h"
#include "XFEMMiscFuncs.h"

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
XFEMSurfacePressure::qRuleOnLine(std::vector<Point> & intersection_points, std::vector<Point> & quad_pts, std::vector<Real> & quad_wts)
{
  Point p1 = intersection_points[0];
  Point p2 = intersection_points[1];

  //number of quadrature points
  unsigned int num_qpoints = 2;

  //quadrature coordinates
  Real xi0 = -std::sqrt(1.0/3.0);
  Real xi1 =  std::sqrt(1.0/3.0);
  
  quad_wts.resize(num_qpoints);
  quad_pts.resize(num_qpoints);

  Real integ_jacobian =  pow((p1 -  p2).size_sq(), 0.5) * 0.5;

  quad_wts[0] = 1.0 * integ_jacobian;
  quad_wts[1] = 1.0 * integ_jacobian;

  quad_pts[0] = (1.0 - xi0) / 2.0 * p1 + (1.0 + xi0) / 2.0 * p2;
  quad_pts[1] = (1.0 - xi1) / 2.0 * p1 + (1.0 + xi1) / 2.0 * p2;
}

void 
XFEMSurfacePressure::qRuleOnSurface(std::vector<Point> & intersection_points, std::vector<Point> & quad_pts, std::vector<Real> & quad_wts)
{
  unsigned nnd_pe = intersection_points.size();
  Point xcrd(0.0, 0.0, 0.0);
  for (unsigned int i = 0; i < intersection_points.size(); ++i)
    xcrd += intersection_points[i];
  xcrd *= (1.0/intersection_points.size());
  
  quad_pts.resize(nnd_pe);
  quad_wts.resize(nnd_pe);

  Real jac = 0.0;

  for (unsigned int j = 0; j < nnd_pe; ++j) // loop all sub-trigs
  {
    std::vector<std::vector<Real> > shape(3, std::vector<Real>(3,0.0));
    std::vector<Point> subtrig_points(3, Point(0.0,0.0,0.0)); // sub-trig nodal coords

    int jplus1(j < nnd_pe-1 ? j+1 : 0);
    subtrig_points[0] = xcrd;
    subtrig_points[1] = intersection_points[j];
    subtrig_points[2] = intersection_points[jplus1];

    std::vector<std::vector<Real> > sg2;
    stdQuadr2D(3, 1, sg2); // get sg2
    for (unsigned int l = 0; l < sg2.size(); ++l) // loop all int pts on a sub-trig
    {
      shapeFunc2D(3, sg2[l], subtrig_points, shape, jac, true); // Get shape
      std::vector<Real> tsg_line(3,0.0);
      for (unsigned int k = 0; k < 3; ++k) // loop sub-trig nodes
      {
        tsg_line[0] += shape[k][2] * subtrig_points[k](0);
        tsg_line[1] += shape[k][2] * subtrig_points[k](1);
        tsg_line[2] += shape[k][2] * subtrig_points[k](2);
      }
      quad_pts[j + l] = Point(tsg_line[0], tsg_line[1], tsg_line[2]);
      quad_wts[j + l] = sg2[l][3] * jac;
    }
  }
}

void
XFEMSurfacePressure::addPoints()
{
  std::vector<const Elem *> cutted_elems;
  _xfem->getCuttedElements(cutted_elems);

  _elem_to_point_to_normal.clear();
  _elem_to_point_to_quadrature_weights.clear();

  std::vector<Point> quad_pts;
  std::vector<Real> quad_wts;

  for (unsigned int i = 0; i < cutted_elems.size(); i++)
  {
    const Elem * elem = cutted_elems[i];

    quad_pts.clear();
    quad_wts.clear();

    std::vector<Point> intersectionPoints;
    Point normal(0.0, 0.0, 0.0);
    _xfem->get_intersection_info(elem, 0, normal, intersectionPoints);

    if (intersectionPoints.size() == 2)
      qRuleOnLine(intersectionPoints, quad_pts, quad_wts);
    else
      qRuleOnSurface(intersectionPoints, quad_pts, quad_wts);

    for (unsigned int j = 0; j < quad_pts.size(); ++j)
    {
      _elem_to_point_to_normal[elem][quad_pts[j]] = -1.0 * normal;

      _elem_to_point_to_quadrature_weights[elem][quad_pts[j]] = quad_wts[j];

      addPoint(elem, quad_pts[j]);
    }
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

  factor *= quad_wts;

  //std::cout << "factor = " << factor << ", normal = " << normal << std::endl;

  return -1.0 * factor * (normal(_component) * _test[_i][_qp]);
}
