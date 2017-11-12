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

#include "LSMForcingFunction.h"

registerMooseObject("MooseApp", LSMForcingFunction);

InputParameters
LSMForcingFunction::validParams()
{
  InputParameters params = Kernel::validParams();
  return params;
}

LSMForcingFunction::LSMForcingFunction(const InputParameters & parameters) :
    Kernel(parameters)
{
}

Real
LSMForcingFunction::f()
{
  Real x = _q_point[_qp](0);
  Real y = _q_point[_qp](1);

  double epsS = 0.05;
  double epsM = 6.0e-5;
  double PI = libMesh::pi;
  double PI2 = pow(libMesh::pi,2.0);
  double PI4 = pow(libMesh::pi,4.0);
  double cos2 = std::cos((PI*x)/3.0)*std::cos((PI*x)/3.0);
  double sin2 = std::sin((PI*x)/3.0)*std::sin((PI*x)/3.0);
  double cos22 = std::cos(PI*y)*std::cos(PI*y);
  double sin22 = std::sin(PI*y)*std::sin(PI*y);
  double dpsidx=(2.0*PI*std::cos((PI*x)/3.0)*std::sin((PI*x)/3.0)*sin22)/3.0;
  double d2psid2x = (2.0*PI2*cos2*sin22)/9.0 - (2.0*PI2*sin2*sin22)/9.0;
  double d2psid2y = 2.0*PI2*cos22*sin2 - 2.0*PI2*sin2*sin22;
  double Lpsi = d2psid2x+d2psid2y;
  double L2psi = (8.0*PI4*cos2*cos22)/9.0 - (80.0*PI4*cos2*sin22)/81.0- (80.0*PI4*cos22*sin2)/9.0 + (728.0*PI4*sin2*sin22)/81.0;

  return epsS*Lpsi - epsM*L2psi + dpsidx;
}

Real
LSMForcingFunction::computeQpResidual()
{
  return _test[_i][_qp] * f();
}
