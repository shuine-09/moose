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

#include "RadialDispAux.h"

template<>
InputParameters validParams<RadialDispAux>()
{
  InputParameters params = validParams<AuxKernel>();


  params.addParam<std::vector<Real> >("center","coordinates of the center of radial displacement");
  params.addRequiredCoupledVar("disp_x", "ddisplacement variable in x direction");
  params.addRequiredCoupledVar("disp_y", "ddisplacement variable in y direction");

  return params;
}

RadialDispAux::RadialDispAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _center(getParam<std::vector<Real> >("center")),
    _disp_x(coupledValue("disp_x")),
    _disp_y(coupledValue("disp_y"))

{
}

Real
RadialDispAux::computeValue()
{

     Point p = (*_current_node);

     double p_0 = p(0)-_center[0];
     double p_1 = p(1)-_center[1];


//compute radius 1
    double rad_1 = sqrt(p_0*p_0+p_1*p_1);
//compute radius 2
    double rad_2 = sqrt((p_0+_disp_x[_qp])*(p_0+_disp_x[_qp])+(p_1+_disp_y[_qp])*(p_1+_disp_y[_qp]));
//compute radial displacement
    double rad = rad_2 - rad_1;


  return  rad;
}
