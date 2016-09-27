/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

#include "ThermalMAMOX.h"
template<>
InputParameters validParams<ThermalMAMOX>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("temp", "Coupled Temperature");
  params.addCoupledVar("porosity", "Coupled Porosity from field variable");
  params.addParam<FunctionName>("fp_function", "", "porosity reducing thermal conductivity function");
  params.addParam<Real>("oxy_to_metal_ratio", 2.0, "Oxygen to metal ratio");
  params.addParam<Real>("Am_content", 0.00, "Weight fraction of Am in MAMOX fuel");
  params.addParam<Real>("Np_content", 0.00, "Weight fraction of Np in MAMOX fuel");
  return params;
}

ThermalMAMOX::ThermalMAMOX(const InputParameters & parameters)
  :Material(parameters),
  _has_temp(isCoupled("temp")),
  _temp(_has_temp ? coupledValue("temp") : _zero),
  _has_porosity(isCoupled("porosity")),
  _porosity(_has_porosity ? coupledValue("porosity") : _zero),
  _fp_function(getParam<FunctionName>("fp_function") == "" ? NULL : &getFunction("fp_function")),
  _thermal_conductivity(declareProperty<Real>("thermal_conductivity")),
  _thermal_conductivity_dT(declareProperty<Real>("thermal_conductivity_dT")),
  _specific_heat(declareProperty<Real>("specific_heat")),
  _stoech_dev(2.0 - getParam<Real>("oxy_to_metal_ratio")),  //deviation from stoechiometry
  _Am_content(getParam<Real>("Am_content")),
  _Np_content(getParam<Real>("Np_content"))
{
}

void
ThermalMAMOX::computeQpProperties()
{
  //    Moose::out << "porosity = " << _porosity[_qp] << std::endl;
  //    Moose::out << "stoech_dev = " << _stoech_dev << std::endl;
  //    Moose::out << "Am = " << _Am_content << std::endl;
  //    Moose::out << "Np = " << _Np_content << std::endl;
  Real fp(1);
  if(_fp_function)
    fp = _fp_function->value(_porosity[_qp], _q_point[_qp]);

  Real term1 = (1-_porosity[_qp])/(1+0.5*_porosity[_qp]);
  Real term2 = 1/((2.713*_stoech_dev+3.583e-1*_Am_content+6.317e-2*_Np_content+1.595e-2)+(-2.625*_stoech_dev+2.493)*1e-4*_temp[_qp]);
  Real term3 = (1.541e11/pow(_temp[_qp],2.5))*exp(-1.522e4/_temp[_qp]);

  if(_fp_function)
    _thermal_conductivity[_qp] = fp*term2*term3;
  else
    _thermal_conductivity[_qp] = term1*term2+term3;
  //    Moose::out << "term1 = " << term1 << std::endl;
  //    Moose::out << "term2 = " << term2 << std::endl;
  //    Moose::out << "term3 = " << term3 << std::endl;
  //    _thermal_conductivity[_qp] = (1/((2.713*_stoech_dev+3.583e-1*_Am_content+6.317e-2*_Np_content+1.595e-2)+(-2.625*_stoech_dev+2.493)*1e-4*_temp[_qp])+((1.541e11/pow(_temp[_qp],2.5))*exp(-1.522e4/_temp[_qp])))*fp;
  //Journal of NUCLEAR SCIENCE and TECHNOLOGY, Vol. 48, No. 4, p. 646-653 (2011)
  _thermal_conductivity_dT[_qp] =0;
  _specific_heat[_qp] = 120;// placeholder, not needed for steady state analysis
}
