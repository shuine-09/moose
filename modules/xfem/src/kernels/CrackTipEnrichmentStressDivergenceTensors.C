/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CrackTipEnrichmentStressDivergenceTensors.h"
#include "ElasticityTensorTools.h"

template <>
InputParameters
validParams<CrackTipEnrichmentStressDivergenceTensors>()
{
  InputParameters params = validParams<ALEKernel>();
  params.addClassDescription("Enrich stress divergence kernel for the Cartesian coordinate system");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction the variable "
                                        "this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<unsigned int>("enrichment_component",
                                        "The component of the enrichement functions");
  params.addRequiredCoupledVar("enrichment_displacement",
                               "The string of displacements suitable for the problem statement");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredParam<UserObjectName>("crack_front_definition",
                                          "The CrackFrontDefinition user object name");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

CrackTipEnrichmentStressDivergenceTensors::CrackTipEnrichmentStressDivergenceTensors(
    const InputParameters & parameters)
  : ALEKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _component(getParam<unsigned int>("component")),
    _enrichment_component(getParam<unsigned int>("enrichment_component")),
    _nenrich_disp(coupledComponents("enrichment_displacement")),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _B(4),
    _dBX(4),
    _dBx(4)
{
  _enrich_disp_var.resize(_nenrich_disp);
  for (unsigned int i = 0; i < _nenrich_disp; ++i)
    _enrich_disp_var[i] = coupled("enrichment_displacement", i);

  if (_nenrich_disp == 8)
    _BI.resize(4); // QUAD4
  else if (_nenrich_disp == 12)
    _BI.resize(8); // HEX8

  for (unsigned int i = 0; i < _BI.size(); ++i)
    _BI[i].resize(4);
}

void
CrackTipEnrichmentStressDivergenceTensors::prepareCrackTipEnrichementFunctionAtNode()
{
  for (unsigned int i = 0; i < _BI.size(); ++i)
  {
    _crack_front_definition->calculateRThetaToCrackFront(*(_current_elem->get_node(i)), _r, _theta);

    Real st = std::sin(_theta);
    Real st2 = std::sin(_theta / 2.0);
    Real ct2 = std::cos(_theta / 2.0);
    Real sr = std::sqrt(_r);

    _BI[i][0] = sr * st2;
    _BI[i][1] = sr * ct2;
    _BI[i][2] = sr * st2 * st;
    _BI[i][3] = sr * ct2 * st;
  }
}

Real
CrackTipEnrichmentStressDivergenceTensors::computeQpResidual()
{
  prepareCrackTipEnrichementFunctionAtNode();

  unsigned int crack_front_point_index =
      _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], _r, _theta);

  Real st = std::sin(_theta);
  Real ct = std::cos(_theta);
  Real st2 = std::sin(_theta / 2.0);
  Real ct2 = std::cos(_theta / 2.0);
  Real st15 = std::sin(1.5 * _theta);
  Real ct15 = std::cos(1.5 * _theta);
  Real sr = std::sqrt(_r);

  _B[0] = sr * st2;
  _B[1] = sr * ct2;
  _B[2] = sr * st2 * st;
  _B[3] = sr * ct2 * st;

  _dBx[0](0) = -0.5 / sr * st2;
  _dBx[0](1) = 0.5 / sr * ct2;
  _dBx[0](2) = 0.0;
  _dBx[1](0) = 0.5 / sr * ct2;
  _dBx[1](1) = 0.5 / sr * st2;
  _dBx[1](2) = 0.0;
  _dBx[2](0) = -0.5 / sr * st15 * st;
  _dBx[2](1) = 0.5 / sr * (st2 + st15 * ct);
  _dBx[2](2) = 0.0;
  _dBx[3](0) = -0.5 / sr * ct15 * st;
  _dBx[3](1) = 0.5 / sr * (ct2 + ct15 * ct);
  _dBx[3](2) = 0.0;

  for (unsigned int i = 0; i < 4; ++i)
    _dBX[i] = _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(_dBx[i],
                                                                          crack_front_point_index);

  RealVectorValue grad_B(_dBX[_enrichment_component]);

  return _stress[_qp].row(_component) *
         (_grad_test[_i][_qp] * (_B[_enrichment_component] - _BI[_i][_enrichment_component]) +
          _test[_i][_qp] * grad_B);
}

Real
CrackTipEnrichmentStressDivergenceTensors::computeQpJacobian()
{
  prepareCrackTipEnrichementFunctionAtNode();

  unsigned int crack_front_point_index =
      _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], _r, _theta);

  Real st = std::sin(_theta);
  Real ct = std::cos(_theta);
  Real st2 = std::sin(_theta / 2.0);
  Real ct2 = std::cos(_theta / 2.0);
  Real st15 = std::sin(1.5 * _theta);
  Real ct15 = std::cos(1.5 * _theta);
  Real sr = std::sqrt(_r);

  _B[0] = sr * st2;
  _B[1] = sr * ct2;
  _B[2] = sr * st2 * st;
  _B[3] = sr * ct2 * st;

  _dBx[0](0) = -0.5 / sr * st2;
  _dBx[0](1) = 0.5 / sr * ct2;
  _dBx[0](2) = 0.0;
  _dBx[1](0) = 0.5 / sr * ct2;
  _dBx[1](1) = 0.5 / sr * st2;
  _dBx[1](2) = 0.0;
  _dBx[2](0) = -0.5 / sr * st15 * st;
  _dBx[2](1) = 0.5 / sr * (st2 + st15 * ct);
  _dBx[2](2) = 0.0;
  _dBx[3](0) = -0.5 / sr * ct15 * st;
  _dBx[3](1) = 0.5 / sr * (ct2 + ct15 * ct);
  _dBx[3](2) = 0.0;

  for (unsigned int i = 0; i < 4; ++i)
    _dBX[i] = _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(_dBx[i],
                                                                          crack_front_point_index);

  RealVectorValue grad_B(_dBX[_enrichment_component]);

  RealVectorValue grad_test =
      _grad_test[_i][_qp] * (_B[_enrichment_component] - _BI[_i][_enrichment_component]) +
      _test[_i][_qp] * grad_B;
  RealVectorValue grad_phi =
      _grad_phi[_j][_qp] * (_B[_enrichment_component] - _BI[_j][_enrichment_component]) +
      _phi[_j][_qp] * grad_B;

  return ElasticityTensorTools::elasticJacobian(
      _Jacobian_mult[_qp], _component, _component, grad_test, grad_phi);
}

Real
CrackTipEnrichmentStressDivergenceTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component = 0;
  unsigned int coupled_enrichment_component = 0;
  bool active(false);

  for (unsigned int i = 0; i < _enrich_disp_var.size(); ++i)
  {
    if (jvar == _enrich_disp_var[i])
    {
      coupled_component = i / 4;
      coupled_enrichment_component = i % 4;
      active = true;
    }
  }

  if (active)
  {
    prepareCrackTipEnrichementFunctionAtNode();

    unsigned int crack_front_point_index =
        _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], _r, _theta);

    Real st = std::sin(_theta);
    Real ct = std::cos(_theta);
    Real st2 = std::sin(_theta / 2.0);
    Real ct2 = std::cos(_theta / 2.0);
    Real st15 = std::sin(1.5 * _theta);
    Real ct15 = std::cos(1.5 * _theta);
    Real sr = std::sqrt(_r);

    _B[0] = sr * st2;
    _B[1] = sr * ct2;
    _B[2] = sr * st2 * st;
    _B[3] = sr * ct2 * st;

    _dBx[0](0) = -0.5 / sr * st2;
    _dBx[0](1) = 0.5 / sr * ct2;
    _dBx[0](2) = 0.0;
    _dBx[1](0) = 0.5 / sr * ct2;
    _dBx[1](1) = 0.5 / sr * st2;
    _dBx[1](2) = 0.0;
    _dBx[2](0) = -0.5 / sr * st15 * st;
    _dBx[2](1) = 0.5 / sr * (st2 + st15 * ct);
    _dBx[2](2) = 0.0;
    _dBx[3](0) = -0.5 / sr * ct15 * st;
    _dBx[3](1) = 0.5 / sr * (ct2 + ct15 * ct);
    _dBx[3](2) = 0.0;

    for (unsigned int i = 0; i < 4; ++i)
      _dBX[i] = _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(
          _dBx[i], crack_front_point_index);

    RealVectorValue grad_B_test(_dBX[_enrichment_component]);
    RealVectorValue grad_B_phi(_dBX[coupled_enrichment_component]);

    RealVectorValue grad_test =
        _grad_test[_i][_qp] * (_B[_enrichment_component] - _BI[_i][_enrichment_component]) +
        _test[_i][_qp] * grad_B_test;
    RealVectorValue grad_phi = _grad_phi[_j][_qp] * (_B[coupled_enrichment_component] -
                                                     _BI[_j][coupled_enrichment_component]) +
                               _phi[_j][_qp] * grad_B_phi;

    return ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], _component, coupled_component, grad_test, grad_phi);
  }

  return 0;
}
