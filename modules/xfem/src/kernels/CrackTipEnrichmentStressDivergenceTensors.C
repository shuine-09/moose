/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CrackTipEnrichmentStressDivergenceTensors.h"
#include "Material.h"
#include "MooseMesh.h"
#include "ElasticityTensorTools.h"
#include "libmesh/quadrature.h"
#include "XFEM.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "MooseUtils.h"

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
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "enrichment_displacement", "The string of displacements suitable for the problem statement");
  params.addRequiredCoupledVar("enrichment_displacement_var",
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
    _nl_vnames(getParam<std::vector<NonlinearVariableName>>("enrichment_displacement")),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition"))
{
  _enrich_disp_var.resize(_nl_vnames.size());
  for (unsigned int i = 0; i < _nl_vnames.size(); ++i)
    _enrich_disp_var[i] = coupled("enrichment_displacement_var", i);

  FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);

  if (fe_problem == NULL)
    mooseError(
        "Problem casting _subproblem to FEProblem in CrackTipEnrichmentStressDivergenceTensors");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in CrackTipEnrichmentStressDivergenceTensors");
}

Real
CrackTipEnrichmentStressDivergenceTensors::computeQpResidual()
{
  // calculate the near-tip enrichement function
  std::vector<Real> B;
  std::vector<RealVectorValue> dBX, dBx;
  std::vector<std::vector<Real>> BI;

  B.resize(4);
  dBX.resize(4);
  dBx.resize(4);
  BI.resize(4);

  Real r, theta;

  for (unsigned int i = 0; i < BI.size(); ++i)
  {
    BI[i].resize(4);

    _crack_front_definition->calculateRThetaToCrackFront(
        *(_current_elem->get_node(i)), 0, r, theta);

    Real st = std::sin(theta);
    Real st2 = std::sin(theta / 2.0);
    Real ct2 = std::cos(theta / 2.0);
    Real sr = std::sqrt(r);

    BI[i][0] = sr * st2;
    BI[i][1] = sr * ct2;
    BI[i][2] = sr * st2 * st;
    BI[i][3] = sr * ct2 * st;
  }

  _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], 0, r, theta);

  Real st = std::sin(theta);
  Real ct = std::cos(theta);
  Real st2 = std::sin(theta / 2.0);
  Real ct2 = std::cos(theta / 2.0);
  Real st15 = std::sin(1.5 * theta);
  Real ct15 = std::cos(1.5 * theta);
  Real sr = std::sqrt(r);

  B[0] = sr * st2;
  B[1] = sr * ct2;
  B[2] = sr * st2 * st;
  B[3] = sr * ct2 * st;

  dBx[0](0) = -0.5 / sr * st2;
  dBx[0](1) = 0.5 / sr * ct2;
  dBx[0](2) = 0.0;
  dBx[1](0) = 0.5 / sr * ct2;
  dBx[1](1) = 0.5 / sr * st2;
  dBx[1](2) = 0.0;
  dBx[2](0) = -0.5 / sr * st15 * st;
  dBx[2](1) = 0.5 / sr * (st2 + st15 * ct);
  dBx[2](2) = 0.0;
  dBx[3](0) = -0.5 / sr * ct15 * st;
  dBx[3](1) = 0.5 / sr * (ct2 + ct15 * ct);
  dBx[3](2) = 0.0;

  for (unsigned int i = 0; i < 4; ++i)
    dBX[i] = _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dBx[i],
                                                                         0); // TODO: point index

  RealVectorValue grad_B(dBX[_enrichment_component]);

  return _stress[_qp].row(_component) *
         (_grad_test[_i][_qp] * (B[_enrichment_component] - BI[_i][_enrichment_component]) +
          _test[_i][_qp] * grad_B);
}

Real
CrackTipEnrichmentStressDivergenceTensors::computeQpJacobian()
{
  // calculate the near-tip enrichement function
  std::vector<Real> B;
  std::vector<RealVectorValue> dBX, dBx;
  std::vector<std::vector<Real>> BI;

  B.resize(4);
  dBX.resize(4);
  dBx.resize(4);
  BI.resize(4);

  Real r, theta;

  for (unsigned int i = 0; i < BI.size(); ++i)
  {
    BI[i].resize(4);

    _crack_front_definition->calculateRThetaToCrackFront(
        *(_current_elem->get_node(i)), 0, r, theta);

    Real st = std::sin(theta);
    Real st2 = std::sin(theta / 2.0);
    Real ct2 = std::cos(theta / 2.0);
    Real sr = std::sqrt(r);

    BI[i][0] = sr * st2;
    BI[i][1] = sr * ct2;
    BI[i][2] = sr * st2 * st;
    BI[i][3] = sr * ct2 * st;
  }

  _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], 0, r, theta);

  Real st = std::sin(theta);
  Real ct = std::cos(theta);
  Real st2 = std::sin(theta / 2.0);
  Real ct2 = std::cos(theta / 2.0);
  Real st15 = std::sin(1.5 * theta);
  Real ct15 = std::cos(1.5 * theta);
  Real sr = std::sqrt(r);

  B[0] = sr * st2;
  B[1] = sr * ct2;
  B[2] = sr * st2 * st;
  B[3] = sr * ct2 * st;

  dBx[0](0) = -0.5 / sr * st2;
  dBx[0](1) = 0.5 / sr * ct2;
  dBx[0](2) = 0.0;
  dBx[1](0) = 0.5 / sr * ct2;
  dBx[1](1) = 0.5 / sr * st2;
  dBx[1](2) = 0.0;
  dBx[2](0) = -0.5 / sr * st15 * st;
  dBx[2](1) = 0.5 / sr * (st2 + st15 * ct);
  dBx[2](2) = 0.0;
  dBx[3](0) = -0.5 / sr * ct15 * st;
  dBx[3](1) = 0.5 / sr * (ct2 + ct15 * ct);
  dBx[3](2) = 0.0;

  for (unsigned int i = 0; i < 4; ++i)
    dBX[i] = _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dBx[i],
                                                                         0); // TODO: point index

  RealVectorValue grad_B(dBX[_enrichment_component]);

  RealVectorValue grad_test =
      _grad_test[_i][_qp] * (B[_enrichment_component] - BI[_i][_enrichment_component]) +
      _test[_i][_qp] * grad_B;
  RealVectorValue grad_phi =
      _grad_phi[_j][_qp] * (B[_enrichment_component] - BI[_j][_enrichment_component]) +
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
      coupled_component = i % 2;
      coupled_enrichment_component = i / 2;
      active = true;
    }
  }

  if (active)
  {
    // calculate the near-tip enrichement function
    std::vector<Real> B;
    std::vector<RealVectorValue> dBX, dBx;
    std::vector<std::vector<Real>> BI;

    B.resize(4);
    dBX.resize(4);
    dBx.resize(4);
    BI.resize(4);

    Real r, theta;

    for (unsigned int i = 0; i < BI.size(); ++i)
    {
      BI[i].resize(4);

      _crack_front_definition->calculateRThetaToCrackFront(
          *(_current_elem->get_node(i)), 0, r, theta);

      Real st = std::sin(theta);
      Real st2 = std::sin(theta / 2.0);
      Real ct2 = std::cos(theta / 2.0);
      Real sr = std::sqrt(r);

      BI[i][0] = sr * st2;
      BI[i][1] = sr * ct2;
      BI[i][2] = sr * st2 * st;
      BI[i][3] = sr * ct2 * st;
    }

    _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], 0, r, theta);

    Real st = std::sin(theta);
    Real ct = std::cos(theta);
    Real st2 = std::sin(theta / 2.0);
    Real ct2 = std::cos(theta / 2.0);
    Real st15 = std::sin(1.5 * theta);
    Real ct15 = std::cos(1.5 * theta);
    Real sr = std::sqrt(r);

    B[0] = sr * st2;
    B[1] = sr * ct2;
    B[2] = sr * st2 * st;
    B[3] = sr * ct2 * st;

    dBx[0](0) = -0.5 / sr * st2;
    dBx[0](1) = 0.5 / sr * ct2;
    dBx[0](2) = 0.0;
    dBx[1](0) = 0.5 / sr * ct2;
    dBx[1](1) = 0.5 / sr * st2;
    dBx[1](2) = 0.0;
    dBx[2](0) = -0.5 / sr * st15 * st;
    dBx[2](1) = 0.5 / sr * (st2 + st15 * ct);
    dBx[2](2) = 0.0;
    dBx[3](0) = -0.5 / sr * ct15 * st;
    dBx[3](1) = 0.5 / sr * (ct2 + ct15 * ct);
    dBx[3](2) = 0.0;

    for (unsigned int i = 0; i < 4; ++i)
      dBX[i] = _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dBx[i],
                                                                           0); // TODO: point index

    RealVectorValue grad_B_test(dBX[_enrichment_component]);
    RealVectorValue grad_B_phi(dBX[coupled_enrichment_component]);

    RealVectorValue grad_test =
        _grad_test[_i][_qp] * (B[_enrichment_component] - BI[_i][_enrichment_component]) +
        _test[_i][_qp] * grad_B_test;
    RealVectorValue grad_phi = _grad_phi[_j][_qp] * (B[coupled_enrichment_component] -
                                                     BI[_j][coupled_enrichment_component]) +
                               _phi[_j][_qp] * grad_B_phi;

    return ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], _component, coupled_component, grad_test, grad_phi);
  }

  return 0;
}
