/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeCrackTipEnrichmentSmallStrain.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/string_to_enum.h"

template <>
InputParameters
validParams<ComputeCrackTipEnrichmentSmallStrain>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<std::vector<NonlinearVariableName>>("enrichment_displacement",
                                                              "The enrichment displacement");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define multiple "
                               "mechanics material systems on the same block, i.e. for multiple "
                               "phases");
  params.addRequiredParam<UserObjectName>("crack_front_definition",
                                          "The CrackFrontDefinition user object name");
  return params;
}

ComputeCrackTipEnrichmentSmallStrain::ComputeCrackTipEnrichmentSmallStrain(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _enrich_disp(3),
    _grad_enrich_disp(3),
    _enrich_variable(4),
    _enrich_strain(declareProperty<RankTwoTensor>("enrich_strain")),
    _phi(_assembly.phi()),
    _grad_phi(_assembly.gradPhi()),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mechanical_strain(declareProperty<RankTwoTensor>(_base_name + "mechanical_strain")),
    _total_strain(declareProperty<RankTwoTensor>(_base_name + "total_strain")),
    _eigenstrain(getDefaultMaterialProperty<RankTwoTensor>(_base_name + "stress_free_strain")),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _B(4),
    _dBX(4),
    _dBx(4),
    _BI(4)
{
  for (unsigned int i = 0; i < _enrich_variable.size(); ++i)
    _enrich_variable[i].resize(3); // TODO 3D

  const std::vector<NonlinearVariableName> & nl_vnames =
      getParam<std::vector<NonlinearVariableName>>("enrichment_displacement");
  NonlinearSystem & nl = _fe_problem.getNonlinearSystem();

  unsigned int num_enrichment_funcs = 4; // TODO!!!
  for (unsigned int i = 0; i < num_enrichment_funcs;
       ++i) // TODO : total 4 enrichment functions per node along one direction
  {
    _enrich_variable[i][0] = &(nl.getVariable(0, nl_vnames[i * _ndisp]));
    _enrich_variable[i][1] = &(nl.getVariable(0, nl_vnames[i * _ndisp + 1]));
  }

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError(
        "The number of variables supplied in 'displacements' must match the mesh dimension.");

  // fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValue("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }

  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
  }

  for (unsigned int i = 0; i < _BI.size(); ++i)
    _BI[i].resize(4);

  _nl = &(_fe_problem.getNonlinearSystem());
}

void
ComputeCrackTipEnrichmentSmallStrain::computeQpProperties()
{
  _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], 0, _r, _theta);

  if (MooseUtils::absoluteFuzzyEqual(_r, 0.0))
    mooseError("ComputeCrackTipEnrichmentSmallStrain: the distance between a point and the crack "
               "tip/front is zero.");

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
                                                                          0); // TODO: point index

  _sln = _nl->currentSolution();

  for (unsigned int m = 0; m < _ndisp; ++m) // TODO: 3D
  {
    _enrich_disp[m] = 0.0;
    _grad_enrich_disp[m].zero();
    for (unsigned int i = 0; i < _current_elem->n_nodes(); ++i)
    {
      Node * node_i = _current_elem->get_node(i);
      for (unsigned int j = 0; j < 4; ++j)
      {
        dof_id_type dof = node_i->dof_number(_nl->number(), _enrich_variable[j][m]->number(), 0);
        Real soln = (*_sln)(dof);
        _enrich_disp[m] += (*_fe_phi)[i][_qp] * (_B[j] - _BI[i][j]) * soln;
        RealVectorValue grad_B(_dBX[j]);
        _grad_enrich_disp[m] +=
            ((*_fe_dphi)[i][_qp] * (_B[j] - _BI[i][j]) + (*_fe_phi)[i][_qp] * grad_B) * soln;
      }
    }
  }

  RankTwoTensor grad_tensor_enrich(
      _grad_enrich_disp[0], _grad_enrich_disp[1], _grad_enrich_disp[2]);

  _enrich_strain[_qp] = (grad_tensor_enrich + grad_tensor_enrich.transpose()) / 2.0;

  RankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);

  _total_strain[_qp] = (grad_tensor + grad_tensor.transpose()) / 2.0;

  _total_strain[_qp] += _enrich_strain[_qp];

  _mechanical_strain[_qp] = _total_strain[_qp];

  // Remove the Eigen strain
  _mechanical_strain[_qp] -= _eigenstrain[_qp];
}

void
ComputeCrackTipEnrichmentSmallStrain::initQpStatefulProperties()
{
  _enrich_strain[_qp].zero();
  _mechanical_strain[_qp].zero();
  _total_strain[_qp].zero();
}

void
ComputeCrackTipEnrichmentSmallStrain::computeProperties()
{
  FEType fe_type(Utility::string_to_enum<Order>("first"),
                 Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule(_qrule);
  _fe_phi = &(fe->get_phi());
  _fe_dphi = &(fe->get_dphi());
  fe->reinit(_current_elem);

  for (unsigned int i = 0; i < _BI.size(); ++i)
  {
    _crack_front_definition->calculateRThetaToCrackFront(
        *(_current_elem->get_node(i)), 0, _r, _theta);

    Real st = std::sin(_theta);
    Real st2 = std::sin(_theta / 2.0);
    Real ct2 = std::cos(_theta / 2.0);
    Real sr = std::sqrt(_r);

    _BI[i][0] = sr * st2;
    _BI[i][1] = sr * ct2;
    _BI[i][2] = sr * st2 * st;
    _BI[i][3] = sr * ct2 * st;
  }

  Material::computeProperties();
}
