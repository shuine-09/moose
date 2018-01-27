/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeLinearElasticStressInterface.h"
#include "MooseMesh.h"
#include "Assembly.h"

template <>
InputParameters
validParams<ComputeLinearElasticStressInterface>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", "List of eigenstrains to be applied in this strain calculation");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<Real>("bulk_modulus", "The bulk modulus for the material.");
  params.addParam<Real>("lambda", "Lame's first constant for the material.");
  params.addParam<Real>("poissons_ratio", "Poisson's ratio for the material.");
  params.addParam<Real>("shear_modulus", "The shear modulus of the material.");
  params.addParam<Real>("youngs_modulus", "Young's modulus of the material.");
  return params;
}

ComputeLinearElasticStressInterface::ComputeLinearElasticStressInterface(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mechanical_strain(declareProperty<RankTwoTensor>(_base_name + "mechanical_strain")),
    _total_strain(declareProperty<RankTwoTensor>(_base_name + "total_strain")),
    _eigenstrain_names(getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
    _eigenstrains(_eigenstrain_names.size()),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _youngs_modulus(getParam<Real>("youngs_modulus")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress"))
{
  for (unsigned int i = 0; i < _eigenstrains.size(); ++i)
  {
    _eigenstrain_names[i] = _base_name + _eigenstrain_names[i];
    _eigenstrains[i] = &getMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
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

  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The strain calculator needs to run on the undisplaced mesh.");

  _Cijkl.fillSymmetricIsotropicEandNu(_youngs_modulus, _poissons_ratio);
}

void
ComputeLinearElasticStressInterface::initQpStatefulProperties()
{
  _mechanical_strain[_qp].zero();
  _total_strain[_qp].zero();
  _stress[_qp].zero();
}

void
ComputeLinearElasticStressInterface::resetQpProperties()
{
  _mechanical_strain[_qp].zero();
  _total_strain[_qp].zero();
  _stress[_qp].zero();
}

void
ComputeLinearElasticStressInterface::computeQpProperties()
{
  RankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);

  _total_strain[_qp] = (grad_tensor + grad_tensor.transpose()) / 2.0;

  _mechanical_strain[_qp] = _total_strain[_qp];

  // Remove the Eigen strain
  for (auto es : _eigenstrains)
    _mechanical_strain[_qp] -= (*es)[_qp];
  // stress = C * e
  _stress[_qp] = _Cijkl * _mechanical_strain[_qp];
}
