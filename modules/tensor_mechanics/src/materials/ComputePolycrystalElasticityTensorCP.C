/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputePolycrystalElasticityTensorCP.h"
#include "RotationTensor.h"
#include "GrainTracker.h"
#include "EulerAngleProvider.h"
#include "Conversion.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<ComputePolycrystalElasticityTensorCP>()
{
  InputParameters params = validParams<ComputeElasticityTensorBase>();
  params.addClassDescription("compute an evolving elasticity tensor coupled to a grain growth phase field model.");
  params.addRequiredParam<UserObjectName>("euler_angle_provider", "Name of Euler angle provider user object");
  params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
  params.addRequiredParam<UserObjectName>("graintracker_object", "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<unsigned int>("grain_num", "Number of initial grains that will be modeled");
  params.addRequiredParam<std::vector<Real> >("C_ijkl", "Vector containing elastic constants for fill method");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputePolycrystalElasticityTensorCP::ComputePolycrystalElasticityTensorCP(const InputParameters & parameters) :
    ComputeElasticityTensorBase(parameters),
    _C_unrotated(getParam<std::vector<Real> >("C_ijkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _grain_tracker(getUserObject<GrainTracker>("graintracker_object")),
    _grain_num(getParam<unsigned int>("grain_num")),
    _nop(coupledComponents("v")),
    _vals(_nop),
    _elastic_tensor_cp(declareProperty<std::vector<RankFourTensor> >("elastic_tensor_cp")),
    _C_rotated(_grain_num),
    _rot_trans(_grain_num),
    _crysrot(declareProperty<std::vector<RankTwoTensor> >("crysrot"))
{
  // Read in Euler angles from the Euler angle provider
  if (_euler.getGrainNum() < _grain_num)
    mooseError("Euler angle provider has too few angles.");

  // Loop over grains
  for (unsigned int grn = 0; grn < _grain_num; ++grn)
  {
    // Rotate one elasticity tensor for each grain
    RotationTensor R(_euler.getEulerAngles(grn));
    _rot_trans[grn] = R.transpose();
    _C_rotated[grn] = _C_unrotated;
    _C_rotated[grn].rotate(_rot_trans[grn]);
  }

  // Loop over variables (ops)
  for (unsigned int op = 0; op < _nop; ++op)
    // Initialize variables
    _vals[op] = &coupledValue("v", op);
}

void
ComputePolycrystalElasticityTensorCP::computeQpElasticityTensor()
{
  // Initialize local elasticity tnesor and sum of h
  RankFourTensor local_elasticity_tensor;

  Real sum_h = 0.0;

  // Get list of active order parameters from grain tracker
  const std::vector<std::pair<unsigned int, unsigned int> > & active_ops = _grain_tracker.getElementalValues(_current_elem->id());

  unsigned int n_active_ops= active_ops.size();

  _crysrot[_qp].resize(n_active_ops);

  _elastic_tensor_cp[_qp].resize(n_active_ops);

  if (n_active_ops < 1 && _t_step > 0)
    mooseError("No active order parameters");

  // Calculate elasticity tensor
  for (unsigned int op = 0; op<n_active_ops; ++op)
  {
    unsigned int grn_index = active_ops[op].first;

    // Second position contains the order parameter index
    unsigned int op_index = active_ops[op].second;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)))/2.0;

    // Sum all rotated elasticity tensors
    local_elasticity_tensor += _C_rotated[op_index] * h;
    sum_h += h;
  }

  Real tol = 1.0e-10;
  if (sum_h < tol)
    sum_h = tol;

  local_elasticity_tensor /= sum_h;

  // Fill in the matrix stiffness material property
  _elasticity_tensor[_qp] = local_elasticity_tensor;

  for (unsigned int op = 0; op < n_active_ops; ++op)
  {
    unsigned int grn_index = active_ops[op].first;
    //unsigned int op_index = active_ops[op].second;

    _crysrot[_qp][op] = _rot_trans[grn_index]; //TODO

    // Fill in material property
    _elastic_tensor_cp[_qp][op] = _C_rotated[grn_index];
  }
}
