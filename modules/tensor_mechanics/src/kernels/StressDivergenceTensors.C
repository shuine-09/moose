/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceTensors.h"

// MOOSE includes
#include "ElasticityTensorTools.h"
#include "Material.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"

// libMesh includes
#include "libmesh/quadrature.h"
#include "XFEM.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"


template <>
InputParameters
validParams<StressDivergenceTensors>()
{
  InputParameters params = validParams<ALEKernel>();
  params.addClassDescription("Stress divergence kernel for the Cartesian coordinate system");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addCoupledVar("temperature", "The temperature");
  params.addParam<std::string>("base_name", "Material property base name");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<bool>(
      "use_finite_deform_jacobian", false, "Jacobian for corotational finite strain");
  params.addParam<bool>("volumetric_locking_correction",
                        false,
                        "Set to false to turn off volumetric locking correction");
  return params;
}

StressDivergenceTensors::StressDivergenceTensors(const InputParameters & parameters)
  : ALEKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _use_finite_deform_jacobian(getParam<bool>("use_finite_deform_jacobian")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _temp_coupled(isCoupled("temperature")),
    _temp_var(_temp_coupled ? coupled("temperature") : 0),
    _avg_grad_test(_test.size(), std::vector<Real>(3, 0.0)),
    _avg_grad_phi(_phi.size(), std::vector<Real>(3, 0.0)),
    _volumetric_locking_correction(getParam<bool>("volumetric_locking_correction"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension.");

  if (_use_finite_deform_jacobian)
  {
    _deformation_gradient =
        &getMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient");
    _deformation_gradient_old =
        &getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient");
    _rotation_increment = &getMaterialProperty<RankTwoTensor>(_base_name + "rotation_increment");
  }

  // Error if volumetic locking correction is turned on for 1D problems
  if (_ndisp == 1 && _volumetric_locking_correction)
    mooseError("Volumetric locking correction should be set to false for 1-D problems.");
}

void
StressDivergenceTensors::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_XYZ)
    mooseError(
        "The coordinate system in the Problem block must be set to XYZ for cartesian geometries.");
  
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in EnrichStressDivergenceTensors");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in EnrichStressDivergenceTensors");

}

void
StressDivergenceTensors::computeResidual()
{
<<<<<<< HEAD
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  if (_volumetric_locking_correction)
    computeAverageGradientTest();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); ++_i)
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _save_in)
      var->sys().solution().add_vector(_local_re, var->dofIndices());
  }
=======
  Point tip_edge((_current_elem->point(0))(0), 1.0, 0);
  Point tip(0.5, 1.0, 0);
  Real elem_h = std::sqrt(_current_elem->volume());
  Point tip_split(0.5 - elem_h, 1.0, 0);
  if(_current_elem->contains_point(tip))
  {
    DenseVector<Number> & re = _assembly.residualBlock(_var.number());
    _local_re.resize(re.size());
    _local_re.zero();

    std::vector<Point> intersectionPoints;

    intersectionPoints.push_back(_current_elem->point(0));
    intersectionPoints.push_back(_current_elem->point(1));
    intersectionPoints.push_back(_current_elem->point(2));
    intersectionPoints.push_back(_current_elem->point(3));
    intersectionPoints.push_back(tip_edge);
    
    std::vector<Point> q_points;
    std::vector<Real> weights;

    _xfem->getXFEMqRuleOnSurface(intersectionPoints, tip, q_points, weights);

    FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

    fe_problem->reinitElemPhys(_current_elem, q_points, 0);

    _subproblem.prepareShapes(_var.number(), 0);

    fe_problem->reinitMaterials(_current_elem->subdomain_id(), 0, false);

    for (_i = 0; _i < _test.size(); _i++)
      for (_qp = 0; _qp < q_points.size(); _qp++)
        _local_re(_i) += weights[_qp] * computeQpResidual();

    re += _local_re;  
  }
  /*
  else if (_current_elem->contains_point(tip_split))
  {
    DenseVector<Number> & re = _assembly.residualBlock(_var.number());
    _local_re.resize(re.size());
    _local_re.zero();

    Point top(0.5- elem_h, 1.01, 0.0);

    FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

    const Elem * undisplaced_elem  = NULL;
    if(fe_problem->getDisplacedProblem() != NULL)
      undisplaced_elem = fe_problem->getDisplacedProblem()->refMesh().elemPtr(_current_elem->id());
    else
      undisplaced_elem = _current_elem;

    bool flag = _xfem->flagQpointInside(undisplaced_elem, top);


    Point embed_pt1(_current_elem->point(0)(0), 1.0, 0);
    Point embed_pt2(_current_elem->point(1)(0), 1.0, 0);

    std::cout << "embed_pt1 = " << embed_pt1 << std::endl;
    std::cout << "embed_pt2 = " << embed_pt2 << std::endl;

    if (flag > 0.5) //top cut element
    {
      std::vector<Point> intersectionPoints;
      intersectionPoints.push_back(embed_pt1);
      intersectionPoints.push_back(embed_pt2);
      intersectionPoints.push_back(_current_elem->point(2));
      intersectionPoints.push_back(_current_elem->point(3));

      std::vector<Point> q_points;
      std::vector<Real> weights;

      _xfem->getXFEMqRuleOnSurface(intersectionPoints, q_points, weights);

      fe_problem->reinitElemPhys(_current_elem, q_points, 0);

      _subproblem.prepareShapes(_var.number(), 0);

      fe_problem->reinitMaterials(_current_elem->subdomain_id(), 0, false);

      for (_i = 0; _i < _test.size(); _i++)
        for (_qp = 0; _qp < q_points.size(); _qp++)
          _local_re(_i) += weights[_qp] * computeQpResidual();

      re += _local_re;
    }
    else
    {
      std::vector<Point> intersectionPoints;
      intersectionPoints.push_back(_current_elem->point(0));
      intersectionPoints.push_back(_current_elem->point(1));

      intersectionPoints.push_back(embed_pt2);
      intersectionPoints.push_back(embed_pt1);

      std::vector<Point> q_points;
      std::vector<Real> weights;

      _xfem->getXFEMqRuleOnSurface(intersectionPoints, q_points, weights);

      fe_problem->reinitElemPhys(_current_elem, q_points, 0);

      _subproblem.prepareShapes(_var.number(), 0);

      fe_problem->reinitMaterials(_current_elem->subdomain_id(), 0, false);

      for (_i = 0; _i < _test.size(); _i++)
        for (_qp = 0; _qp < q_points.size(); _qp++)
          _local_re(_i) += weights[_qp] * computeQpResidual();

      re += _local_re;

    }
  }
  */
  else
    Kernel::computeResidual();
}

Real
StressDivergenceTensors::computeQpResidual()
{
  Real residual = _stress[_qp].row(_component) * _grad_test[_i][_qp];
  // volumetric locking correction
  if (_volumetric_locking_correction)
    residual += _stress[_qp].trace() / 3.0 *
                (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component));

  return residual;
}

void
StressDivergenceTensors::computeJacobian()
{
  /*
  if (_use_finite_deform_jacobian)
  {
    _finite_deform_Jacobian_mult.resize(_qrule->n_points());

    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      computeFiniteDeformJacobian();

    ALEKernel::computeJacobian();
  }
  else
  {
    if (_volumetric_locking_correction)
    {
      computeAverageGradientTest();
      computeAverageGradientPhi();
    }
    Kernel::computeJacobian();
<<<<<<< HEAD
  }
=======
    */
  Point tip_edge((_current_elem->point(0))(0), 1.0, 0);
  Point tip(0.5, 1.0, 0);
  Real elem_h = std::sqrt(_current_elem->volume());
  Point tip_split(0.5 - elem_h, 1.0, 0);

  if(_current_elem->contains_point(tip))
  {
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
    _local_ke.resize(ke.m(), ke.n());
    _local_ke.zero();

    std::vector<Point> intersectionPoints;

    intersectionPoints.push_back(_current_elem->point(0));
    intersectionPoints.push_back(_current_elem->point(1));
    intersectionPoints.push_back(_current_elem->point(2));
    intersectionPoints.push_back(_current_elem->point(3));
    intersectionPoints.push_back(tip_edge);

    std::vector<Point> q_points;
    std::vector<Real> weights;

    _xfem->getXFEMqRuleOnSurface(intersectionPoints, tip, q_points, weights);

    FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

    fe_problem->reinitElemPhys(_current_elem, q_points, 0);

    _subproblem.prepareShapes(_var.number(), 0);

    fe_problem->reinitMaterials(_current_elem->subdomain_id(), 0, false);


    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < q_points.size(); _qp++)
          _local_ke(_i, _j) += weights[_qp] * computeQpJacobian();

    ke += _local_ke;
  }
  /*
  else if (_current_elem->contains_point(tip_split))
  {
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
    _local_ke.resize(ke.m(), ke.n());
    _local_ke.zero();

    Point top(0.5- elem_h, 1.01, 0.0);

    FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

    const Elem * undisplaced_elem  = NULL;
    if(fe_problem->getDisplacedProblem() != NULL)
      undisplaced_elem = fe_problem->getDisplacedProblem()->refMesh().elemPtr(_current_elem->id());
    else
      undisplaced_elem = _current_elem;

    bool flag = _xfem->flagQpointInside(undisplaced_elem, top);


    Point embed_pt1(_current_elem->point(0)(0), 1.0, 0);
    Point embed_pt2(_current_elem->point(1)(0), 1.0, 0);

    if (flag > 0.5) //top cut element
    {
      std::vector<Point> intersectionPoints;
      intersectionPoints.push_back(embed_pt1);
      intersectionPoints.push_back(embed_pt2);
      intersectionPoints.push_back(_current_elem->point(2));
      intersectionPoints.push_back(_current_elem->point(3));

      std::vector<Point> q_points;
      std::vector<Real> weights;

      _xfem->getXFEMqRuleOnSurface(intersectionPoints, q_points, weights);

      fe_problem->reinitElemPhys(_current_elem, q_points, 0);

      _subproblem.prepareShapes(_var.number(), 0);

      fe_problem->reinitMaterials(_current_elem->subdomain_id(), 0, false);
      
      for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _phi.size(); _j++)
          for (_qp = 0; _qp < q_points.size(); _qp++)
            _local_ke(_i, _j) += weights[_qp] * computeQpJacobian();

      ke += _local_ke;

    }
    else
    {
      std::vector<Point> intersectionPoints;
      intersectionPoints.push_back(_current_elem->point(0));
      intersectionPoints.push_back(_current_elem->point(1));

      intersectionPoints.push_back(embed_pt2);
      intersectionPoints.push_back(embed_pt1);

      std::vector<Point> q_points;
      std::vector<Real> weights;

      _xfem->getXFEMqRuleOnSurface(intersectionPoints, q_points, weights);

      fe_problem->reinitElemPhys(_current_elem, q_points, 0);

      _subproblem.prepareShapes(_var.number(), 0);

      fe_problem->reinitMaterials(_current_elem->subdomain_id(), 0, false);

      for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _phi.size(); _j++)
          for (_qp = 0; _qp < q_points.size(); _qp++)
            _local_ke(_i, _j) += weights[_qp] * computeQpJacobian();

      ke += _local_ke;

    }
  }
  */
  else
    Kernel::computeJacobian();

}

void
StressDivergenceTensors::computeOffDiagJacobian(unsigned int jvar)
{
  if (_use_finite_deform_jacobian)
  {
    _finite_deform_Jacobian_mult.resize(_qrule->n_points());

    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      computeFiniteDeformJacobian();

    ALEKernel::computeOffDiagJacobian(jvar);
  }
  else
  {
    if (_volumetric_locking_correction)
    {
      computeAverageGradientPhi();
      computeAverageGradientTest();
    }
    Kernel::computeOffDiagJacobian(jvar);
  }
}

Real
StressDivergenceTensors::computeQpJacobian()
{
  if (_use_finite_deform_jacobian)
    return ElasticityTensorTools::elasticJacobian(_finite_deform_Jacobian_mult[_qp],
                                                  _component,
                                                  _component,
                                                  _grad_test[_i][_qp],
                                                  _grad_phi_undisplaced[_j][_qp]);

  Real sum_C3x3 = _Jacobian_mult[_qp].sum3x3();
  RealGradient sum_C3x1 = _Jacobian_mult[_qp].sum3x1();

  Real jacobian = 0.0;
  // B^T_i * C * B_j
  jacobian += ElasticityTensorTools::elasticJacobian(
      _Jacobian_mult[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);

  if (_volumetric_locking_correction)
  {
    // jacobian = Bbar^T_i * C * Bbar_j where Bbar = B + Bvol
    // jacobian = B^T_i * C * B_j + Bvol^T_i * C * Bvol_j +  Bvol^T_i * C * B_j + B^T_i * C * Bvol_j

    // Bvol^T_i * C * Bvol_j
    jacobian += sum_C3x3 * (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) *
                (_avg_grad_phi[_j][_component] - _grad_phi[_j][_qp](_component)) / 9.0;

    // B^T_i * C * Bvol_j
    jacobian += sum_C3x1(_component) * _grad_test[_i][_qp](_component) *
                (_avg_grad_phi[_j][_component] - _grad_phi[_j][_qp](_component)) / 3.0;

    // Bvol^T_i * C * B_j
    RankTwoTensor phi;
    if (_component == 0)
    {
      phi(0, 0) = _grad_phi[_j][_qp](0);
      phi(0, 1) = phi(1, 0) = _grad_phi[_j][_qp](1);
      phi(0, 2) = phi(2, 0) = _grad_phi[_j][_qp](2);
    }
    else if (_component == 1)
    {
      phi(1, 1) = _grad_phi[_j][_qp](1);
      phi(0, 1) = phi(1, 0) = _grad_phi[_j][_qp](0);
      phi(1, 2) = phi(2, 1) = _grad_phi[_j][_qp](2);
    }
    else if (_component == 2)
    {
      phi(2, 2) = _grad_phi[_j][_qp](2);
      phi(0, 2) = phi(2, 0) = _grad_phi[_j][_qp](0);
      phi(1, 2) = phi(2, 1) = _grad_phi[_j][_qp](1);
    }

    jacobian += (_Jacobian_mult[_qp] * phi).trace() *
                (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) / 3.0;
  }
  return jacobian;
}

Real
StressDivergenceTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
  // off-diagonal Jacobian with respect to a coupled displacement component
  for (unsigned int coupled_component = 0; coupled_component < _ndisp; ++coupled_component)
    if (jvar == _disp_var[coupled_component])
    {
      if (_use_finite_deform_jacobian)
        return ElasticityTensorTools::elasticJacobian(_finite_deform_Jacobian_mult[_qp],
                                                      _component,
                                                      coupled_component,
                                                      _grad_test[_i][_qp],
                                                      _grad_phi_undisplaced[_j][_qp]);

      const Real sum_C3x3 = _Jacobian_mult[_qp].sum3x3();
      const RealGradient sum_C3x1 = _Jacobian_mult[_qp].sum3x1();
      Real jacobian = 0.0;

      // B^T_i * C * B_j
      jacobian += ElasticityTensorTools::elasticJacobian(_Jacobian_mult[_qp],
                                                         _component,
                                                         coupled_component,
                                                         _grad_test[_i][_qp],
                                                         _grad_phi[_j][_qp]);

      if (_volumetric_locking_correction)
      {
        // jacobian = Bbar^T_i * C * Bbar_j where Bbar = B + Bvol
        // jacobian = B^T_i * C * B_j + Bvol^T_i * C * Bvol_j +  Bvol^T_i * C * B_j + B^T_i * C *
        // Bvol_j

        // Bvol^T_i * C * Bvol_j
        jacobian += sum_C3x3 * (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) *
                    (_avg_grad_phi[_j][coupled_component] - _grad_phi[_j][_qp](coupled_component)) /
                    9.0;

        // B^T_i * C * Bvol_j
        jacobian += sum_C3x1(_component) * _grad_test[_i][_qp](_component) *
                    (_avg_grad_phi[_j][coupled_component] - _grad_phi[_j][_qp](coupled_component)) /
                    3.0;

        // Bvol^T_i * C * B_i
        RankTwoTensor phi;
        for (unsigned int i = 0; i < 3; ++i)
          phi(coupled_component, i) = _grad_phi[_j][_qp](i);

        jacobian += (_Jacobian_mult[_qp] * phi).trace() *
                    (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) / 3.0;
      }

      return jacobian;
    }

  // off-diagonal Jacobian with respect to a coupled temperature variable
  if (_temp_coupled && jvar == _temp_var)
  {
    // return _d_stress_dT[_qp].rowDot(_component, _grad_test[_i][_qp]) * _phi[_j][_qp];
    return 0.0;
  }

  return 0.0;
}

void
StressDivergenceTensors::computeFiniteDeformJacobian()
{
  const RankTwoTensor I(RankTwoTensor::initIdentity);
  const RankFourTensor II_ijkl = I.mixedProductIkJl(I);

  // Bring back to unrotated config
  const RankTwoTensor unrotated_stress =
      (*_rotation_increment)[_qp].transpose() * _stress[_qp] * (*_rotation_increment)[_qp];

  // Incremental deformation gradient Fhat
  const RankTwoTensor Fhat =
      (*_deformation_gradient)[_qp] * (*_deformation_gradient_old)[_qp].inverse();
  const RankTwoTensor Fhatinv = Fhat.inverse();

  const RankTwoTensor rot_times_stress = (*_rotation_increment)[_qp] * unrotated_stress;
  const RankFourTensor dstress_drot =
      I.mixedProductIkJl(rot_times_stress) + I.mixedProductJkIl(rot_times_stress);
  const RankFourTensor rot_rank_four =
      (*_rotation_increment)[_qp].mixedProductIkJl((*_rotation_increment)[_qp]);
  const RankFourTensor drot_dUhatinv = Fhat.mixedProductIkJl(I);

  const RankTwoTensor A = I - Fhatinv;

  // Ctilde = Chat^-1 - I
  const RankTwoTensor Ctilde = A * A.transpose() - A - A.transpose();
  const RankFourTensor dCtilde_dFhatinv =
      -I.mixedProductIkJl(A) - I.mixedProductJkIl(A) + II_ijkl + I.mixedProductJkIl(I);

  // Second order approximation of Uhat - consistent with strain increment definition
  // const RankTwoTensor Uhat = I - 0.5 * Ctilde - 3.0/8.0 * Ctilde * Ctilde;

  RankFourTensor dUhatinv_dCtilde =
      0.5 * II_ijkl - 1.0 / 8.0 * (I.mixedProductIkJl(Ctilde) + Ctilde.mixedProductIkJl(I));
  RankFourTensor drot_dFhatinv = drot_dUhatinv * dUhatinv_dCtilde * dCtilde_dFhatinv;

  drot_dFhatinv -= Fhat.mixedProductIkJl((*_rotation_increment)[_qp].transpose());
  _finite_deform_Jacobian_mult[_qp] = dstress_drot * drot_dFhatinv;

  const RankFourTensor dstrain_increment_dCtilde =
      -0.5 * II_ijkl + 0.25 * (I.mixedProductIkJl(Ctilde) + Ctilde.mixedProductIkJl(I));
  _finite_deform_Jacobian_mult[_qp] +=
      rot_rank_four * _Jacobian_mult[_qp] * dstrain_increment_dCtilde * dCtilde_dFhatinv;
  _finite_deform_Jacobian_mult[_qp] += Fhat.mixedProductJkIl(_stress[_qp]);

  const RankFourTensor dFhat_dFhatinv = -Fhat.mixedProductIkJl(Fhat.transpose());
  const RankTwoTensor dJ_dFhatinv = dFhat_dFhatinv.innerProductTranspose(Fhat.ddet());

  // Component from Jacobian derivative
  _finite_deform_Jacobian_mult[_qp] += _stress[_qp].outerProduct(dJ_dFhatinv);

  // Derivative of Fhatinv w.r.t. undisplaced coordinates
  const RankTwoTensor Finv = (*_deformation_gradient)[_qp].inverse();
  const RankFourTensor dFhatinv_dGradu = -Fhatinv.mixedProductIkJl(Finv.transpose());
  _finite_deform_Jacobian_mult[_qp] = _finite_deform_Jacobian_mult[_qp] * dFhatinv_dGradu;
}

void
StressDivergenceTensors::computeAverageGradientTest()
{
  // Calculate volume averaged value of shape function derivative
  _avg_grad_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
  {
    _avg_grad_test[_i].resize(3);
    _avg_grad_test[_i][_component] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_grad_test[_i][_component] += _grad_test[_i][_qp](_component) * _JxW[_qp] * _coord[_qp];

    _avg_grad_test[_i][_component] /= _current_elem_volume;
  }
}

void
StressDivergenceTensors::computeAverageGradientPhi()
{
  // Calculate volume average derivatives for phi
  _avg_grad_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
  {
    _avg_grad_phi[_i].resize(3);
    for (unsigned int component = 0; component < _mesh.dimension(); ++component)
    {
      _avg_grad_phi[_i][component] = 0.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _avg_grad_phi[_i][component] += _grad_phi[_i][_qp](component) * _JxW[_qp] * _coord[_qp];

      _avg_grad_phi[_i][component] /= _current_elem_volume;
    }
  }
}
