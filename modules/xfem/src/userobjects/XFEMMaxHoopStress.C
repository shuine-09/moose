/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This userobject calculates the configurational force
//
#include "XFEMMaxHoopStress.h"
#include "libmesh/fe_interface.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"
#include "XFEM.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"
#include "libmesh/utility.h"

template <>
InputParameters
validParams<XFEMMaxHoopStress>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<Real>("radius_inner", "Inner radius for volume integral domain");
  params.addParam<Real>("radius_outer", "Outer radius for volume integral domain");
  params.addParam<Real>("thermal_expansion", 0.0, "Thermal expansion");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>("poissons_ratio", "Poisson's ratio for the material.");
  params.addParam<Real>("youngs_modulus", "Young's modulus of the material.");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addCoupledVar("temp",
                       "The temperature (optional). Must be provided to correctly compute "
                       "stress intensity factors in models with thermal strain gradients.");
  params.addParam<BoundaryName>("intersecting_boundary", "Boundary intersected by ends of crack.");
  params.addParam<PostprocessorName>("average_h", "Postprocessor that gives average element size");
  params.addParam<Real>("critical_k", 0.0, "Critical hoop stress.");
  params.addParam<bool>("use_weibull", false, "Use weibull distribution to initiate crack?");
  return params;
}

XFEMMaxHoopStress::XFEMMaxHoopStress(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _ndisp(coupledComponents("displacements")),
    // _stress(hasMaterialProperty<RankTwoTensor>("stress")
    //             ? &getMaterialPropertyByName<RankTwoTensor>("stress")
    //             : nullptr),
    // _strain(hasMaterialProperty<RankTwoTensor>("elastic_strain")
    //             ? &getMaterialPropertyByName<RankTwoTensor>("elastic_strain")
    //             : nullptr),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain")),
    _grad_disp(3),
    _has_temp(isCoupled("temp")),
    _grad_temp(_has_temp ? coupledGradient("temp") : _grad_zero),
    _qp(0),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _mesh(_subproblem.mesh()),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _youngs_modulus(getParam<Real>("youngs_modulus")),
    _postprocessor(isParamValid("average_h") ? &getPostprocessorValue("average_h") : NULL),
    _critical_k(getParam<Real>("critical_k")),
    _use_weibull(getParam<bool>("use_weibull")),
    _weibull(_use_weibull ? &getMaterialProperty<Real>("weibull") : nullptr)
{
  // if (!hasMaterialProperty<RankTwoTensor>("stress"))
  //   mooseError("InteractionIntegral Error: RankTwoTensor material property 'stress' not found. "
  //              "This may be because solid mechanics system is being used to calculate a
  //              SymmTensor "
  //              "'stress' material property. To use interaction integral calculation with solid "
  //              "mechanics application, please set 'solid_mechanics = true' in the DomainIntegral
  //              " "block.");

  // if (!hasMaterialProperty<RankTwoTensor>("elastic_strain"))
  //   mooseError("InteractionIntegral Error: RankTwoTensor material property 'elastic_strain' not "
  //              "found. This may be because solid mechanics system is being used to calculate a "
  //              "SymmTensor 'elastic_strain' material property. To use interaction integral "
  //              "calculation with solid mechanics application, please set 'solid_mechanics = true'
  //              " "in the DomainIntegral block.");

  _fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);
  if (_fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMaxHoopStress");

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem->getXFEM());

  if (isParamValid("radius_inner") && isParamValid("radius_outer"))
  {
    _radius_inner = getParam<Real>("radius_inner");
    _radius_outer = getParam<Real>("radius_outer");
  }
  else
    mooseError("XFEMMaxHoopStress error: must set radius.");

  if (isParamValid("intersecting_boundary"))
    _intersecting_boundary_name = getParam<BoundaryName>("intersecting_boundary");

  // plane strain
  _kappa = 3 - 4 * _poissons_ratio;
  _shear_modulus = _youngs_modulus / (2 * (1 + _poissons_ratio));

  _K_factor = 0.5 * _youngs_modulus / (1 - std::pow(_poissons_ratio, 2));

  _sif_mode.push_back(SIF_MODE(0));
  _sif_mode.push_back(SIF_MODE(1));

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("InteractionIntegral Error: number of variables supplied in 'displacements' must "
               "match the mesh dimension.");

  // fetch gradients of coupled variables
  for (unsigned int i = 0; i < _ndisp; ++i)
    _grad_disp[i] = &coupledGradient("displacements", i);

  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
    _grad_disp[i] = &_grad_zero;
}

void
XFEMMaxHoopStress::computeAuxFields(const SIF_MODE sif_mode,
                                    RankTwoTensor & aux_stress,
                                    RankTwoTensor & grad_disp)
{
  RealVectorValue k(0.0);
  if (sif_mode == KI)
    k(0) = 1.0;
  else if (sif_mode == KII)
    k(1) = 1.0;

  Real t = _theta;
  Real t2 = _theta / 2.0;
  Real tt2 = 3.0 * _theta / 2.0;
  Real st = std::sin(t);
  Real ct = std::cos(t);
  Real st2 = std::sin(t2);
  Real ct2 = std::cos(t2);
  Real stt2 = std::sin(tt2);
  Real ctt2 = std::cos(tt2);
  Real ct2sq = Utility::pow<2>(ct2);
  Real ct2cu = Utility::pow<3>(ct2);
  Real sqrt2PiR = std::sqrt(2.0 * libMesh::pi * _r);

  // Calculate auxiliary stress tensor
  aux_stress.zero();

  aux_stress(0, 0) =
      1.0 / sqrt2PiR * (k(0) * ct2 * (1.0 - st2 * stt2) - k(1) * st2 * (2.0 + ct2 * ctt2));
  aux_stress(1, 1) = 1.0 / sqrt2PiR * (k(0) * ct2 * (1.0 + st2 * stt2) + k(1) * st2 * ct2 * ctt2);
  aux_stress(0, 1) = 1.0 / sqrt2PiR * (k(0) * ct2 * st2 * ctt2 + k(1) * ct2 * (1.0 - st2 * stt2));
  aux_stress(0, 2) = -1.0 / sqrt2PiR * k(2) * st2;
  aux_stress(1, 2) = 1.0 / sqrt2PiR * k(2) * ct2;
  // plane stress
  // Real s33 = 0;
  // plane strain
  aux_stress(2, 2) = _poissons_ratio * (aux_stress(0, 0) + aux_stress(1, 1));

  aux_stress(1, 0) = aux_stress(0, 1);
  aux_stress(2, 0) = aux_stress(0, 2);
  aux_stress(2, 1) = aux_stress(1, 2);

  // Calculate x1 derivative of auxiliary displacements
  grad_disp.zero();

  grad_disp(0, 0) = k(0) / (4.0 * _shear_modulus * sqrt2PiR) *
                        (ct * ct2 * _kappa + ct * ct2 - 2.0 * ct * ct2cu + st * st2 * _kappa +
                         st * st2 - 6.0 * st * st2 * ct2sq) +
                    k(1) / (4.0 * _shear_modulus * sqrt2PiR) *
                        (ct * st2 * _kappa + ct * st2 + 2.0 * ct * st2 * ct2sq - st * ct2 * _kappa +
                         3.0 * st * ct2 - 6.0 * st * ct2cu);

  grad_disp(0, 1) = k(0) / (4.0 * _shear_modulus * sqrt2PiR) *
                        (ct * st2 * _kappa + ct * st2 - 2.0 * ct * st2 * ct2sq - st * ct2 * _kappa -
                         5.0 * st * ct2 + 6.0 * st * ct2cu) +
                    k(1) / (4.0 * _shear_modulus * sqrt2PiR) *
                        (-ct * ct2 * _kappa + 3.0 * ct * ct2 - 2.0 * ct * ct2cu -
                         st * st2 * _kappa + 3.0 * st * st2 - 6.0 * st * st2 * ct2sq);

  grad_disp(0, 2) = k(2) / (_shear_modulus * sqrt2PiR) * (st2 * ct - ct2 * st);
}

void
XFEMMaxHoopStress::initialize()
{
  _crack_front_points.clear();
  _crack_front_directions.clear();
  _elem_id_crack_tip.clear();
  _integral_values.clear();

  _weibull_at_tip.clear();

  _num_crack_front_points = _xfem->numberCrackTips();

  _integral_values.resize(_num_crack_front_points * 2);

  _weibull_at_tip.resize(_num_crack_front_points);

  for (unsigned int i = 0; i < _num_crack_front_points * 2; i++)
    _integral_values[i] = 0.0;

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
    _weibull_at_tip[i] = 0.0;

  _xfem->getCrackTipOriginDirection(
      _elem_id_crack_tip, _crack_front_points, _crack_front_directions);

  for (unsigned int i = 0; i < _crack_front_points.size(); i++)
  {
    std::cout << "MAX Hoop Stress: crack_front_points " << _crack_front_points[i] << std::endl;
    std::cout << "MAX Hoop Stress: crack_front_directions " << _crack_front_directions[i]
              << std::endl;
  }
}

std::vector<Real>
XFEMMaxHoopStress::computeIntegrals()
{
  FEType fe_type(Utility::string_to_enum<Order>("first"),
                 Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule(_qrule);

  // The values of the shape functions at the quadrature points
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  fe->reinit(_current_elem);

  std::vector<Real> sums(_num_crack_front_points * 2);
  for (unsigned int i = 0; i < _num_crack_front_points * 2; i++)
    sums[i] = 0.0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    std::vector<Real> QpIntegrals = computeQpIntegrals(phi, dphi);
    for (unsigned int i = 0; i < _num_crack_front_points * 2; i++)
    {
      sums[i] += QpIntegrals[i] * _JxW[_qp] * _coord[_qp];
    }
  }
  return sums;
}

Real
XFEMMaxHoopStress::calcQValue(Point & node, Point & crack_front)
{
  Point dist_to_crack_front_vector = node - crack_front;
  Real dist_to_crack_front = std::pow(dist_to_crack_front_vector.size_sq(), 0.5);

  Real q = 1.0;
  if (dist_to_crack_front > _radius_inner && dist_to_crack_front < _radius_outer)
    q = (_radius_outer - dist_to_crack_front) / (_radius_outer - _radius_inner);
  else if (dist_to_crack_front >= _radius_outer)
    q = 0.0;

  return q;
}

RealVectorValue
XFEMMaxHoopStress::rotateToCrackFrontCoords(const RealVectorValue vector) const
{
  return _rot_matrix * vector;
}

RealVectorValue
XFEMMaxHoopStress::rotateToCrackFrontCoords(const Point point) const
{
  RealVectorValue vector(point(0), point(1), point(2));
  return _rot_matrix * vector;
}

RankTwoTensor
XFEMMaxHoopStress::rotateToCrackFrontCoords(const RankTwoTensor tensor) const
{
  RankTwoTensor tmp_tensor(tensor);
  tmp_tensor.rotate(_rot_matrix);
  return tmp_tensor;
}

void
XFEMMaxHoopStress::calcRTheta(Point & p, Point & crack_front_point, Point & crack_direction)
{
  RealVectorValue tangent_direction;
  tangent_direction(2) = 1.0;
  _crack_plane_normal = tangent_direction.cross(crack_direction);
  _rot_matrix.fillRow(0, crack_direction);
  _rot_matrix.fillRow(1, _crack_plane_normal);
  _rot_matrix(2, 2) = 1.0;

  Point closest_point(0.0);
  RealVectorValue crack_front_point_rot = rotateToCrackFrontCoords(crack_front_point);

  RealVectorValue crack_front_edge = rotateToCrackFrontCoords(tangent_direction);

  Point p_rot = rotateToCrackFrontCoords(p);
  p_rot = p_rot - crack_front_point_rot;

  RealVectorValue closest_point_to_p = p_rot;

  // Find r, the distance between the qp and the crack front
  RealVectorValue r_vec = p_rot;
  _r = r_vec.size();

  // Find theta, the angle between r and the crack front plane
  RealVectorValue crack_plane_normal = rotateToCrackFrontCoords(_crack_plane_normal);
  Real p_to_plane_dist = std::abs(closest_point_to_p * crack_plane_normal);

  // Determine if p is above or below the crack plane
  Real y_local = p_rot(1) - closest_point(1);

  // Determine if p is in front of or behind the crack front
  RealVectorValue p2(p_rot);
  p2(1) = 0;
  RealVectorValue p2_vec = p2 - closest_point;
  Real ahead = crack_front_edge(2) * p2_vec(0) - crack_front_edge(0) * p2_vec(2);

  Real x_local(0);
  if (ahead >= 0)
    x_local = 1;
  else
    x_local = -1;

  // Calculate theta based on in which quadrant in the crack front coordinate
  // system the qp is located
  if (_r > 0)
  {
    Real theta_quadrant1(0.0);
    if (MooseUtils::absoluteFuzzyEqual(_r, p_to_plane_dist, 1e-10))
      theta_quadrant1 = 0.5 * libMesh::pi;
    else if (p_to_plane_dist > _r)
      mooseError(
          "Invalid distance p_to_plane_dist in CrackFrontDefinition::calculateRThetaToCrackFront");
    else
      theta_quadrant1 = std::asin(p_to_plane_dist / _r);

    if (x_local >= 0 && y_local >= 0)
      _theta = theta_quadrant1;

    else if (x_local < 0 && y_local >= 0)
      _theta = libMesh::pi - theta_quadrant1;

    else if (x_local < 0 && y_local < 0)
      _theta = -(libMesh::pi - theta_quadrant1);

    else if (x_local >= 0 && y_local < 0)
      _theta = -theta_quadrant1;
  }
  else if (_r == 0)
    _theta = 0;
  else
    mooseError("Invalid distance r in XFEMMaxHoopStress::calculateRTheta");
}

std::vector<Real>
XFEMMaxHoopStress::computeQpIntegrals(const std::vector<std::vector<Real>> & N_shape_func,
                                      const std::vector<std::vector<RealGradient>> & dN_shape_func)
{
  std::vector<Real> QpIntegrals(_num_crack_front_points * 2);
  for (unsigned int i = 0; i < _num_crack_front_points * 2; i++)
    QpIntegrals[i] = 0.0;

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    unsigned int n_nodes = _current_elem->n_nodes();
    std::vector<Real> q_nodes(n_nodes, 0.0);
    RealVectorValue grad_of_scalar_q(0.0, 0.0, 0.0);
    Real scalar_q = 0.0;

    Point crack_front_point = _crack_front_points[i];
    Point crack_direction = _crack_front_directions[i];

    Point p = _q_point[_qp];

    calcRTheta(p, crack_front_point, crack_direction);

    // calculate Q function at finite element node
    for (unsigned int j = 0; j < n_nodes; j++)
    {
      Real q = calcQValue((*_current_elem->get_node(j)), crack_front_point);
      dof_id_type node_id = _current_elem->get_node(j)->id();
      BoundaryID intersecting_boundary_id = _mesh.getBoundaryID(_intersecting_boundary_name);
      if (_mesh.isBoundaryNode(node_id, intersecting_boundary_id))
        q = 0.0;
      q_nodes[j] = q;
    }

    // calcuate the Q function and its gradient at quadrature point
    for (unsigned int j = 0; j < n_nodes; j++)
    {
      grad_of_scalar_q(0) += q_nodes[j] * dN_shape_func[j][_qp](0);
      grad_of_scalar_q(1) += q_nodes[j] * dN_shape_func[j][_qp](1);
      grad_of_scalar_q(2) += q_nodes[j] * dN_shape_func[j][_qp](2);
      scalar_q += q_nodes[j] * N_shape_func[j][_qp];
    }

    RealVectorValue grad_q = grad_of_scalar_q;

    // Calculate auxiliary stress tensor at current qp
    for (std::vector<SIF_MODE>::iterator it = _sif_mode.begin(); it != _sif_mode.end(); ++it)
    {
      RankTwoTensor aux_stress;
      RankTwoTensor aux_du;

      if (*it == KI)
        computeAuxFields(*it, aux_stress, aux_du);

      else if (*it == KII)
        computeAuxFields(*it, aux_stress, aux_du);

      // In the crack front coordinate system, the crack direction is (1,0,0)
      RealVectorValue crack_direction_local(0.0);
      crack_direction_local(0) = 1.0;

      RankTwoTensor grad_disp((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);

      // Rotate stress, strain, and displacement to crack front coordinate system
      RealVectorValue grad_q_cf = rotateToCrackFrontCoords(grad_q);
      RankTwoTensor grad_disp_cf = rotateToCrackFrontCoords(grad_disp);
      RankTwoTensor stress_cf = rotateToCrackFrontCoords((_stress)[_qp]);
      RankTwoTensor strain_cf = rotateToCrackFrontCoords((_strain)[_qp]);

      RankTwoTensor dq;
      dq(0, 0) = crack_direction_local(0) * grad_q_cf(0);
      dq(0, 1) = crack_direction_local(0) * grad_q_cf(1);
      dq(0, 2) = crack_direction_local(0) * grad_q_cf(2);

      // Calculate interaction integral terms

      // Term1 = stress * x1-derivative of aux disp * dq
      RankTwoTensor tmp1 = dq * stress_cf;
      Real term1 = aux_du.doubleContraction(tmp1);

      // Term2 = aux stress * x1-derivative of disp * dq
      RankTwoTensor tmp2 = dq * aux_stress;
      Real term2 = grad_disp_cf(0, 0) * tmp2(0, 0) + grad_disp_cf(1, 0) * tmp2(0, 1) +
                   grad_disp_cf(2, 0) * tmp2(0, 2);

      // Term3 = aux stress * strain * dq_x   (= stress * aux strain * dq_x)
      Real term3 = dq(0, 0) * aux_stress.doubleContraction(strain_cf);

      /*
      RealVectorValue J_vec(0);

      for (unsigned int j=0; j<3; ++j)
      {
        Real dthermstrain_dx = _temp_grad[_qp](j) * _thermal_expansion;
        J_vec(j) = _aux_stress.tr()*dthermstrain_dx;
      }

      Real eq_thermal = 0.0;

      for (unsigned int j = 0; j < 3; j++)
        eq_thermal += crack_direction_local(j)*scalar_q*J_vec(j);

      */

      RealVectorValue grad_temp_cf = rotateToCrackFrontCoords(_grad_temp[_qp]);

      Real eq_thermal = 0.0;
      Real aux_stress_trace = aux_stress(0, 0) + aux_stress(1, 1) + aux_stress(2, 2);
      eq_thermal = scalar_q * aux_stress_trace * _thermal_expansion * grad_temp_cf(0);

      Real eq = term1 + term2 - term3 + eq_thermal;

      QpIntegrals[i * 2 + int(*it)] = eq;
    }
  }
  return QpIntegrals;
}

void
XFEMMaxHoopStress::execute()
{
  // if (_postprocessor)
  // {
  //   _radius_inner = 1.5 * *_postprocessor;
  //   _radius_outer = 3.5 * *_postprocessor;
  // }

  std::vector<Real> comp_integ = computeIntegrals();
  for (unsigned int i = 0; i < _num_crack_front_points * 2; i++)
    _integral_values[i] += comp_integ[i];

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    if (_current_elem == _elem_id_crack_tip[i])
    {
      if (_use_weibull)
        _weibull_at_tip[i] = (*_weibull)[0];
      else
        _weibull_at_tip[i] = 1.0;
      break;
    }
  }
}

void
XFEMMaxHoopStress::threadJoin(const UserObject & y)
{
  const XFEMMaxHoopStress & pps = static_cast<const XFEMMaxHoopStress &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points * 2; i++)
    _integral_values[i] += pps._integral_values[i];
}

void
XFEMMaxHoopStress::finalize()
{
  _xfem->clearDoesCrackGrowth();
  gatherSum(_weibull_at_tip);
  gatherSum(_integral_values);

  for (unsigned int i = 0; i < _num_crack_front_points * 2; i++)
    _integral_values[i] *= _K_factor;

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Real KI = _integral_values[i * 2];
    Real KII = _integral_values[i * 2 + 1];

    Real effective_K = std::sqrt(KI * KI + KII * KII);

    bool does_elem_crack = false;
    if (_use_weibull)
      does_elem_crack = (effective_K > _weibull_at_tip[i] * _critical_k) && (KI > 0.0);
    else
      does_elem_crack = effective_K > _critical_k;

    std::cout << "effective_K = " << effective_K << ", critical_k = " << _critical_k
              << ", weibull_tip = " << _weibull_at_tip[i]
              << " does elem crack = " << does_elem_crack << std::endl;
    _xfem->updateDoesCrackGrowth(_elem_id_crack_tip[i], does_elem_crack);
  }

  _xfem->clearCrackGrowthDirection();

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Real KI = _integral_values[i * 2];
    Real KII = _integral_values[i * 2 + 1];
    std::cout << "KI = " << KI << std::endl;
    std::cout << "KII = " << KII << std::endl;
    Real theta1 = 2 * std::atan(0.25 * (KI / KII + std::sqrt(pow(KI / KII, 2.0) + 8.0)));
    Real theta2 = 2 * std::atan(0.25 * (KI / KII - std::sqrt(pow(KI / KII, 2.0) + 8.0)));

    Real hoop_stress1 = KI * (3 * std::cos(theta1 * 0.5) + std::cos(theta1 * 1.5)) +
                        KII * (-3.0 * std::sin(theta1 * 0.5) - 3.0 * std::sin(1.5 * theta1));
    Real hoop_stress2 = KI * (3 * std::cos(theta2 * 0.5) + std::cos(theta2 * 1.5)) +
                        KII * (-3.0 * std::sin(theta2 * 0.5) - 3.0 * std::sin(1.5 * theta2));

    Real theta = 0.0;

    std::cout << "hoop_stress1 = " << hoop_stress1 << ", theta1 = " << theta1 / libMesh::pi * 180.0
              << std::endl;
    std::cout << "hoop_stress2 = " << hoop_stress2 << ", theta2 = " << theta2 / libMesh::pi * 180.0
              << std::endl;

    if (hoop_stress1 > hoop_stress2)
      theta = theta1;
    else
      theta = theta2;

    std::cout << "theta = " << theta / libMesh::pi * 180.0 << std::endl;

    Point crack_front_point = _crack_front_points[i];
    Point crack_direction = _crack_front_directions[i];

    std::cout << "crack_direction = " << crack_direction << std::endl;

    Real omega = std::atan2(crack_direction(1), crack_direction(0));

    std::cout << "omega = " << omega / libMesh::pi * 180 << std::endl;

    Point direction(std::cos(omega + theta), std::sin(omega + theta), 0.0);

    _xfem->updateCrackGrowthDirection(_elem_id_crack_tip[i], direction);
    std::cout << "MAXHOOPSTRESS crack front point = " << crack_front_point << std::endl;
    std::cout << "MAXHOOPSTRESS crack front index (" << i << ") : direction  = " << direction
              << std::endl;
  }
}
