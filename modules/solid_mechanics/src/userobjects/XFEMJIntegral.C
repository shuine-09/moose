/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This userobject calculates the J-Integral
//
#include "XFEMJIntegral.h"
#include "libmesh/fe_interface.h"

template<>
InputParameters validParams<XFEMJIntegral>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredParam<UserObjectName>("crack_front_definition","The CrackFrontDefinition user object name");
  //params.addParam<unsigned int>("crack_front_point_index","The index of the point on the crack front corresponding to this q function");
  params.addParam<Real>("radius_inner", "Inner radius for volume integral domain");
  params.addParam<Real>("radius_outer", "Outer radius for volume integral domain");
  params.addParam<bool>("convert_J_to_K",false,"Convert J-integral to stress intensity factor K.");
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.addParam<Real>("poissons_ratio","Poisson's ratio");
  params.addParam<Real>("youngs_modulus","Young's modulus of the material.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

XFEMJIntegral::XFEMJIntegral(const InputParameters & parameters):
    ElementUserObject(parameters),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    //_has_crack_front_point_index(isParamValid("crack_front_point_index")),
    //_crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _treat_as_2d(false),
    _Eshelby_tensor(getMaterialProperty<ColumnMajorMatrix>("Eshelby_tensor")),
    _J_thermal_term_vec(hasMaterialProperty<RealVectorValue>("J_thermal_term_vec")?
                        &getMaterialProperty<RealVectorValue>("J_thermal_term_vec"):
                        NULL),
    _convert_J_to_K(getParam<bool>("convert_J_to_K")),
    _has_symmetry_plane(isParamValid("symmetry_plane")),
    _poissons_ratio(isParamValid("poissons_ratio") ? getParam<Real>("poissons_ratio") : 0),
    _youngs_modulus(isParamValid("youngs_modulus") ? getParam<Real>("youngs_modulus") : 0),
    _qp(0),
    _mesh(_subproblem.mesh())
{

  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMarkerUserObject");
  _xfem = fe_problem->get_xfem();

  if (isParamValid("radius_inner") && isParamValid("radius_outer"))
  {
    _radius_inner = getParam<Real>("radius_inner");
    _radius_outer = getParam<Real>("radius_outer");
  }
  else
    mooseError("DomainIntegral error: must set radius_inner and radius_outer.");
}

void
XFEMJIntegral::initialize()
{
  _J_thermal_term_vec = hasMaterialProperty<RealVectorValue>("J_thermal_term_vec") ? &getMaterialProperty<RealVectorValue>("J_thermal_term_vec") : NULL;

  _treat_as_2d = _crack_front_definition->treatAs2D();

  if (_treat_as_2d)
  {
    if (_has_crack_front_point_index)
    {
      mooseWarning("crack_front_point_index ignored because CrackFrontDefinition is set to treat as 2D");
    }
  }
  if (_convert_J_to_K && (!isParamValid("youngs_modulus") || !isParamValid("poissons_ratio")))
    mooseError("youngs_modulus and poissons_ratio must be specified if convert_J_to_K = true");


  _num_crack_front_points = _crack_front_definition->getNumCrackFrontPoints();
  _integral_values.resize(_num_crack_front_points);

  std::cout << "number of crack front points = " << _num_crack_front_points << std::endl;
}

DenseVector<Real>
XFEMJIntegral::computeIntegrals()
{
  FEType fe_type(Utility::string_to_enum<Order>("first"),Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule (_qrule);

  // The values of the shape functions at the quadrature points
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  fe->reinit (_current_elem); 
  
  DenseVector<Real> sums(_num_crack_front_points);
  DenseVector<Real> QpIntegrals(_num_crack_front_points);

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
  {
    Real JxW_and_coord = _JxW[_qp]*_coord[_qp];
    QpIntegrals = computeQpIntegrals(phi,dphi);
    QpIntegrals *= JxW_and_coord;  
    sums += QpIntegrals;
  }
  return sums;
}

Real
XFEMJIntegral::calcQValue(Node* node)
{
  Real dist_to_crack_front;
  Real dist_along_tangent;
  projectToFrontAtPoint(dist_to_crack_front,dist_along_tangent,*node);
  bool is_point_on_intersecting_boundary = _crack_front_definition->isPointWithIndexOnIntersectingBoundary(_crack_front_point_index);

  Real q = 1.0;
  if ( dist_to_crack_front > _radius_inner &&
       dist_to_crack_front < _radius_outer)
    q = (_radius_outer - dist_to_crack_front) /
      (_radius_outer - _radius_inner);
  else if ( dist_to_crack_front >= _radius_outer)
    q = 0.0;

  if (q > 0.0)
  {
    Real tangent_multiplier = 1.0;
    if (!_treat_as_2d)
    {   
      const Real forward_segment_length = _crack_front_definition->getCrackFrontForwardSegmentLength(_crack_front_point_index);
      const Real backward_segment_length = _crack_front_definition->getCrackFrontBackwardSegmentLength(_crack_front_point_index);

      if (dist_along_tangent >= 0.0)
      { 
        if (forward_segment_length > 0.0)
          tangent_multiplier = 1.0 - dist_along_tangent/forward_segment_length;
      }
      else
      {
        if (backward_segment_length > 0.0)
          tangent_multiplier = 1.0 + dist_along_tangent/backward_segment_length;
      }
    }

    tangent_multiplier=std::max(tangent_multiplier,0.0);
    tangent_multiplier=std::min(tangent_multiplier,1.0);

    //Set to zero if a node is on a designated free surface and its crack front node is not.
    if (_crack_front_definition->isNodeOnIntersectingBoundary(node) &&
        !is_point_on_intersecting_boundary)
      tangent_multiplier=0.0;

    q *= tangent_multiplier;
  }
  return q;
}

DenseVector<Real>
XFEMJIntegral::computeQpIntegrals(const std::vector<std::vector<Real> > & N_shape_func, const std::vector<std::vector<RealGradient> > & dN_shape_func)
{
  DenseVector<Real> QpIntegrals(_num_crack_front_points);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    _crack_front_point_index = i; // at i th crack front front point

    unsigned int n_nodes = _current_elem->n_nodes();
    std::vector<Real> q_nodes(n_nodes,0.0);
    RealVectorValue grad_of_scalar_q(0.0,0.0,0.0);
    Real scalar_q = 0.0;

    // calculate Q function at finite element node
    for (unsigned i = 0; i < n_nodes; i++)
    {
      Real q = calcQValue(_current_elem->get_node(i));
      q_nodes[i] = q;
    }

    // calcuate the Q function and its gradient at quadrature point
    for (unsigned i = 0; i < n_nodes; i++)
    {
      grad_of_scalar_q(0) += q_nodes[i] * dN_shape_func[i][_qp](0);
      grad_of_scalar_q(1) += q_nodes[i] * dN_shape_func[i][_qp](1);
      grad_of_scalar_q(2) += q_nodes[i] * dN_shape_func[i][_qp](2);
      scalar_q += q_nodes[i] * N_shape_func[i][_qp];
    }


    ColumnMajorMatrix grad_of_vector_q;
    const RealVectorValue& crack_direction = _crack_front_definition->getCrackDirection(_crack_front_point_index);
    grad_of_vector_q(0,0) = crack_direction(0)*grad_of_scalar_q(0);
    grad_of_vector_q(0,1) = crack_direction(0)*grad_of_scalar_q(1);
    grad_of_vector_q(0,2) = crack_direction(0)*grad_of_scalar_q(2);
    grad_of_vector_q(1,0) = crack_direction(1)*grad_of_scalar_q(0);
    grad_of_vector_q(1,1) = crack_direction(1)*grad_of_scalar_q(1);
    grad_of_vector_q(1,2) = crack_direction(1)*grad_of_scalar_q(2);
    grad_of_vector_q(2,0) = crack_direction(2)*grad_of_scalar_q(0);
    grad_of_vector_q(2,1) = crack_direction(2)*grad_of_scalar_q(1);
    grad_of_vector_q(2,2) = crack_direction(2)*grad_of_scalar_q(2);

    Real eq = _Eshelby_tensor[_qp].doubleContraction(grad_of_vector_q);

    //Thermal component
    Real eq_thermal = 0.0;
    if (_J_thermal_term_vec)
    {
      for (unsigned int i = 0; i < 3; i++)
        eq_thermal += crack_direction(i)*scalar_q*(*_J_thermal_term_vec)[_qp](i);
    }

    Real q_avg_seg = 1.0;
    if (!_crack_front_definition->treatAs2D())
    {
      q_avg_seg = (_crack_front_definition->getCrackFrontForwardSegmentLength(_crack_front_point_index) +
                   _crack_front_definition->getCrackFrontBackwardSegmentLength(_crack_front_point_index)) / 2.0;
    }

    Real etot = -eq + eq_thermal;
    etot /= q_avg_seg;
    QpIntegrals(i) = etot;
  }

  return QpIntegrals;
}

void 
XFEMJIntegral::execute()
{
  DenseVector<Real> comp_integ(_num_crack_front_points);
  comp_integ = computeIntegrals();
  _integral_values += comp_integ;
}

void
XFEMJIntegral::threadJoin(const UserObject & y)
{
  const XFEMJIntegral & pps = static_cast<const XFEMJIntegral &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
    _integral_values(i) += pps._integral_values(i);
}

void 
XFEMJIntegral::finalize()
{
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    gatherSum(_integral_values(i));
    if (_has_symmetry_plane)
      _integral_values(i) *= 2.0;

    Real sign = (_integral_values(i) > 0.0) ? 1.0 : ((_integral_values(i) < 0.0) ? -1.0: 0.0);
    if (_convert_J_to_K)
      _integral_values(i) = sign * std::sqrt(std::abs(_integral_values(i)) * _youngs_modulus / (1 - std::pow(_poissons_ratio,2)));
  }

  for (unsigned i = 0; i < _num_crack_front_points; i++)
  {
    _map_jintegral.insert(std::pair<unsigned int, Real>(i, _integral_values(i)));
    std::cout << "crack front index (" << i << ") : J integral = " << _integral_values(i) << std::endl; 
  } 
}

void
XFEMJIntegral::projectToFrontAtPoint(Real & dist_to_front, Real & dist_along_tangent, Point p)
{
  const Point *crack_front_point = _crack_front_definition->getCrackFrontPoint(_crack_front_point_index);
  
  const RealVectorValue &crack_front_tangent =
    _crack_front_definition->getCrackFrontTangent(_crack_front_point_index);

  RealVectorValue crack_node_to_current_node = p - *crack_front_point;
  dist_along_tangent = crack_node_to_current_node * crack_front_tangent;
  RealVectorValue projection_point = *crack_front_point + dist_along_tangent * crack_front_tangent;
  RealVectorValue axis_to_current_node = p - projection_point;
  dist_to_front = axis_to_current_node.size();
}
