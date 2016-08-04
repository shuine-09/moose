/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This userobject calculates the configurational force
//
#include "XFEMEnergyReleaseRate.h"
#include "libmesh/fe_interface.h"

template<>
InputParameters validParams<XFEMEnergyReleaseRate>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<Real>("radius_inner", "Inner radius for volume integral domain");
  params.addParam<Real>("radius_outer", "Outer radius for volume integral domain");
  params.addParam<Real>("critial_energy_release_rate", 0.0, "Critical energy release rate.");
  params.addParam<bool>("use_weibull", false,"Use weibull distribution to initiate crack?");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

XFEMEnergyReleaseRate::XFEMEnergyReleaseRate(const InputParameters & parameters):
    ElementUserObject(parameters),
    _Eshelby_tensor(getMaterialProperty<ColumnMajorMatrix>("Eshelby_tensor")),
    _J_thermal_term_vec(hasMaterialProperty<RealVectorValue>("J_thermal_term_vec")?
                        &getMaterialProperty<RealVectorValue>("J_thermal_term_vec"):
                        NULL),
    _weibull_eta(getMaterialProperty<Real>("weibull_eta")),
    _qp(0),
    _critical_energy_release_rate(getParam<Real>("critial_energy_release_rate")),
    _mesh(_subproblem.mesh()),
     _use_weibull(getParam<bool>("use_weibull"))
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
XFEMEnergyReleaseRate::initialize()
{
  _J_thermal_term_vec = hasMaterialProperty<RealVectorValue>("J_thermal_term_vec") ? &getMaterialProperty<RealVectorValue>("J_thermal_term_vec") : NULL;

  _crack_front_points.clear();
  _crack_directions.clear();
  _elem_id_crack_tip.clear();
  _integral_values.clear();

  _weibull_at_tip.clear();

  _num_crack_front_points = _xfem->num_crack_tips();

  _integral_values.resize(_num_crack_front_points);

  _weibull_at_tip.resize(_num_crack_front_points);
  
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    _integral_values[i] = 0.0;
    _weibull_at_tip[i] = 0.0;
  }

  _xfem->get_crack_tip_origin_and_direction(_elem_id_crack_tip, _crack_front_points, _crack_directions);
}

std::vector<Real>
XFEMEnergyReleaseRate::computeIntegrals()
{
  FEType fe_type(Utility::string_to_enum<Order>("first"),Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule (_qrule);

  // The values of the shape functions at the quadrature points
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  fe->reinit (_current_elem); 
  
  std::vector<Real> sums(_num_crack_front_points);
  for (signed int i = 0; i < _num_crack_front_points; i++)
    sums[i] = 0.0;

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
  {
    std::vector<Real> QpIntegrals = computeQpIntegrals(phi,dphi);
    for (unsigned int i = 0; i < _num_crack_front_points; i++){
      sums[i] += QpIntegrals[i] * _JxW[_qp]*_coord[_qp];
    }
  }
  return sums;
}

Real
XFEMEnergyReleaseRate::calcQValue(Point & node, Point & crack_front)
{
  Point dist_to_crack_front_vector = node - crack_front;
  Real dist_to_crack_front = std::pow(dist_to_crack_front_vector.size_sq(),0.5);
  
  Real q = 1.0;
  if ( dist_to_crack_front > _radius_inner &&
       dist_to_crack_front < _radius_outer)
    q = (_radius_outer - dist_to_crack_front) /
      (_radius_outer - _radius_inner);
  else if ( dist_to_crack_front >= _radius_outer)
    q = 0.0;

  return q;
}

std::vector<Real>
XFEMEnergyReleaseRate::computeQpIntegrals(const std::vector<std::vector<Real> > & N_shape_func, const std::vector<std::vector<RealGradient> > & dN_shape_func)
{
  std::vector<Real> QpIntegrals(_num_crack_front_points);
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
    QpIntegrals[i] = 0.0;
  
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    _crack_front_point_index = i; // at i th crack front front point

    Point crack_front = _crack_front_points[i];

    unsigned int n_nodes = _current_elem->n_nodes();
    std::vector<Real> q_nodes(n_nodes,0.0);
    RealVectorValue grad_of_scalar_q(0.0,0.0,0.0);
    Real scalar_q = 0.0;

    // calculate Q function at finite element node
    for (unsigned i = 0; i < n_nodes; i++)
    {
      Real q = calcQValue(*(_current_elem->get_node(i)), crack_front);
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

    //if (grad_of_scalar_q.size_sq() < 1.0e-10)
    //  continue;

    ColumnMajorMatrix grad_of_vector_q;
    Point crack_direction = _crack_directions[i];
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

    if (_J_thermal_term_vec != NULL)
    {
      for (unsigned int i = 0; i < 3; i++)
        eq_thermal += crack_direction(i) * scalar_q * (*_J_thermal_term_vec)[_qp](i);
    }

    QpIntegrals[i] = -eq + eq_thermal;
  }

  return QpIntegrals;
}

void 
XFEMEnergyReleaseRate::execute()
{
  std::vector<Real> comp_integ = computeIntegrals();
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    _integral_values[i] += comp_integ[i];
    if (_current_elem == _elem_id_crack_tip[i])
    {
      _weibull_at_tip[i] = _weibull_eta[0];
      break;
    }
  }
}

void
XFEMEnergyReleaseRate::threadJoin(const UserObject & y)
{
  const XFEMEnergyReleaseRate & pps = static_cast<const XFEMEnergyReleaseRate &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
    _integral_values[i] += pps._integral_values[i];
}

void 
XFEMEnergyReleaseRate::finalize()
{
  _xfem->clear_doesElemCrackEnergyReleaseRate();
  gatherSum(_integral_values);
  gatherSum(_weibull_at_tip);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    //std::cout << "crack front index (" << i << ") : J_Integral  = " << _integral_values[i] << std::endl;
    //std::cout << "crack extension direction =  " << _crack_directions[i] << std::endl;
    std::cout << "energy release rate = " << _integral_values[i] << std::endl;
    bool does_elem_crack = false;
    if (_use_weibull)
      does_elem_crack = _integral_values[i] > _critical_energy_release_rate * _weibull_at_tip[i];
    else
      does_elem_crack = _integral_values[i] > _critical_energy_release_rate;
    std::cout << "critical_energy_release_rate = " << _critical_energy_release_rate << ", weibull_tip = " << _weibull_at_tip[i] << " does elem crack = " << does_elem_crack << std::endl;
    _xfem->update_doesElemCrackEnergyReleaseRate(_elem_id_crack_tip[i], does_elem_crack);
  }
}
