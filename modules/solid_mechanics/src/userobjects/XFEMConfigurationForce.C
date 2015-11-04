/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This userobject calculates the configurational force
//
#include "XFEMConfigurationForce.h"
#include "libmesh/fe_interface.h"

template<>
InputParameters validParams<XFEMConfigurationForce>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<Real>("radius_inner", "Inner radius for volume integral domain");
  params.addParam<Real>("radius_outer", "Outer radius for volume integral domain");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

XFEMConfigurationForce::XFEMConfigurationForce(const InputParameters & parameters):
    ElementUserObject(parameters),
    _Eshelby_tensor(getMaterialProperty<ColumnMajorMatrix>("Eshelby_tensor")),
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
XFEMConfigurationForce::initialize()
{ 
  _crack_front_points.clear();
  _elem_id_crack_tip.clear();
  _integral_values.clear();

  _num_crack_front_points = _xfem->num_crack_tips();

  _integral_values.resize(_num_crack_front_points*3);
  
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    _integral_values[i] = 0.0;

  _xfem->get_crack_tip_origin(_elem_id_crack_tip, _crack_front_points);
  std::cout << "number of crack front points = " << _num_crack_front_points << std::endl;

  for (unsigned int i = 0; i < _num_crack_front_points; i ++)
  {
    std::cout << "crack front points[" << i << "] = " << _crack_front_points[i] << std::endl;
  }
}

std::vector<Real>
XFEMConfigurationForce::computeIntegrals()
{
  FEType fe_type(Utility::string_to_enum<Order>("first"),Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule (_qrule);

  // The values of the shape functions at the quadrature points
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  fe->reinit (_current_elem); 
  
  std::vector<Real> sums(_num_crack_front_points*3);
  for (signed int i = 0; i < _num_crack_front_points*3; i++)
    sums[i] = 0.0;

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
  {
    std::vector<Real> QpIntegrals = computeQpIntegrals(phi,dphi);
    for (unsigned int i = 0; i < _num_crack_front_points*3; i++){
      sums[i] += QpIntegrals[i] * _JxW[_qp]*_coord[_qp];
    }
  }
  return sums;
}

Real
XFEMConfigurationForce::calcQValue(Point & node, Point & crack_front)
{
  Point dist_to_crack_front_vector = node - crack_front;
  Real dist_to_crack_front = std::pow(dist_to_crack_front_vector.size_sq(),0.5);
  Real q = 1.0;
  //if ( dist_to_crack_front > _radius_inner &&
  //     dist_to_crack_front < _radius_outer){
  //  q = (_radius_outer - dist_to_crack_front) /
  //    (_radius_outer - _radius_inner);
  //}
  //else if ( dist_to_crack_front >= _radius_outer)
  //  q = 0.0;

  if(dist_to_crack_front > _radius_outer)
    q = -1.0;
  else
    q = 1.0;
  
  return q;
}

std::vector<Real>
XFEMConfigurationForce::computeQpIntegrals(const std::vector<std::vector<Real> > & N_shape_func, const std::vector<std::vector<RealGradient> > & dN_shape_func)
{
  std::vector<Real> QpIntegrals(_num_crack_front_points*3);
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    QpIntegrals[i] = 0.0;


  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    unsigned int n_nodes = _current_elem->n_nodes();
    std::vector<Real> q_nodes(n_nodes,0.0);
    RealVectorValue grad_of_scalar_q(0.0,0.0,0.0);

    Point crack_front = _crack_front_points[i];

    // calculate Q function at finite element node
    /*
    for (unsigned int i = 0; i < n_nodes; i++)
    {
      Real q = calcQValue(*(_current_elem->get_node(i)), crack_front);
      q_nodes[i] = q;
    }

    // calcuate the Q function and its gradient at quadrature point
    for (unsigned int i = 0; i < n_nodes; i++)
    {
      grad_of_scalar_q(0) += q_nodes[i] * dN_shape_func[i][_qp](0);
      grad_of_scalar_q(1) += q_nodes[i] * dN_shape_func[i][_qp](1);
      grad_of_scalar_q(2) += q_nodes[i] * dN_shape_func[i][_qp](2);
    }
  
    ColumnMajorMatrix Jvec = _Eshelby_tensor[_qp]*grad_of_scalar_q;

    QpIntegrals[i*3]   = Jvec(0,0);
    QpIntegrals[i*3+1] = Jvec(1,0);
    QpIntegrals[i*3+2] = Jvec(2,0);
  */
    ColumnMajorMatrix Jvec(3,1);
    Jvec.zero();
    for (unsigned int i = 0; i < n_nodes; i++)
    {
      Real q = calcQValue(*(_current_elem->get_node(i)), crack_front);
      RealVectorValue grad_of_N(0.0,0.0,0.0);
      if (q > 0.0)
      {
        grad_of_N(0) = dN_shape_func[i][_qp](0);
        grad_of_N(1) = dN_shape_func[i][_qp](1);
        grad_of_N(2) = dN_shape_func[i][_qp](2);
        Jvec += _Eshelby_tensor[_qp]*grad_of_N;
      }
    }

    QpIntegrals[i*3]   = Jvec(0,0);
    QpIntegrals[i*3+1] = Jvec(1,0);
    QpIntegrals[i*3+2] = Jvec(2,0);
  }

  return QpIntegrals;
}

void 
XFEMConfigurationForce::execute()
{
  std::vector<Real> comp_integ = computeIntegrals();
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    _integral_values[i] += comp_integ[i];
}

void
XFEMConfigurationForce::threadJoin(const UserObject & y)
{
  const XFEMConfigurationForce & pps = static_cast<const XFEMConfigurationForce &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    _integral_values[i] += pps._integral_values[i];
}

void 
XFEMConfigurationForce::finalize()
{
  //_communicator.set_union(_elem_id_crack_tip);

  _xfem->clear_crack_propagation_direction();
  gatherSum(_integral_values);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Point direction(_integral_values[i*3], _integral_values[i*3+1], _integral_values[i*3+2]);
    std::cout << "direction = " << direction << std::endl;
    direction /= pow(direction.size_sq(),0.5);
    direction *= -1.0; // crack propagations in the direction of the inverse of the crack-driving force
    _xfem->update_crack_propagation_direction(_elem_id_crack_tip[i], direction);
    std::cout << "crack front index (" << i << ") : configuration force  = " << direction << std::endl; 
  } 
}

