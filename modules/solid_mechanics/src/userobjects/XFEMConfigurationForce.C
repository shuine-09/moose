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
#include "DisplacedProblem.h"

template<>
InputParameters validParams<XFEMConfigurationForce>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<Real>("radius", "Radius for volume integral domain");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

XFEMConfigurationForce::XFEMConfigurationForce(const InputParameters & parameters):
    ElementUserObject(parameters),
    _Eshelby_tensor(getMaterialProperty<ColumnMajorMatrix>("Eshelby_tensor")),
    _J_thermal_term_vec(hasMaterialProperty<RealVectorValue>("J_thermal_term_vec")?
                        &getMaterialProperty<RealVectorValue>("J_thermal_term_vec"):
                        NULL),
    _qp(0),
    _mesh(_subproblem.mesh())
{
  _fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (_fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMarkerUserObject");
  _xfem = _fe_problem->get_xfem();

  if (isParamValid("radius"))
    _radius = getParam<Real>("radius");
  else
    mooseError("DomainIntegral error: must set radius.");
}

void
XFEMConfigurationForce::initialize()
{
  _J_thermal_term_vec = hasMaterialProperty<RealVectorValue>("J_thermal_term_vec") ? &getMaterialProperty<RealVectorValue>("J_thermal_term_vec"):NULL;
 
  _crack_front_points.clear();
  _elem_id_crack_tip.clear();
  _integral_values.clear();

  _num_crack_front_points = _xfem->num_crack_tips();

  _integral_values.resize(_num_crack_front_points*3);
  
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    _integral_values[i] = 0.0;

  _xfem->get_crack_tip_origin(_elem_id_crack_tip, _crack_front_points);

  for (unsigned int i = 0; i <  _crack_front_points.size(); i++)
    std::cout << "CONFIGURATIONFORCE: crack_front_points " << _crack_front_points[i] << std::endl; 
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
  Real q = 0.0;
  if(dist_to_crack_front > _radius)
    q = -1.0;
  else
    q = 1.0;
  
  return q;
}

bool
XFEMConfigurationForce::isWithin(Point & crack_front)
{
  bool is_within = false;
  unsigned int n_nodes = _current_elem->n_nodes();
  Real pi;
  for (unsigned int i = 0; i < n_nodes; i++)
  { 
    pi = calcQValue(*(_current_elem->get_node(i)), crack_front);
    if (pi > 0.0)
    {
      is_within = true;
      return is_within;
    }
  }
  return is_within;
}


bool
XFEMConfigurationForce::isIntersect(Point & crack_front)
{
  bool is_intersect = false;
  unsigned int n_nodes = _current_elem->n_nodes();
  Real pi, pj;
  for (unsigned int i = 0; i < n_nodes; i++)
  { 
    pi = calcQValue(*(_current_elem->get_node(i)), crack_front);
    for (unsigned int j = i; j < n_nodes; j++)
    { 
      pj = calcQValue(*(_current_elem->get_node(j)), crack_front);
      if (pi*pj < 0.0)
      {
        is_intersect = true;
        return is_intersect;
      }
    } 
  }
  return is_intersect;
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

    ColumnMajorMatrix Jvec(3,1);
    Jvec.zero();

    //Contour approach
    //if ( ! isIntersect(crack_front))
    //  continue;

    //Domain approach
    if (!isWithin(crack_front))
      continue;

    //const Elem * undisplaced_elem  = NULL;
    //if(_fe_problem->getDisplacedProblem() != NULL)
    //  undisplaced_elem = _fe_problem->getDisplacedProblem()->refMesh().elem(_current_elem->id());
    //else
    //  undisplaced_elem = _current_elem;

    for (unsigned int i = 0; i < n_nodes; i++)
    {
      //Real flag = _xfem->flag_qp_inside(undisplaced_elem, *(_current_elem->get_node(i))); //inside (flag = 1) or ouside (flag = 0) real domain
      //if ( flag < 0.5)
      //  continue;

      Real q = calcQValue(*(_current_elem->get_node(i)), crack_front);
      RealVectorValue grad_of_N(0.0,0.0,0.0);
      if (q > 0.0)
      {
        grad_of_N(0) = dN_shape_func[i][_qp](0);
        grad_of_N(1) = dN_shape_func[i][_qp](1);
        grad_of_N(2) = dN_shape_func[i][_qp](2);
        Jvec += _Eshelby_tensor[_qp]*grad_of_N;
        //Thermal component
        if (_J_thermal_term_vec)
        {
          Jvec(0,0) -= (*_J_thermal_term_vec)[_qp](0) * N_shape_func[i][_qp];
          Jvec(1,0) -= (*_J_thermal_term_vec)[_qp](1) * N_shape_func[i][_qp];
          Jvec(2,0) -= (*_J_thermal_term_vec)[_qp](2) * N_shape_func[i][_qp];
        }
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
  _xfem->clear_crack_propagation_direction();
  gatherSum(_integral_values);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Point direction(_integral_values[i*3], _integral_values[i*3+1], _integral_values[i*3+2]);
    std::cout << "direction = " << direction << std::endl;
    if (direction.size_sq() > 1.0e-20)
      direction /= pow(direction.size_sq(),0.5);
    direction *= -1.0; // crack propagations in the direction of the inverse of the crack-driving force
    _xfem->update_crack_propagation_direction(_elem_id_crack_tip[i], direction);
    std::cout << "crack front index (" << i << ") : configuration force  = " << direction << std::endl; 
  } 
}

