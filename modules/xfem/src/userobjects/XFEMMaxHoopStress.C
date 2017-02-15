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

template<>
InputParameters validParams<XFEMMaxHoopStress>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<Real>("radius_inner", "Inner radius for volume integral domain");
  params.addParam<Real>("radius_outer", "Outer radius for volume integral domain");
  params.addParam<Real>("thermal_expansion", 0.0, "Thermal expansion");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>("poissons_ratio", "Poisson's ratio for the material.");
  params.addParam<Real>("youngs_modulus", "Young's modulus of the material.");
  params.addCoupledVar("disp_x", "The x displacement");
  params.addCoupledVar("disp_y", "The y displacement");
  params.addCoupledVar("disp_z", "The z displacement");
  params.addCoupledVar("temp", "Coupled Temperature");
  params.addParam<BoundaryName>("intersecting_boundary", "Boundary intersected by ends of crack.");
  params.addParam<PostprocessorName>("average_h", "Postprocessor that gives average element size");
  return params;
}

XFEMMaxHoopStress::XFEMMaxHoopStress(const InputParameters & parameters):
    ElementUserObject(parameters),
    _Eshelby_tensor(getMaterialProperty<ColumnMajorMatrix>("Eshelby_tensor")),
    _J_thermal_term_vec(hasMaterialProperty<RealVectorValue>("J_thermal_term_vec")?
                        &getMaterialProperty<RealVectorValue>("J_thermal_term_vec"):
                        NULL),
    _stress(getMaterialPropertyByName<SymmTensor>("stress")),
    _strain(getMaterialPropertyByName<SymmTensor>("elastic_strain")),
    _grad_disp_x(coupledGradient("disp_x")),
    _grad_disp_y(coupledGradient("disp_y")),
    _grad_disp_z(parameters.get<SubProblem *>("_subproblem")->mesh().dimension() == 3 ? coupledGradient("disp_z") : _grad_zero), 
    _temp_grad(coupledGradient("temp")),
    _qp(0),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _mesh(_subproblem.mesh()),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _youngs_modulus(getParam<Real>("youngs_modulus")),
    _postprocessor( isParamValid("average_h") ? &getPostprocessorValue("average_h") : NULL )
{
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
 
  //plane strain
  _kappa = 3 - 4*_poissons_ratio;
  _shear_modulus = _youngs_modulus / (2*(1 + _poissons_ratio));

  _K_factor = 0.5 * _youngs_modulus / (1 - std::pow(_poissons_ratio,2));

  _sif_mode.push_back(SIF_MODE(0));
  _sif_mode.push_back(SIF_MODE(1));
}

void
XFEMMaxHoopStress::computeAuxFields(const SIF_MODE sif_mode, ColumnMajorMatrix & stress, ColumnMajorMatrix & disp, ColumnMajorMatrix & grad_disp, ColumnMajorMatrix & strain)
{

  RealVectorValue k(0);
  if (sif_mode == KI)
    k(0) = 1;

  else if (sif_mode == KII)
    k(1) = 1;

  Real t = _theta;
  Real t2 = _theta/2;
  Real tt2 = 3*_theta/2;
  Real st = std::sin(t);
  Real ct = std::cos(t);
  Real st2 = std::sin(t2);
  Real ct2 = std::cos(t2);
  Real stt2 = std::sin(tt2);
  Real ctt2 = std::cos(tt2);
  Real ctsq = std::pow(ct,2);
  Real ct2sq = std::pow(ct2,2);
  Real ct2cu = std::pow(ct2,3);
  Real sqrt2PiR = std::sqrt(2*libMesh::pi*_r);

  //Calculate auxiliary stress tensor
  Real s11 = 1 / sqrt2PiR * ( k(0) * ct2 * (1 - st2 * stt2)
                                  - k(1) * st2 * (2 + ct2 * ctt2));
  Real s22 = 1 / sqrt2PiR * ( k(0) * ct2 * (1 + st2 * stt2)
                                  + k(1) * st2 * ct2 * ctt2);
  Real s12 = 1 / sqrt2PiR * ( k(0) * ct2 * st2 * ctt2
                                  + k(1) * ct2 * (1 - st2*stt2));
  Real s13 = -1 / sqrt2PiR * k(2) * st2;
  Real s23 =  1 / sqrt2PiR * k(2) * ct2;
  //plain stress
  //Real s33 = 0;
  //plain strain
  Real s33 = _poissons_ratio * (s11 + s22);

  stress(0,0) = s11;
  stress(0,1) = s12;
  stress(0,2) = s13;
  stress(1,0) = s12;
  stress(1,1) = s22;
  stress(1,2) = s23;
  stress(2,0) = s13;
  stress(2,1) = s23;
  stress(2,2) = s33;

  //Calculate x1 derivative of auxiliary stress tensor
  Real ds111 = k(0) / (2*_r*sqrt2PiR)
    * ( - ct*ct2 + ct*ct2*st2*stt2 + st*st2 - st*stt2
        + 2*st*ct2*ct2*stt2 + 3*st*ct2*st2*ctt2 )
    + k(1) / (2*_r*sqrt2PiR)
    * ( 2*st2*ct + ct*st2*ct2*ctt2 + 2*st*ct2 - st*ctt2
        + 2*st*ct2*ct2*ctt2 - 3*st*st2*ct2*stt2 );
  Real ds121 = k(0) / (2*_r*sqrt2PiR)
    * ( - ct*ct2*st2*ctt2 + st*ctt2
        - 2*st*ct2*ct2*ctt2 + 3*st*st2*ct2*stt2 )
    + k(1) / (2*_r*sqrt2PiR)
    * ( - ct*ct2 + ct*ct2*st2*stt2 + st*st2 - st*stt2
        + 2*st*ct2*ct2*stt2 + 3*st*ct2*st2*ctt2 );
  Real ds131 = k(2) / (2*_r*sqrt2PiR) * ( st2*ct + ct2*st);
  Real ds221 = k(0) / (2*_r*sqrt2PiR)
    * ( - ct*ct2 - ct*ct2*st2*stt2 + st*st2 + st*stt2
        - 2*st*ct2*ct2*stt2 - 3*st*ct2*st2*ctt2 )
    + k(1) / (2*_r*sqrt2PiR)
    * ( - ct*ct2*st2*ctt2 + st*ctt2
        - 2*st*ct2*ct2*ctt2 + 3*st*st2*ct2*stt2 );
  Real ds231 = k(2) / (2*_r*sqrt2PiR) * ( - ct2*ct + st2*st );
  Real ds331 = _poissons_ratio * (ds111 + ds221);

  //Calculate auxiliary displacements
  disp(0,0) = 1 / (2*_shear_modulus) * std::sqrt(_r/(2*libMesh::pi)) * (k(0) * ct2 * (_kappa - 1 + 2*st2*st2)
                                                      + k(1) * st2 * (_kappa + 1 + 2*ct2*ct2));
  disp(0,1) = 1 / (2*_shear_modulus) * std::sqrt(_r/(2*libMesh::pi)) * (k(0) * st2 * (_kappa + 1 - 2*ct2*ct2)
                                                      - k(1) * ct2 * (_kappa - 1 - 2*st2*st2));
  disp(0,2) = 1 /   _shear_modulus * std::sqrt(_r/(2*libMesh::pi)) * k(2) * st2*st2;

  //Calculate x1 derivative of auxiliary displacements
  Real du11 = k(0)/(4 * _shear_modulus * sqrt2PiR)
    * ( ct*ct2*_kappa + ct*ct2 - 2*ct*ct2cu
        + st*st2*_kappa + st*st2 - 6*st*st2*ct2sq )
    + k(1)/(4 * _shear_modulus * sqrt2PiR)
    * ( ct*st2*_kappa + ct*st2 + 2*ct*st2*ct2sq
        - st*ct2*_kappa + 3*st*ct2 - 6*st*ct2cu );

  Real du21 =  k(0)/(4 * _shear_modulus * sqrt2PiR)
    * ( ct*st2*_kappa + ct*st2 - 2*ct*st2*ct2sq
        - st*ct2*_kappa - 5*st*ct2 + 6*st*ct2cu )
    + k(1)/(4 * _shear_modulus * sqrt2PiR)
    * ( - ct*ct2*_kappa + 3*ct*ct2 - 2*ct*ct2cu
        - st*st2*_kappa + 3*st*st2 - 6*st*st2*ct2sq );

  Real du31 = k(2) / (_shear_modulus * sqrt2PiR)
    * ( st2*ct - ct2*st);

  grad_disp(0,0) = du11;
  grad_disp(0,1) = du21;
  grad_disp(0,2) = du31;

  //Calculate second derivatives of displacements (u,1i)
  //only needed for inhomogenous materials
  Real du111 = k(0)/(8 * _shear_modulus * _r * sqrt2PiR)
    * ( - 2*ctsq*ct2*_kappa
        + ct2*_kappa
        + 10*ctsq*ct2 - 11*ct2
        - 12*ctsq*ct2cu
        + 14*ct2cu
        - 2*ct*st2*st*_kappa
        - 2*ct*st2*st
        + 12*ct*st2*st*ct2sq )
    + k(1)/(8 * _shear_modulus * _r * sqrt2PiR )
    * ( - 2*ctsq*st2*_kappa
        - 6*ctsq*st2
        + 12*ctsq*st2*ct2sq
        + 2*ct*ct2*st*_kappa
        - 6*ct*ct2*st
        + 12*ct*ct2cu*st
        + st2*_kappa + 5*st2
        - 14*st2*ct2sq );
  Real du112 = k(0)/(8 * _shear_modulus * _r * sqrt2PiR)
    * ( - 2*ct*ct2*st*_kappa
        + 10*ct*ct2*st
        - 12*st*ct*ct2cu
        - st2*_kappa
        + 2*ctsq*st2*_kappa
        - st2 + 2*ctsq*st2
        + 6*st2*ct2sq
        - 12*ctsq*st2*ct2sq )
    + k(2)/(8 * _shear_modulus * _r * sqrt2PiR )
    * ( - 2*ct*st2*st*_kappa
        - 6*ct*st2*st
        + 12*ct*st2*st*ct2sq
        + ct2*_kappa + 6*ctsq*ct2
        - 2*ctsq*ct2*_kappa
        - 3*ct2 + 6*ct2cu
        - 12*ctsq*ct2cu );
  Real du113 = 0.0;

  Real du211 = k(0)/(8 * _shear_modulus * _r * sqrt2PiR)
    * ( - 2*ctsq*ct2*_kappa
        + ct2*_kappa
        + 10*ctsq*ct2 - 11*ct2
        - 12*ctsq*ct2cu
        + 14*ct2cu
        - 2*ct*st2*st*_kappa
        - 2*ct*st2*st
        + 12*ct*st2*st*ct2sq )
    + k(1)/(8 * _shear_modulus * _r * sqrt2PiR )
    * ( - 2*ctsq*st2*_kappa
        - 6*ctsq*st2
        + 12*ctsq*st2*ct2sq
        + 2*ct*ct2*st*_kappa
        - 6*ct*ct2*st
        + 12*ct*ct2cu*st
        + st2*_kappa + 5*st2
        - 14*st2*ct2sq );

  Real du212 = k(0)/(8 * _shear_modulus * _r * sqrt2PiR)
    *( - 2*ct*st2*st*_kappa
       + 2*ct*st2*st - 6*ct2cu
       - 12*ct*st2*st*ct2sq
       + ct2*_kappa
       - 2*ctsq*ct2*_kappa
       + 5*ct2 - 10*ctsq*ct2
       + 12*ctsq*ct2cu )
    + k(2)/(8 * _shear_modulus * _r * sqrt2PiR )
    * (   2*ct*ct2*st*_kappa
          + 6*ct*ct2*st
          + st2*_kappa
          - 2*ctsq*st2*_kappa
          - 3*st2 + 6*ctsq*st2
          + 6*st2*ct2sq
          - 12*ctsq*st2*ct2sq
          - 12*ct*ct2cu*st );

  Real du213 = 0.0;

  Real du311 = k(2)/(2 * _shear_modulus * _r * sqrt2PiR)
    * ( - 2*ctsq*st2 + st2
        + 2*ct*ct2*st );

  Real du312 = k(2)/(2 * _shear_modulus * _r * sqrt2PiR)
    * ( - 2*ct*st2*st + ct2
        - 2*ctsq*ct2 );

  Real du313 = 0.0;

  //Calculate auxiliary strains
  strain(0,0) = (1 / _youngs_modulus) * (s11 - _poissons_ratio * (s22 + s33));
  strain(1,1) = (1 / _youngs_modulus) * (s22 - _poissons_ratio * (s11 + s33));
  strain(2,2) = (1 / _youngs_modulus) * (s33 - _poissons_ratio * (s11 + s22));
  strain(0,1) = (1 / _shear_modulus)  * s12;
  strain(1,0) = (1 / _shear_modulus)  * s12;
  strain(1,2) = (1 / _shear_modulus)  * s23;
  strain(2,1) = (1 / _shear_modulus)  * s23;
  strain(0,2) = (1 / _shear_modulus)  * s13;
  strain(2,0) = (1 / _shear_modulus)  * s13;
}

void
XFEMMaxHoopStress::initialize()
{
  _crack_front_points.clear();
  _crack_front_directions.clear();
  _elem_id_crack_tip.clear();
  _integral_values.clear();

  _num_crack_front_points = _xfem->numberCrackTips();

  _integral_values.resize(_num_crack_front_points*2);
  
  for (unsigned int i = 0; i < _num_crack_front_points*2; i++)
    _integral_values[i] = 0.0;

  _xfem->getCrackTipOriginDirection(_elem_id_crack_tip, _crack_front_points, _crack_front_directions);

  for (unsigned int i = 0; i <  _crack_front_points.size(); i++)
  {
    std::cout << "MAX Hoop Stress: crack_front_points " << _crack_front_points[i] << std::endl; 
    std::cout << "MAX Hoop Stress: crack_front_directions " << _crack_front_directions[i] << std::endl; 
  }
}

std::vector<Real>
XFEMMaxHoopStress::computeIntegrals()
{
  FEType fe_type(Utility::string_to_enum<Order>("first"),Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule (_qrule);

  // The values of the shape functions at the quadrature points
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  fe->reinit (_current_elem); 
  
  std::vector<Real> sums(_num_crack_front_points*2);
  for (unsigned int i = 0; i < _num_crack_front_points*2; i++)
    sums[i] = 0.0;

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
  {
    std::vector<Real> QpIntegrals = computeQpIntegrals(phi,dphi);
    for (unsigned int i = 0; i < _num_crack_front_points*2; i++){
      sums[i] += QpIntegrals[i] * _JxW[_qp]*_coord[_qp];
    }
  }
  return sums;
}

Real
XFEMMaxHoopStress::calcQValue(Point & node, Point & crack_front)
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

RealVectorValue
XFEMMaxHoopStress::rotateToCrackFrontCoords(const RealVectorValue vector) const
{
  ColumnMajorMatrix vec3x1;
  vec3x1 = _rot_mat * vector;
  RealVectorValue vec;
  vec(0) = vec3x1(0,0);
  vec(1) = vec3x1(1,0);
  vec(2) = vec3x1(2,0);
  return vec;
}

RealVectorValue
XFEMMaxHoopStress::rotateToCrackFrontCoords(const Point point) const
{
  RealVectorValue vector(point(0), point(1), point(2));
  ColumnMajorMatrix vec3x1;
  vec3x1 = _rot_mat * vector;
  RealVectorValue vec;
  vec(0) = vec3x1(0,0);
  vec(1) = vec3x1(1,0);
  vec(2) = vec3x1(2,0);
  return vec;
}


ColumnMajorMatrix
XFEMMaxHoopStress::rotateToCrackFrontCoords(const SymmTensor tensor) const
{
  ColumnMajorMatrix tensor_CMM;
  tensor_CMM(0,0) = tensor.xx();
  tensor_CMM(0,1) = tensor.xy();
  tensor_CMM(0,2) = tensor.xz();
  tensor_CMM(1,0) = tensor.xy();
  tensor_CMM(1,1) = tensor.yy();
  tensor_CMM(1,2) = tensor.yz();
  tensor_CMM(2,0) = tensor.xz();
  tensor_CMM(2,1) = tensor.yz();
  tensor_CMM(2,2) = tensor.zz();

  ColumnMajorMatrix tmp = _rot_mat * tensor_CMM;
  ColumnMajorMatrix rotT = _rot_mat.transpose();
  ColumnMajorMatrix rotated_tensor = tmp * rotT;

  return rotated_tensor;
}

ColumnMajorMatrix
XFEMMaxHoopStress::rotateToCrackFrontCoords(const ColumnMajorMatrix tensor) const
{
  ColumnMajorMatrix tmp = _rot_mat * tensor;
  ColumnMajorMatrix rotT = _rot_mat.transpose();
  ColumnMajorMatrix rotated_tensor = tmp * rotT;

  return rotated_tensor;
}


void 
XFEMMaxHoopStress::calcRTheta(Point & p, Point & crack_front_point, Point & crack_direction)
{
  RealVectorValue tangent_direction;
  tangent_direction(2) = 1.0;
  _crack_plane_normal = tangent_direction.cross(crack_direction);
  _rot_mat(0,0) = crack_direction(0);
  _rot_mat(0,1) = crack_direction(1);
  _rot_mat(0,2) = crack_direction(2);
  _rot_mat(1,0) = _crack_plane_normal(0);
  _rot_mat(1,1) = _crack_plane_normal(1);
  _rot_mat(1,2) = _crack_plane_normal(2);
  _rot_mat(2,0) = 0.0;
  _rot_mat(2,1) = 0.0;
  _rot_mat(2,2) = 0.0;
  _rot_mat(2,2) = 1.0;

  Point closest_point(0.0);
  RealVectorValue crack_front_point_rot = rotateToCrackFrontCoords(crack_front_point);

  RealVectorValue crack_front_edge = rotateToCrackFrontCoords(tangent_direction);

  Point p_rot = rotateToCrackFrontCoords(p);
  p_rot = p_rot - crack_front_point_rot;

  RealVectorValue closest_point_to_p = p_rot;

  //Find r, the distance between the qp and the crack front
  RealVectorValue r_vec = p_rot;
  _r = r_vec.size();

  //Find theta, the angle between r and the crack front plane
  RealVectorValue crack_plane_normal = rotateToCrackFrontCoords(_crack_plane_normal);
  Real p_to_plane_dist = std::abs(closest_point_to_p*crack_plane_normal);

  //Determine if p is above or below the crack plane
  Real y_local = p_rot(1) - closest_point(1);

  //Determine if p is in front of or behind the crack front
  RealVectorValue p2(p_rot);
  p2(1) = 0;
  RealVectorValue p2_vec = p2 - closest_point;
  Real ahead = crack_front_edge(2) * p2_vec(0) - crack_front_edge(0) * p2_vec(2);

  Real x_local(0);
  if (ahead >= 0)
    x_local = 1;
  else
    x_local = -1;

  //Calculate theta based on in which quadrant in the crack front coordinate
  //system the qp is located
  if (_r > 0)
  {
    Real theta_quadrant1(0.0);
    if (MooseUtils::absoluteFuzzyEqual(_r, p_to_plane_dist, 1e-10))
      theta_quadrant1 = 0.5*libMesh::pi;
    else if (p_to_plane_dist > _r)
      mooseError("Invalid distance p_to_plane_dist in CrackFrontDefinition::calculateRThetaToCrackFront");
    else
      theta_quadrant1 = std::asin(p_to_plane_dist/_r);

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
XFEMMaxHoopStress::computeQpIntegrals(const std::vector<std::vector<Real> > & N_shape_func, const std::vector<std::vector<RealGradient> > & dN_shape_func)
{
  std::vector<Real> QpIntegrals(_num_crack_front_points*2);
  for (unsigned int i = 0; i < _num_crack_front_points*2; i++)
    QpIntegrals[i] = 0.0;


  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    unsigned int n_nodes = _current_elem->n_nodes();
    std::vector<Real> q_nodes(n_nodes,0.0);
    RealVectorValue grad_of_scalar_q(0.0,0.0,0.0);
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

    //Calculate auxiliary stress tensor at current qp
    for (std::vector<SIF_MODE>::iterator it = _sif_mode.begin(); it != _sif_mode.end(); ++it)
    {
      if (*it == KI)
        computeAuxFields(*it, _aux_stress, _aux_disp, _aux_grad_disp, _aux_strain);

      else if (*it == KII)
        computeAuxFields(*it, _aux_stress, _aux_disp, _aux_grad_disp, _aux_strain);


      //In the crack front coordinate system, the crack direction is (1,0,0)
      RealVectorValue crack_direction_local(0.0);
      crack_direction_local(0) = 1.0;

      ColumnMajorMatrix aux_du;
      aux_du(0,0) = _aux_grad_disp(0,0);
      aux_du(0,1) = _aux_grad_disp(0,1);
      aux_du(0,2) = _aux_grad_disp(0,2);

      ColumnMajorMatrix stress;
      stress(0,0) = _stress[_qp].xx();
      stress(0,1) = _stress[_qp].xy();
      stress(0,2) = _stress[_qp].xz();
      stress(1,0) = _stress[_qp].xy();
      stress(1,1) = _stress[_qp].yy();
      stress(1,2) = _stress[_qp].yz();
      stress(2,0) = _stress[_qp].xz();
      stress(2,1) = _stress[_qp].yz();
      stress(2,2) = _stress[_qp].zz();

      ColumnMajorMatrix strain;
      strain(0,0) = _strain[_qp].xx();
      strain(0,1) = _strain[_qp].xy();
      strain(0,2) = _strain[_qp].xz();
      strain(1,0) = _strain[_qp].xy();
      strain(1,1) = _strain[_qp].yy();
      strain(1,2) = _strain[_qp].yz();
      strain(2,0) = _strain[_qp].xz();
      strain(2,1) = _strain[_qp].yz();
      strain(2,2) = _strain[_qp].zz();

      ColumnMajorMatrix grad_disp;
      grad_disp(0,0) = _grad_disp_x[_qp](0);
      grad_disp(0,1) = _grad_disp_x[_qp](1);
      grad_disp(0,2) = _grad_disp_x[_qp](2);
      grad_disp(1,0) = _grad_disp_y[_qp](0);
      grad_disp(1,1) = _grad_disp_y[_qp](1);
      grad_disp(1,2) = _grad_disp_y[_qp](2);
      grad_disp(2,0) = _grad_disp_z[_qp](0);
      grad_disp(2,1) = _grad_disp_z[_qp](1);
      grad_disp(2,2) = _grad_disp_z[_qp](2);

      //Rotate stress, strain, and displacement to crack front coordinate system
      RealVectorValue grad_q_cf = rotateToCrackFrontCoords(grad_q);
      ColumnMajorMatrix grad_disp_cf = rotateToCrackFrontCoords(grad_disp);
      ColumnMajorMatrix stress_cf = rotateToCrackFrontCoords(stress);
      ColumnMajorMatrix strain_cf = rotateToCrackFrontCoords(strain);

      ColumnMajorMatrix dq;
      dq(0,0) = crack_direction_local(0)*grad_q_cf(0);
      dq(0,1) = crack_direction_local(0)*grad_q_cf(1);
      dq(0,2) = crack_direction_local(0)*grad_q_cf(2);

      //Calculate interaction integral terms

      // Term1 = stress * x1-derivative of aux disp * dq
      ColumnMajorMatrix tmp1 = dq * stress_cf;
      Real term1 = aux_du.doubleContraction(tmp1);

      // Term2 = aux stress * x1-derivative of disp * dq
      ColumnMajorMatrix tmp2 = dq * _aux_stress;
      Real term2 = grad_disp_cf(0,0)*tmp2(0,0)+grad_disp_cf(1,0)*tmp2(0,1)+grad_disp_cf(2,0)*tmp2(0,2);

      // Term3 = aux stress * strain * dq_x   (= stress * aux strain * dq_x)
      Real term3 = dq(0,0) * _aux_stress.doubleContraction(strain_cf);

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

      RealVectorValue grad_temp_cf = rotateToCrackFrontCoords(_temp_grad[_qp]);

      Real eq_thermal = 0.0;
      Real aux_stress_trace = _aux_stress(0,0) + _aux_stress(1,1) + _aux_stress(2,2);
      eq_thermal = scalar_q * aux_stress_trace * _thermal_expansion * grad_temp_cf(0);

      Real eq = term1 + term2 - term3 + eq_thermal;

      QpIntegrals[i*2 + int(*it)]   = eq;
    }
  }
  return QpIntegrals;
}

void 
XFEMMaxHoopStress::execute()
{
  if (_postprocessor)
  {
    _radius_inner = 1.5 *  *_postprocessor;
    _radius_outer = 3.5 *  *_postprocessor;
  } 

  std::vector<Real> comp_integ = computeIntegrals();
  for (unsigned int i = 0; i < _num_crack_front_points*2; i++)
    _integral_values[i] += comp_integ[i];
}

void
XFEMMaxHoopStress::threadJoin(const UserObject & y)
{
  const XFEMMaxHoopStress & pps = static_cast<const XFEMMaxHoopStress &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points*2; i++)
    _integral_values[i] += pps._integral_values[i];
}

void 
XFEMMaxHoopStress::finalize()
{
  _xfem->clearCrackGrowthDirection();
  gatherSum(_integral_values);

  for (unsigned int i = 0; i < _num_crack_front_points*2; i++)
    _integral_values[i] *= _K_factor;

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Real KI = _integral_values[i*2];
    Real KII = _integral_values[i*2+1];
    std::cout << "KI = " << KI << std::endl;
    std::cout << "KII = " << KII << std::endl;
    Real theta1 = 2*std::atan(0.25*(KI/KII + std::sqrt(pow(KI/KII,2.0) + 8.0)));
    Real theta2 = 2*std::atan(0.25*(KI/KII - std::sqrt(pow(KI/KII,2.0) + 8.0)));

    Real hoop_stress1 = KI * (3 * std::cos(theta1*0.5) + std::cos(theta1*1.5)) + KII * (-3.0*std::sin(theta1*0.5) - 3.0 * std::sin(1.5*theta1)); 
    Real hoop_stress2 = KI * (3 * std::cos(theta2*0.5) + std::cos(theta2*1.5)) + KII * (-3.0*std::sin(theta2*0.5) - 3.0 * std::sin(1.5*theta2)); 

    Real theta = 0.0;

    std::cout << "hoop_stress1 = " << hoop_stress1 << ", theta1 = " << theta1 / libMesh::pi * 180.0 << std::endl;
    std::cout << "hoop_stress2 = " << hoop_stress2 << ", theta2 = " << theta2 / libMesh::pi * 180.0 << std::endl;

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

    Point direction(std::cos(omega+theta), std::sin(omega+theta), 0.0);
    
    _xfem->updateCrackGrowthDirection(_elem_id_crack_tip[i], direction);
    std::cout << "MAXHOOPSTRESS crack front point = " << crack_front_point << std::endl;
    std::cout << "MAXHOOPSTRESS crack front index (" << i << ") : direction  = " << direction << std::endl; 
  } 
}

