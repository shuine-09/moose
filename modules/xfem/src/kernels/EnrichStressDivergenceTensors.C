/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "EnrichStressDivergenceTensors.h"
#include "Material.h"
#include "MooseMesh.h"
#include "ElasticityTensorTools.h"
#include "libmesh/quadrature.h"
#include "XFEM.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"


template<>
InputParameters validParams<EnrichStressDivergenceTensors>()
{
  InputParameters params = validParams<ALEKernel>();
  params.addClassDescription("Enrich stress divergence kernel for the Cartesian coordinate system");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<unsigned int>("enrichment_component", "The component of the enrichement functions");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("enrichment_displacement", "The string of displacements suitable for the problem statement");
  params.addRequiredCoupledVar("enrichment_displacement_var", "The string of displacements suitable for the problem statement");
  params.addParam<std::string>("base_name", "Material property base name");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}


EnrichStressDivergenceTensors::EnrichStressDivergenceTensors(const InputParameters & parameters) :
    ALEKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _component(getParam<unsigned int>("component")),
    _enrichment_component(getParam<unsigned int>("enrichment_component")),
    _nl_vnames(getParam<std::vector<NonlinearVariableName> >("enrichment_displacement"))
{
  _enrich_disp_var.resize(_nl_vnames.size());
  for (unsigned int i = 0; i < _nl_vnames.size(); ++i)
    _enrich_disp_var[i] = coupled("enrichment_displacement_var", i);

  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in EnrichStressDivergenceTensors");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in EnrichStressDivergenceTensors");
}

void
EnrichStressDivergenceTensors::computeResidual()
{
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

void
EnrichStressDivergenceTensors::computeJacobian()
{
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

Real
EnrichStressDivergenceTensors::computeQpResidual()
{
  // calculate the near-tip enrichement function
  std::vector<Real> B, Bx, By, Br, Bt, Bxl, Byl; 
  // Br  : derivative w.r.t r  
  // Bt  : derivative w.r.t theta
  // Bxl : derivative w.r.t local x (crack tip coordinate)
  // Byl : derivative w.r.t local y (crack tip coordinate)
  // Bx  : derivative w.r.t global x
  // By  : derivative w.r.t global y
  B.resize(4);
  Bx.resize(4);
  By.resize(4);
  Br.resize(4);
  Bt.resize(4);
  Bxl.resize(4);
  Byl.resize(4);

  std::vector<std::vector<Real> > BI;
  BI.resize(4);
  for (unsigned int i = 0; i < BI.size(); ++i)
  {
    BI[i].resize(4);

    Point crack_tip(0.5, 1.0, 0); //crack tip is at (0.5, 0.5, 0)
    Node * node_i = _current_elem->get_node(i);

    Real x_to_tip = (*node_i)(0) - crack_tip(0);
    Real y_to_tip = (*node_i)(1) - crack_tip(1);

    Real alpha = 0.0; // crack direction

    Real x_local =  std::cos(alpha) * x_to_tip + std::sin(alpha) * y_to_tip;
    Real y_local = -std::sin(alpha) * x_to_tip + std::cos(alpha) * y_to_tip;

    Real r = std::sqrt(x_local * x_local + y_local * y_local);

    if (r < 0.001)
      r = 0.001;

    Real theta = std::atan2(y_local, x_local);

    BI[i][0] = std::sqrt(r) * std::sin(theta / 2.0);
    BI[i][1] = std::sqrt(r) * std::cos(theta / 2.0);
    BI[i][2] = std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
    BI[i][3] = std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);
  }


  Point crack_tip(0.5, 1.0, 0); //crack tip is at (0.5, 0.5, 0)
  Point q_pt = _q_point[_qp];

  Real x_to_tip = q_pt(0) - crack_tip(0);
  Real y_to_tip = q_pt(1) - crack_tip(1);

  Real alpha = 0.0; // crack direction
  
  Real x_local =  std::cos(alpha) * x_to_tip + std::sin(alpha) * y_to_tip;
  Real y_local = -std::sin(alpha) * x_to_tip + std::cos(alpha) * y_to_tip;
 
  Real r = std::sqrt(x_local * x_local + y_local * y_local);

  if (r < 0.001)
    r = 0.001;

  Real theta = std::atan2(y_local, x_local);

  B[0] = std::sqrt(r) * std::sin(theta / 2.0);
  B[1] = std::sqrt(r) * std::cos(theta / 2.0);
  B[2] = std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
  B[3] = std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);


  Br[0] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0);
  Bt[0] = std::sqrt(r) / 2.0 * std::cos(theta / 2.0);
  Br[1] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
  Bt[1] = -std::sqrt(r) / 2.0 * std::sin(theta / 2.0);
  Br[2] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
  Bt[2] = std::sqrt(r) * (0.5 * std::cos(theta/2.0) * std::sin(theta) + std::sin(theta / 2.0) * std::cos(theta));
  Br[2] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);
  Bt[2] = std::sqrt(r) * (-0.5 * std::sin(theta/2.0) * std::sin(theta) + std::cos(theta / 2.0) * std::cos(theta));

  //Real r_xl = std::cos(theta);
  //Real r_yl = std::sin(theta);
  //Real t_xl = -std::sin(theta) / r;
  //Real t_yl = std::cos(theta) / r;

  Bxl[0] = -0.5 / std::sqrt(r) * std::sin(theta / 2.0);
  Byl[0] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
  Bxl[1] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
  Byl[1] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0);
  Bxl[2] = -0.5 / std::sqrt(r) * std::sin(1.5 * theta) * std::sin(theta);
  Byl[2] = 0.5 / std::sqrt(r) * (std::sin(theta / 2.0) + std::sin(1.5 * theta) * std::cos(theta));
  Bxl[3] = -0.5 / std::sqrt(r) * std::cos(1.5 * theta) * std::sin(theta);
  Byl[3] = 0.5 / std::sqrt(r) * (std::cos(theta / 2.0) + std::cos(1.5 * theta) * std::cos(theta));

  for (unsigned int i = 0; i < 4; ++i)
  {
    Bx[i] = Bxl[i] * std::cos(alpha) - Byl[i] * std::sin(alpha);
    By[i] = Bxl[i] * std::sin(alpha) + Byl[i] * std::cos(alpha);
    //std::cout << "Bx[" << i << "] = " << Bx[i] << std::endl;
    //std::cout << "By[" << i << "] = " << By[i] << std::endl;
  }

  RealVectorValue grad_B(Bx[_enrichment_component], By[_enrichment_component], 0.0);

  return _stress[_qp].row(_component) * (_grad_test[_i][_qp] * (B[_enrichment_component] - BI[_i][_enrichment_component]) + _test[_i][_qp] * grad_B);

}

Real
EnrichStressDivergenceTensors::computeQpJacobian()
{
  // calculate the near-tip enrichement function
  std::vector<Real> B, Bx, By, Br, Bt, Bxl, Byl; 
  // Br  : derivative w.r.t r  
  // Bt  : derivative w.r.t theta
  // Bxl : derivative w.r.t local x (crack tip coordinate)
  // Byl : derivative w.r.t local y (crack tip coordinate)
  // Bx  : derivative w.r.t global x
  // By  : derivative w.r.t global y
  B.resize(4);
  Bx.resize(4);
  By.resize(4);
  Br.resize(4);
  Bt.resize(4);
  Bxl.resize(4);
  Byl.resize(4);

  std::vector<std::vector<Real> > BI;
  BI.resize(4);
  for (unsigned int i = 0; i < BI.size(); ++i)
  {
    BI[i].resize(4);

    Point crack_tip(0.5, 1.0, 0); //crack tip is at (0.5, 0.5, 0)
    Node * node_i = _current_elem->get_node(i);

    Real x_to_tip = (*node_i)(0) - crack_tip(0);
    Real y_to_tip = (*node_i)(1) - crack_tip(1);

    Real alpha = 0.0; // crack direction

    Real x_local =  std::cos(alpha) * x_to_tip + std::sin(alpha) * y_to_tip;
    Real y_local = -std::sin(alpha) * x_to_tip + std::cos(alpha) * y_to_tip;

    Real r = std::sqrt(x_local * x_local + y_local * y_local);

    if (r < 0.001)
      r = 0.001;

    Real theta = std::atan2(y_local, x_local);

    BI[i][0] = std::sqrt(r) * std::sin(theta / 2.0);
    BI[i][1] = std::sqrt(r) * std::cos(theta / 2.0);
    BI[i][2] = std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
    BI[i][3] = std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);
  }


  Point crack_tip(0.5, 1.0, 0); //crack tip is at (0.5, 0.5, 0)
  Point q_pt = _q_point[_qp];

  Real x_to_tip = q_pt(0) - crack_tip(0);
  Real y_to_tip = q_pt(1) - crack_tip(1);

  Real alpha = 0.0; // crack direction
  
  Real x_local =  std::cos(alpha) * x_to_tip + std::sin(alpha) * y_to_tip;
  Real y_local = -std::sin(alpha) * x_to_tip + std::cos(alpha) * y_to_tip;
 
  Real r = std::sqrt(x_local * x_local + y_local * y_local);

  if (r < 0.001)
    r = 0.001;

  Real theta = std::atan2(y_local, x_local);

  B[0] = std::sqrt(r) * std::sin(theta / 2.0);
  B[1] = std::sqrt(r) * std::cos(theta / 2.0);
  B[2] = std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
  B[3] = std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);


  Br[0] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0);
  Bt[0] = std::sqrt(r) / 2.0 * std::cos(theta / 2.0);
  Br[1] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
  Bt[1] = -std::sqrt(r) / 2.0 * std::sin(theta / 2.0);
  Br[2] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
  Bt[2] = std::sqrt(r) * (0.5 * std::cos(theta/2.0) * std::sin(theta) + std::sin(theta / 2.0) * std::cos(theta));
  Br[2] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);
  Bt[2] = std::sqrt(r) * (-0.5 * std::sin(theta/2.0) * std::sin(theta) + std::cos(theta / 2.0) * std::cos(theta));

  //Real r_xl = std::cos(theta);
  //Real r_yl = std::sin(theta);
  //Real t_xl = -std::sin(theta) / r;
  //Real t_yl = std::cos(theta) / r;

  Bxl[0] = -0.5 / std::sqrt(r) * std::sin(theta / 2.0);
  Byl[0] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
  Bxl[1] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
  Byl[1] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0);
  Bxl[2] = -0.5 / std::sqrt(r) * std::sin(1.5 * theta) * std::sin(theta);
  Byl[2] = 0.5 / std::sqrt(r) * (std::sin(theta / 2.0) + std::sin(1.5 * theta) * std::cos(theta));
  Bxl[3] = -0.5 / std::sqrt(r) * std::cos(1.5 * theta) * std::sin(theta);
  Byl[3] = 0.5 / std::sqrt(r) * (std::cos(theta / 2.0) + std::cos(1.5 * theta) * std::cos(theta));

  for (unsigned int i = 0; i < 4; ++i)
  {
    Bx[i] = Bxl[i] * std::cos(alpha) - Byl[i] * std::sin(alpha);
    By[i] = Bxl[i] * std::sin(alpha) + Byl[i] * std::cos(alpha);
  }

  RealVectorValue grad_B(Bx[_enrichment_component], By[_enrichment_component], 0.0);
  RealVectorValue grad_test = _grad_test[_i][_qp] * (B[_enrichment_component] - BI[_i][_enrichment_component]) + _test[_i][_qp] * grad_B;
  RealVectorValue grad_phi = _grad_phi[_j][_qp] * (B[_enrichment_component] - BI[_j][_enrichment_component]) + _phi[_j][_qp] * grad_B;

  return ElasticityTensorTools::elasticJacobian(_Jacobian_mult[_qp], _component, _component, grad_test, grad_phi);
}

Real
EnrichStressDivergenceTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component = 0;
  unsigned int coupled_enrichment_component = 0;
  bool active(false);

  for (unsigned int i = 0; i < _enrich_disp_var.size(); ++i)
  {
    if (jvar == _enrich_disp_var[i])
    {
      coupled_component = i % 2;
      coupled_enrichment_component = i / 2;
      active = true;
    }
  }

  if (active)
  {
    // calculate the near-tip enrichement function
    std::vector<Real> B, Bx, By, Br, Bt, Bxl, Byl; 
    // Br  : derivative w.r.t r  
    // Bt  : derivative w.r.t theta
    // Bxl : derivative w.r.t local x (crack tip coordinate)
    // Byl : derivative w.r.t local y (crack tip coordinate)
    // Bx  : derivative w.r.t global x
    // By  : derivative w.r.t global y
    B.resize(4);
    Bx.resize(4);
    By.resize(4);
    Br.resize(4);
    Bt.resize(4);
    Bxl.resize(4);
    Byl.resize(4);

    std::vector<std::vector<Real> > BI;
    BI.resize(4);
    for (unsigned int i = 0; i < BI.size(); ++i)
    {
      BI[i].resize(4);

      Point crack_tip(0.5, 1.0, 0); //crack tip is at (0.5, 0.5, 0)
      Node * node_i = _current_elem->get_node(i);

      Real x_to_tip = (*node_i)(0) - crack_tip(0);
      Real y_to_tip = (*node_i)(1) - crack_tip(1);

      Real alpha = 0.0; // crack direction

      Real x_local =  std::cos(alpha) * x_to_tip + std::sin(alpha) * y_to_tip;
      Real y_local = -std::sin(alpha) * x_to_tip + std::cos(alpha) * y_to_tip;

      Real r = std::sqrt(x_local * x_local + y_local * y_local);

      if (r < 0.001)
        r = 0.001;

      Real theta = std::atan2(y_local, x_local);

      BI[i][0] = std::sqrt(r) * std::sin(theta / 2.0);
      BI[i][1] = std::sqrt(r) * std::cos(theta / 2.0);
      BI[i][2] = std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
      BI[i][3] = std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);
    }


    Point crack_tip(0.5, 1.0, 0); //crack tip is at (0.5, 0.5, 0)
    Point q_pt = _q_point[_qp];

    Real x_to_tip = q_pt(0) - crack_tip(0);
    Real y_to_tip = q_pt(1) - crack_tip(1);

    Real alpha = 0.0; // crack direction

    Real x_local =  std::cos(alpha) * x_to_tip + std::sin(alpha) * y_to_tip;
    Real y_local = -std::sin(alpha) * x_to_tip + std::cos(alpha) * y_to_tip;

    Real r = std::sqrt(x_local * x_local + y_local * y_local);

    if (r < 0.001)
      r = 0.001;

    Real theta = std::atan2(y_local, x_local);

    B[0] = std::sqrt(r) * std::sin(theta / 2.0);
    B[1] = std::sqrt(r) * std::cos(theta / 2.0);
    B[2] = std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
    B[3] = std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);


    Br[0] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0);
    Bt[0] = std::sqrt(r) / 2.0 * std::cos(theta / 2.0);
    Br[1] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
    Bt[1] = -std::sqrt(r) / 2.0 * std::sin(theta / 2.0);
    Br[2] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0) * std::sin(theta);
    Bt[2] = std::sqrt(r) * (0.5 * std::cos(theta/2.0) * std::sin(theta) + std::sin(theta / 2.0) * std::cos(theta));
    Br[2] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);
    Bt[2] = std::sqrt(r) * (-0.5 * std::sin(theta/2.0) * std::sin(theta) + std::cos(theta / 2.0) * std::cos(theta));

    //Real r_xl = std::cos(theta);
    //Real r_yl = std::sin(theta);
    //Real t_xl = -std::sin(theta) / r;
    //Real t_yl = std::cos(theta) / r;

    Bxl[0] = -0.5 / std::sqrt(r) * std::sin(theta / 2.0);
    Byl[0] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
    Bxl[1] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0);
    Byl[1] = 0.5 / std::sqrt(r) * std::sin(theta / 2.0);
    Bxl[2] = -0.5 / std::sqrt(r) * std::sin(1.5 * theta) * std::sin(theta);
    Byl[2] = 0.5 / std::sqrt(r) * (std::sin(theta / 2.0) + std::sin(1.5 * theta) * std::cos(theta));
    Bxl[3] = -0.5 / std::sqrt(r) * std::cos(1.5 * theta) * std::sin(theta);
    Byl[3] = 0.5 / std::sqrt(r) * (std::cos(theta / 2.0) + std::cos(1.5 * theta) * std::cos(theta));

    for (unsigned int i = 0; i < 4; ++i)
    {
      Bx[i] = Bxl[i] * std::cos(alpha) - Byl[i] * std::sin(alpha);
      By[i] = Bxl[i] * std::sin(alpha) + Byl[i] * std::cos(alpha);
    }

    RealVectorValue grad_B_test(Bx[_enrichment_component], By[_enrichment_component], 0.0);
    RealVectorValue grad_B_phi(Bx[coupled_enrichment_component], By[coupled_enrichment_component], 0.0);
    RealVectorValue grad_test = _grad_test[_i][_qp] * (B[_enrichment_component] - BI[_i][_enrichment_component]) + _test[_i][_qp] * grad_B_test;
    RealVectorValue grad_phi = _grad_phi[_j][_qp] * (B[coupled_enrichment_component] - BI[_j][coupled_enrichment_component]) + _phi[_j][_qp] * grad_B_phi;

    return ElasticityTensorTools::elasticJacobian(_Jacobian_mult[_qp], _component, coupled_component,
                                                  grad_test, grad_phi);

  }

    return 0;
}

