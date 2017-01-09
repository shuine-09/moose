/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeEnrichStrain.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<ComputeEnrichStrain>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<std::vector<NonlinearVariableName> >("enrichment_displacement", "The enrichment displacement");
  params.addRequiredCoupledVar("displacements", "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeEnrichStrain::ComputeEnrichStrain(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),
    _enrich_disp(3),
    _grad_enrich_disp(3),
    _enrich_variable(4),
    _enrich_strain(declareProperty<RankTwoTensor>("enrich_strain")),
    _phi(_assembly.phi()),
    _grad_phi(_assembly.gradPhi()),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
    _mechanical_strain(declareProperty<RankTwoTensor>(_base_name + "mechanical_strain")),
    _total_strain(declareProperty<RankTwoTensor>(_base_name + "total_strain")),
    _eigenstrain(getDefaultMaterialProperty<RankTwoTensor>(_base_name + "stress_free_strain"))
{
  for (unsigned int i = 0; i < _enrich_variable.size(); ++i)
    _enrich_variable[i].resize(2); //TODO 3D

  const std::vector<NonlinearVariableName> & nl_vnames = getParam<std::vector<NonlinearVariableName> >("enrichment_displacement");
  NonlinearSystem & nl = _fe_problem.getNonlinearSystem();

  for (unsigned int i = 0; i < 4; ++i) //TODO : total 4 enrichment functions per node along one direction
  {
    _enrich_variable[i][0] = &(nl.getVariable(0, nl_vnames[i * 2]));
    _enrich_variable[i][1] = &(nl.getVariable(0, nl_vnames[i * 2 + 1]));
  }
  
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of variables supplied in 'displacements' must match the mesh dimension.");

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
}

void
ComputeEnrichStrain::computeQpProperties()
{
  FEType fe_type(Utility::string_to_enum<Order>("first"),Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule (_qrule);

  // The values of the shape functions at the quadrature points
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  fe->reinit (_current_elem); 

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
  Bt[2] = std::sqrt(r) * (0.5 * std::cos(theta / 2.0) * std::sin(theta) + std::sin(theta / 2.0) * std::cos(theta));
  Br[2] = 0.5 / std::sqrt(r) * std::cos(theta / 2.0) * std::sin(theta);
  Bt[2] = std::sqrt(r) * (-0.5 * std::sin(theta / 2.0) * std::sin(theta) + std::cos(theta / 2.0) * std::cos(theta));

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

  NonlinearSystem & nl = _fe_problem.getNonlinearSystem();
  const NumericVector<Number> & sln = *nl.currentSolution();

  //if(_current_elem->contains_point(crack_tip))
  if (0)
  {
    std::cout << "q point = " << q_pt << std::endl;
    std::vector<Point> qpoint;
    qpoint.push_back(q_pt);
    std::vector<Point> ref_point;
    FEInterface::inverse_map(2, fe_type, _current_elem, qpoint, ref_point);
    std::cout << "ref q point = " << ref_point[0] << std::endl;
    for (unsigned int i = 0; i < 4; ++i)
    {
      std::cout << "Bx[" << i << "] = " << Bx[i] << std::endl;
      std::cout << "By[" << i << "] = " << By[i] << std::endl;
      std::cout << "phi[" << i << "] = " << phi[i][_qp] << std::endl;
      std::cout << "dphi[" << i << "] = " << dphi[i][_qp] << std::endl;
    }
  }

  for (unsigned int m = 0; m < 2; ++m) //TODO: 3D
  {
    _enrich_disp[m] = 0.0;
    _grad_enrich_disp[m].zero();
    for (unsigned int i = 0; i < _current_elem->n_nodes(); ++i)
    {
      Node * node_i = _current_elem->get_node(i);
      dof_id_type dofx = node_i->dof_number(nl.number(), 0, 0);
      dof_id_type dofy = node_i->dof_number(nl.number(), 1, 0);
      Real soln_localx = sln(dofx);
      Real soln_localy = sln(dofy);
      if(_current_elem->contains_point(crack_tip))
        std::cout << " node(" << i << ") x = " << soln_localx << ", y = " << soln_localy << std::endl;

      for (unsigned int j = 0; j < 4; ++j)
      {
        dof_id_type dof = node_i->dof_number(nl.number(), _enrich_variable[j][m]->number(), 0);
        Real soln_local = sln(dof);
        if(_current_elem->contains_point(crack_tip))
          std::cout << "j = " << j << ", m = " << m << ", _enrich_variable[j][m]->number() = " << _enrich_variable[j][m]->number() << ", dof = " << dof << ", sln = " << soln_local << std::endl;
        _enrich_disp[m] += phi[i][_qp] * (B[j] - BI[i][j])* soln_local;
        RealVectorValue grad_B(Bx[j], By[j], 0.0);
        _grad_enrich_disp[m] += (dphi[i][_qp] * (B[j] - BI[i][j]) + phi[i][_qp] * grad_B) * soln_local;
      }
    }
  }


  _enrich_disp[2] = 0.0; //TODO: 3D
  _grad_enrich_disp[2].zero();

  //strain = (grad_disp + grad_disp^T)/2
  RankTwoTensor grad_tensor_enrich(_grad_enrich_disp[0], _grad_enrich_disp[1], _grad_enrich_disp[2]);

  _enrich_strain[_qp] = (grad_tensor_enrich + grad_tensor_enrich.transpose()) / 2.0;

  /*
  Point tip(0.5, 0.5, 0);
  if(_current_elem->contains_point(tip))
  {
    std::cout << "_qp = " << _qp << ", enrich_strain = " << std::endl;
    _enrich_strain[_qp].print();
  }
  */

  RealVectorValue grad_disp_z(0.0, 0.0, 0.0);

  RankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], grad_disp_z); //TODO 3D!!

  _total_strain[_qp] = (grad_tensor + grad_tensor.transpose()) / 2.0;

  _total_strain[_qp] += _enrich_strain[_qp];

  _mechanical_strain[_qp] = _total_strain[_qp];

  Point tip(0.5, 1.0, 0);
  if(_current_elem->contains_point(tip))
  //if(0)
  {
    //for(unsigned int i = 0; i < 4; i++)
    //  for(unsigned int j = 0; j < 4; j++)
    //    std::cout << "BI[" << i << "][" << j << "] = " << BI[i][j] << std::endl;


    std::cout << "r = " <<  r << ", theta = " << theta << std::endl;
    std::cout << "_qp = " << _qp << ", enrich_strain = " << std::endl;
    _enrich_strain[_qp].print();
    std::cout << "strain = " << std::endl;
    ((grad_tensor + grad_tensor.transpose()) / 2.0).print();
    std::cout << "total strain = " << std::endl;
    _total_strain[_qp].print();
  }


  //Remove the Eigen strain
  _mechanical_strain[_qp] -= _eigenstrain[_qp];
}

void
ComputeEnrichStrain::initQpStatefulProperties()
{
  _enrich_strain[_qp].zero();
  _mechanical_strain[_qp].zero();
  _total_strain[_qp].zero();
}
