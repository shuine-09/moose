/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEM.h"
#include "XFEMCrackTipMeanStress.h"
#include "libmesh/fe_interface.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<XFEMCrackTipMeanStress>()
{
  InputParameters params = validParams<ElementUserObject>();
  params += validParams<MaterialTensorCalculator>();
  params.addRequiredParam<std::string>("tensor", "The material tensor name.");
  params.addParam<Real>("radius", "Radius of gaussian weight function to average out stress");
  params.addParam<Real>("critical_stress", 0.0, "Critical stress.");
  return params;
}

XFEMCrackTipMeanStress::XFEMCrackTipMeanStress(const InputParameters & parameters):
  ElementUserObject(parameters),
  _material_tensor_calculator(parameters),
  _tensor(getMaterialProperty<SymmTensor>(getParam<std::string>("tensor"))),
  _critical_stress(getParam<Real>("critical_stress"))
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMCrackTipMeanStress");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMCrackTipMeanStress");
  if (isNodal())
    mooseError("XFEMCrackTipMeanStress can only be run on an element variable");

  if (isParamValid("radius"))
  {
    _radius = getParam<Real>("radius");
  }
  else
    mooseError("XFEMCrackTipMeanStress error: must set radius.");
}

void
XFEMCrackTipMeanStress::initialize()
{ 
  _crack_front_points.clear();
  _stress_tensor.clear();

  _num_crack_front_points = _xfem->numberCrackTips();
  _stress_tensor.resize(_num_crack_front_points*9);

  _weights.clear();
  _weights.resize(_num_crack_front_points);

  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] = 0.0;

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
    _weights[i] = 0;

  _xfem->getCrackTipOrigin(_elem_id_crack_tip, _crack_front_points);

}

std::vector<Real>
XFEMCrackTipMeanStress::getStressTensor()
{
  std::vector<Real> StressTensor(_num_crack_front_points*9);
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    StressTensor[i] = 0.0;

  unsigned int numqp = _qrule->n_points();

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Point crack_front = _crack_front_points[i];

    const Elem * undisplaced_elem  = NULL;
    FEProblem * _fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
    if(_fe_problem->getDisplacedProblem() != NULL)
      undisplaced_elem = _fe_problem->getDisplacedProblem()->refMesh().elemPtr(_current_elem->id());
    else
      undisplaced_elem = _current_elem;

    for ( unsigned int qp = 0; qp < numqp; ++qp )
    {
      Real flag = _xfem->flagQpointInside(undisplaced_elem, _q_point[qp]); //qp inside (flag = 1) or ouside (flag = 0) real domain
      Point dist_to_crack_front_vector = _q_point[qp] - crack_front;
      Real dist = std::pow(dist_to_crack_front_vector.size_sq(),0.5);
      Real fact = 1.0/(pow(2*libMesh::pi, 1.5) * pow(_radius, 3.0)) * std::exp(-0.5 * pow(dist/_radius,2.0));
      if (dist < _radius * 2 && flag > 0.5)
      {
        StressTensor[i*9+0] += _tensor[qp](0,0) * fact;
        StressTensor[i*9+1] += _tensor[qp](0,1) * fact;
        StressTensor[i*9+2] += _tensor[qp](0,2) * fact;
        StressTensor[i*9+3] += _tensor[qp](1,0) * fact;
        StressTensor[i*9+4] += _tensor[qp](1,1) * fact;
        StressTensor[i*9+5] += _tensor[qp](1,2) * fact;
        StressTensor[i*9+6] += _tensor[qp](2,0) * fact;
        StressTensor[i*9+7] += _tensor[qp](2,1) * fact;
        StressTensor[i*9+8] += _tensor[qp](2,2) * fact;
        _weights[i] += fact;
      }
    }
  }
  return StressTensor;
}

void 
XFEMCrackTipMeanStress::execute()
{
  std::vector<Real> StressTensor = getStressTensor();
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] += StressTensor[i];
}

void
XFEMCrackTipMeanStress::threadJoin(const UserObject & y)
{
  const XFEMCrackTipMeanStress & pps = static_cast<const XFEMCrackTipMeanStress &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] += pps._stress_tensor[i];
}

void 
XFEMCrackTipMeanStress::finalize()
{
  _xfem->clearCrackPropagationDirection();
  _xfem->clearDoesCrackPropagate();

  gatherSum(_stress_tensor);
  gatherSum(_weights);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    RealVectorValue direction;
    SymmTensor average_tensor;
    average_tensor(0,0) = _stress_tensor[i*9+0] / _weights[i];
    average_tensor(0,1) = _stress_tensor[i*9+1] / _weights[i];
    average_tensor(0,2) = _stress_tensor[i*9+2] / _weights[i];
    average_tensor(1,0) = _stress_tensor[i*9+3] / _weights[i];
    average_tensor(1,1) = _stress_tensor[i*9+4] / _weights[i];
    average_tensor(1,2) = _stress_tensor[i*9+5] / _weights[i];
    average_tensor(2,0) = _stress_tensor[i*9+6] / _weights[i];
    average_tensor(2,1) = _stress_tensor[i*9+7] / _weights[i];
    average_tensor(2,2) = _stress_tensor[i*9+8] / _weights[i];
    Real tensor_quantity = _material_tensor_calculator.getTensorQuantity(average_tensor,_q_point[0],direction);
    direction /= pow(direction.size_sq(),0.5);
    Point normal(0.0,0.0,0.0);
    normal(0) = -direction(1);
    normal(1) = direction(0);
    bool does_crack_propagate = (tensor_quantity > _critical_stress);
    
    _xfem->updateDoesCrackPropagate(_elem_id_crack_tip[i], does_crack_propagate);
    _xfem->updateCrackPropagationDirection(_elem_id_crack_tip[i], normal);

    std::cout << "stress = " << tensor_quantity << " crack propagate ? = " << does_crack_propagate << std::endl;
    std::cout << "crack front index (" << i << ") " <<  _crack_front_points[i] << ", average stress direction  = " << normal << std::endl; 
  } 
}
