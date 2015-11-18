/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This userobject calculates the configurational force
//
#include "XFEMMeanStress.h"
#include "libmesh/fe_interface.h"
#include "DisplacedProblem.h"

template<>
InputParameters validParams<XFEMMeanStress>()
{
  InputParameters params = validParams<ElementUserObject>();
  params += validParams<MaterialTensorCalculator>();
  params.addRequiredParam<std::string>("tensor", "The material tensor name.");
  params.addParam<Real>("radius", "Radius for average out stress");
  return params;
}

XFEMMeanStress::XFEMMeanStress(const InputParameters & parameters):
    ElementUserObject(parameters),
    _material_tensor_calculator(parameters),
    _tensor(getMaterialProperty<SymmTensor>(getParam<std::string>("tensor")))
{

  _fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (_fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMarkerUserObject");
  _xfem = _fe_problem->get_xfem();

  if (isParamValid("radius"))
  {
    _radius = getParam<Real>("radius");
  }
  else
    mooseError("MeanStress error: must set radius.");
}

void
XFEMMeanStress::initialize()
{ 
  _crack_front_points.clear();
  _elem_id_crack_tip.clear();
  _stress_tensor.clear();

  _num_crack_front_points = _xfem->num_crack_tips();

  _stress_tensor.resize(_num_crack_front_points*9);
  
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] = 0.0;

  _xfem->get_crack_tip_origin(_elem_id_crack_tip, _crack_front_points);

}

std::vector<Real>
XFEMMeanStress::getStressTensor()
{
  std::vector<Real> StressTensor(_num_crack_front_points*9);
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    StressTensor[i] = 0.0;

  unsigned int numqp = _qrule->n_points();

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Point crack_front = _crack_front_points[i];

    const Elem * undisplaced_elem  = NULL;
    if(_fe_problem->getDisplacedProblem() != NULL)
      undisplaced_elem = _fe_problem->getDisplacedProblem()->refMesh().elem(_current_elem->id());
    else
      undisplaced_elem = _current_elem;
    
    for ( unsigned int qp = 0; qp < numqp; ++qp )
    {
      Real flag = _xfem->flag_qp_inside(undisplaced_elem, _q_point[qp]); //qp inside (flag = 1) or ouside (flag = 0) real domain
      Point dist_to_crack_front_vector = _q_point[qp] - crack_front;
      Real dist = std::pow(dist_to_crack_front_vector.size_sq(),0.5);
      Real fact = (1.0-dist/_radius);
      if (dist < _radius && flag > 0.5)
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
      }
    }
  }
  return StressTensor;
}

void 
XFEMMeanStress::execute()
{
  std::vector<Real> StressTensor = getStressTensor();
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] += StressTensor[i];
}

void
XFEMMeanStress::threadJoin(const UserObject & y)
{
  const XFEMMeanStress & pps = static_cast<const XFEMMeanStress &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] += pps._stress_tensor[i];
}

void 
XFEMMeanStress::finalize()
{
  _xfem->clear_crack_propagation_direction();
  gatherSum(_stress_tensor);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    RealVectorValue direction;
    SymmTensor average_tensor;
    average_tensor(0,0) = _stress_tensor[i*9+0];
    average_tensor(0,1) = _stress_tensor[i*9+1];
    average_tensor(0,2) = _stress_tensor[i*9+2];
    average_tensor(1,0) = _stress_tensor[i*9+3];
    average_tensor(1,1) = _stress_tensor[i*9+4];
    average_tensor(1,2) = _stress_tensor[i*9+5];
    average_tensor(2,0) = _stress_tensor[i*9+6];
    average_tensor(2,1) = _stress_tensor[i*9+7];
    average_tensor(2,2) = _stress_tensor[i*9+8];
    Real tensor_quantity = _material_tensor_calculator.getTensorQuantity(average_tensor,&_q_point[0],direction);
    direction /= pow(direction.size_sq(),0.5);
    Point normal(0.0,0.0,0.0);
    normal(0) = direction(1);
    normal(1) = -direction(0);
    _xfem->update_crack_propagation_direction(_elem_id_crack_tip[i], normal);
    //std::cout << "crack front index (" << i << ") : average stress direction  = " << normal << std::endl; 
  } 
}

