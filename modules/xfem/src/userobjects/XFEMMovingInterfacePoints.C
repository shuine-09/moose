/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "XFEMMovingInterfacePoints.h"
#include "XFEM.h"

template<>
InputParameters validParams<XFEMMovingInterfacePoints>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addClassDescription("XFEM moving interface points class. Override the virtual functions in your class");
  return params;
}

XFEMMovingInterfacePoints::XFEMMovingInterfacePoints(const InputParameters & parameters) :
    GeneralUserObject(parameters)
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMovingInterfacePoints");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMMovingInterfacePoints");

  _save_time_step = 0;
  skip = false;
}

void
XFEMMovingInterfacePoints::finalize()
{
  if (_save_time_step != _t_step)
    skip = false;
  else
    skip = true;

  if (!skip)
  {
    for (unsigned int i = 0; i < _xfem->numberInterfacePoints(); ++i)
    {
      Real interface_quantity = _xfem->getInterfaceQuantity(i);
      std::cout << "WJ: interface_quantiies[" << i << "] = " << interface_quantity << std::endl;

      Point interface_point = _xfem->getInterfacePoint(i);

      interface_point(0) += interface_quantity * 0.05;

      //if (i > 3 && i < 7)
      //  interface_point(0) += interface_quantity * 0.0025;

      _xfem->setInterfacePoint(i, interface_point);
    }

    _save_time_step = _t_step;
  }
}
