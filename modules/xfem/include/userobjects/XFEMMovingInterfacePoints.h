/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMMOVINGINTERFACEPOINTS_H
#define XFEMMOVINGINTERFACEPOINTS_H

#include "GeneralUserObject.h"

class XFEM;

class XFEMMovingInterfacePoints;

template<>
InputParameters validParams<XFEMMovingInterfacePoints>();

class XFEMMovingInterfacePoints : public GeneralUserObject
{
 public:
  XFEMMovingInterfacePoints(const InputParameters & parameters);

  virtual void finalize();
  
  virtual void execute(){};

  virtual void initialize(){};

  protected:

  MooseSharedPointer<XFEM> _xfem;

  bool skip;
  unsigned int _save_time_step;
};

#endif // XFEMMOVINGINTERFACEPOINTS_H
