/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef COMPUTEELEMAUXVARSTHREAD_H
#define COMPUTEELEMAUXVARSTHREAD_H

// libMesh includes
#include "libmesh/elem_range.h"

// MOOSE includes
#include "ThreadedElementLoop.h"
#include "MooseObjectWarehouse.h"
#include "AuxKernel.h"

// Forward declarations
class FEProblem;
class AuxiliarySystem;

class ComputeElemAuxVarsThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  ComputeElemAuxVarsThread(FEProblem & problem, AuxiliarySystem & sys, const MooseObjectWarehouse<AuxKernel> & storage, bool need_materials);
  // Splitting Constructor
  ComputeElemAuxVarsThread(ComputeElemAuxVarsThread & x, Threads::split split);

  virtual ~ComputeElemAuxVarsThread();

  virtual void subdomainChanged();
  virtual void onElement(const Elem *elem);
  virtual void post();

  void join(const ComputeElemAuxVarsThread & /*y*/);

protected:
  AuxiliarySystem & _aux_sys;

  /// Storage object containing active AuxKernel objects
  const MooseObjectWarehouse<AuxKernel> & _storage;

  bool _need_materials;
};

#endif //COMPUTEELEMAUXVARSTHREAD_H
