/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "XFEMMeshCutUserObject.h"
#include "XFEM.h"

#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"

template<>
InputParameters validParams<XFEMMeshCutUserObject>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("XFEM mesh cut base class. Override the virtual functions in your class");
  return params;
}

XFEMMeshCutUserObject::XFEMMeshCutUserObject(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _mesh(_subproblem.mesh())
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMeshCutUserObject");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMMeshCutUserObject");
  if (isNodal())
    mooseError("XFEMMeshCutUserObject can only be run on an element variable");
}

void
XFEMMeshCutUserObject::initialize()
{
  _marked_elems_distance.clear();
  _marked_elems_host_id.clear();
}

void
XFEMMeshCutUserObject::execute()
{
}

void
XFEMMeshCutUserObject::threadJoin(const UserObject &y)
{
  const XFEMMeshCutUserObject &xmcuo = dynamic_cast<const XFEMMeshCutUserObject &>(y);

  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_distance.begin();
        mit != xmcuo._marked_elems_distance.end();
        ++mit )
  {
    _marked_elems_distance[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_host_id.begin();
        mit != xmcuo._marked_elems_host_id.end();
        ++mit )
  {
    _marked_elems_host_id[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

}

void
XFEMMeshCutUserObject::finalize()
{
  _communicator.set_union(_marked_elems_distance);
  _communicator.set_union(_marked_elems_host_id);
  _xfem->clearUOMarkedElems();

  std::map<unsigned int, RealVectorValue >::iterator mit;
  for (mit = _marked_elems_host_id.begin(); mit != _marked_elems_host_id.end(); ++mit)
  {
   _xfem->addUOMarkedElemHostID(mit->first, mit->second);
  }

  for (mit = _marked_elems_distance.begin(); mit != _marked_elems_distance.end(); ++mit)
  {
   _xfem->addUOMarkedElemDistance(mit->first, mit->second);
  }

  _marked_elems_distance.clear();
  _marked_elems_host_id.clear();
}
