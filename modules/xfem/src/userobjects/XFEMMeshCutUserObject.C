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
  _marked_elems_distance2.clear();
  _marked_elems_host_id2.clear();
  _marked_elems_distance3.clear();
  _marked_elems_host_id3.clear();
  _marked_elems_distance4.clear();
  _marked_elems_host_id4.clear();
  _marked_elems_distance5.clear();
  _marked_elems_host_id5.clear();
  _marked_elems_distance6.clear();
  _marked_elems_host_id6.clear();
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

  // 2
  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_distance2.begin();
        mit != xmcuo._marked_elems_distance2.end();
        ++mit )
  {
    _marked_elems_distance2[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_host_id2.begin();
        mit != xmcuo._marked_elems_host_id2.end();
        ++mit )
  {
    _marked_elems_host_id2[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  //3
  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_distance3.begin();
        mit != xmcuo._marked_elems_distance3.end();
        ++mit )
  {
    _marked_elems_distance3[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_host_id3.begin();
        mit != xmcuo._marked_elems_host_id3.end();
        ++mit )
  {
    _marked_elems_host_id3[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  //4
  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_distance4.begin();
        mit != xmcuo._marked_elems_distance4.end();
        ++mit )
  {
    _marked_elems_distance4[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_host_id4.begin();
        mit != xmcuo._marked_elems_host_id4.end();
        ++mit )
  {
    _marked_elems_host_id4[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  //5
  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_distance5.begin();
        mit != xmcuo._marked_elems_distance5.end();
        ++mit )
  {
    _marked_elems_distance5[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_host_id5.begin();
        mit != xmcuo._marked_elems_host_id5.end();
        ++mit )
  {
    _marked_elems_host_id5[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  //6
  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_distance6.begin();
        mit != xmcuo._marked_elems_distance6.end();
        ++mit )
  {
    _marked_elems_distance6[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }

  for ( std::map<unsigned int, RealVectorValue >::const_iterator mit = xmcuo._marked_elems_host_id6.begin();
        mit != xmcuo._marked_elems_host_id6.end();
        ++mit )
  {
    _marked_elems_host_id6[mit->first] = mit->second; //TODO do error checking for duplicates here too
  }
}

void
XFEMMeshCutUserObject::finalize()
{
  _communicator.set_union(_marked_elems_distance);
  _communicator.set_union(_marked_elems_host_id);
  _communicator.set_union(_marked_elems_distance2);
  _communicator.set_union(_marked_elems_host_id2);
  _communicator.set_union(_marked_elems_distance3);
  _communicator.set_union(_marked_elems_host_id3);
  _communicator.set_union(_marked_elems_distance4);
  _communicator.set_union(_marked_elems_host_id4);
  _communicator.set_union(_marked_elems_distance5);
  _communicator.set_union(_marked_elems_host_id5);
  _communicator.set_union(_marked_elems_distance6);
  _communicator.set_union(_marked_elems_host_id6);

  _xfem->clearUOMarkedElems();

  std::map<unsigned int, RealVectorValue >::iterator mit;
  //1
  for (mit = _marked_elems_host_id.begin(); mit != _marked_elems_host_id.end(); ++mit)
  {
   _xfem->addUOMarkedElemHostID(mit->first, mit->second, 0);
  }

  for (mit = _marked_elems_distance.begin(); mit != _marked_elems_distance.end(); ++mit)
  {
   _xfem->addUOMarkedElemDistance(mit->first, mit->second, 0);
  }

  //2
  for (mit = _marked_elems_host_id2.begin(); mit != _marked_elems_host_id2.end(); ++mit)
  {
    _xfem->addUOMarkedElemHostID(mit->first, mit->second, 1);
  }

  for (mit = _marked_elems_distance2.begin(); mit != _marked_elems_distance2.end(); ++mit)
  {
    _xfem->addUOMarkedElemDistance(mit->first, mit->second, 1);
  }

  //3
  for (mit = _marked_elems_host_id3.begin(); mit != _marked_elems_host_id3.end(); ++mit)
  {
    _xfem->addUOMarkedElemHostID(mit->first, mit->second, 2);
  }

  for (mit = _marked_elems_distance3.begin(); mit != _marked_elems_distance3.end(); ++mit)
  {
    _xfem->addUOMarkedElemDistance(mit->first, mit->second, 2);
  }

  //4
  for (mit = _marked_elems_host_id4.begin(); mit != _marked_elems_host_id4.end(); ++mit)
  {
    _xfem->addUOMarkedElemHostID(mit->first, mit->second, 3);
  }

  for (mit = _marked_elems_distance4.begin(); mit != _marked_elems_distance4.end(); ++mit)
  {
    _xfem->addUOMarkedElemDistance(mit->first, mit->second, 3);
  }

  //5
  for (mit = _marked_elems_host_id5.begin(); mit != _marked_elems_host_id5.end(); ++mit)
  {
    _xfem->addUOMarkedElemHostID(mit->first, mit->second, 4);
  }

  for (mit = _marked_elems_distance5.begin(); mit != _marked_elems_distance5.end(); ++mit)
  {
    _xfem->addUOMarkedElemDistance(mit->first, mit->second, 4);
  }

  //6
  for (mit = _marked_elems_host_id6.begin(); mit != _marked_elems_host_id6.end(); ++mit)
  {
    _xfem->addUOMarkedElemHostID(mit->first, mit->second, 5);
  }

  for (mit = _marked_elems_distance6.begin(); mit != _marked_elems_distance6.end(); ++mit)
  {
    _xfem->addUOMarkedElemDistance(mit->first, mit->second, 5);
  }

  _marked_elems_distance.clear();
  _marked_elems_host_id.clear();
  _marked_elems_distance2.clear();
  _marked_elems_host_id2.clear();
  _marked_elems_distance3.clear();
  _marked_elems_host_id3.clear();
  _marked_elems_distance4.clear();
  _marked_elems_host_id4.clear();
  _marked_elems_distance5.clear();
  _marked_elems_host_id5.clear();
  _marked_elems_distance6.clear();
  _marked_elems_host_id6.clear();
}
