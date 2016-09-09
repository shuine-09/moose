/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMMESHCUTBYLEVELSET_H
#define XFEMMESHCUTBYLEVELSET_H

#include "XFEMMeshCutUserObject.h"

class XFEMMeshCutByLevelSet;

template<>
InputParameters validParams<XFEMMeshCutByLevelSet>();

class XFEMMeshCutByLevelSet : public XFEMMeshCutUserObject
{
 public:
  XFEMMeshCutByLevelSet(const InputParameters & parameters);

  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutEdge> & cut_edges);
  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutFace> & cut_faces);

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_edges,
                                     std::vector<CutEdge> & cut_edges);
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_faces,
                                     std::vector<CutFace> & cut_faces);
  
  virtual void execute(); 

 protected:
  
  /// name of level set variable
  NonlinearVariableName _ls_var_name;

  /// level set variable
  MooseVariable & _ls_var;

  SystemBase & _aux_system;
  const NumericVector<Number> * _aux_solution;

};

#endif // XFEMMESHCUTBYLEVELSET_H
