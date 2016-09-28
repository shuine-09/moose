/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMMESHCUTBYFIELD_H
#define XFEMMESHCUTBYFIELD_H

#include "XFEMMeshCutUserObject.h"

class XFEMMeshCutByField;

template<>
InputParameters validParams<XFEMMeshCutByField>();

class XFEMMeshCutByField : public XFEMMeshCutUserObject
{
 public:
  XFEMMeshCutByField(const InputParameters & parameters);

  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutEdge> & cut_edges);
  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutFace> & cut_faces);

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_edges,
                                     std::vector<CutEdge> & cut_edges);
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_faces,
                                     std::vector<CutFace> & cut_faces);
  
  virtual void execute(); 

 protected:
  
  Real _contour;

  /// name of field variable
  NonlinearVariableName _var_name;

  /// variable
  MooseVariable & _var;

  NonlinearSystem & _system;
  const NumericVector<Number> * _solution;

};

#endif // XFEMMESHCUTBYFIELD_H
