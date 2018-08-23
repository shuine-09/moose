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

#ifndef COHESIVEZONEMESH_H
#define COHESIVEZONEMESH_H

#include "MooseMesh.h"

// forward declaration
class CohesiveZoneMesh;

template <>
InputParameters validParams<CohesiveZoneMesh>();

class CohesiveZoneMesh : public MooseMesh
{
public:
  CohesiveZoneMesh(const InputParameters & parameters);
  CohesiveZoneMesh(const CohesiveZoneMesh & other_mesh);
  virtual ~CohesiveZoneMesh(); // empty dtor required for unique_ptr with forward declarations

  virtual void init() override;

  virtual std::unique_ptr<MooseMesh> safeClone() const override;

  virtual void buildMesh() override;

  void read(const std::string & file_name);
  virtual ExodusII_IO * exReader() const override { return _exreader.get(); }

  // Get/Set Filename (for meshes read from a file)
  void setFileName(const std::string & file_name) { _file_name = file_name; }

  virtual std::string getFileName() const override { return _file_name; }

  void buildNodeSupport();
  void buildInterfacialNodes();
  void duplicateNodes();
  void tearElements();
  void addInterfaceBoundary();

protected:
  /// the file_name from whence this mesh came
  std::string _file_name;
  /// Auxiliary object for restart
  std::unique_ptr<ExodusII_IO> _exreader;

private:
  MeshBase * _file_mesh;
  std::map<dof_id_type, std::vector<dof_id_type>> _node_support;
  std::map<dof_id_type, unsigned int> _node_duplicity;
  std::map<dof_id_type, std::vector<dof_id_type>> _duplicated_node;
};

#endif // COHESIVEZONEMESH_H
