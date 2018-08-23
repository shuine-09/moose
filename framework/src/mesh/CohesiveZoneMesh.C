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

#include "CohesiveZoneMesh.h"
#include "Parser.h"
#include "MooseUtils.h"
#include "Moose.h"
#include "MooseApp.h"

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/serial_mesh.h"

registerMooseObject("MooseApp", CohesiveZoneMesh);

template <>
InputParameters
validParams<CohesiveZoneMesh>()
{
  InputParameters params = validParams<MooseMesh>();
  params.addRequiredParam<MeshFileName>("file", "The name of the mesh file to read");
  params.addClassDescription("Read a mesh from a file.");
  return params;
}

CohesiveZoneMesh::CohesiveZoneMesh(const InputParameters & parameters)
  : MooseMesh(parameters), _file_name(getParam<MeshFileName>("file"))
{
  getMesh().set_mesh_dimension(getParam<MooseEnum>("dim"));
}

CohesiveZoneMesh::CohesiveZoneMesh(const CohesiveZoneMesh & other_mesh)
  : MooseMesh(other_mesh), _file_name(other_mesh._file_name)
{
}

CohesiveZoneMesh::~CohesiveZoneMesh() {}

std::unique_ptr<MooseMesh>
CohesiveZoneMesh::safeClone() const
{
  return libmesh_make_unique<CohesiveZoneMesh>(*this);
}

void
CohesiveZoneMesh::init()
{
  // read the temporary mesh from Exodus file
  // std::string file_name = getParam<MeshFileName>("file");
  // MooseUtils::checkFileReadable(file_name);
  // _file_mesh = new SerialMesh(_communicator);
  // ExodusII_IO * exodusII_io = new ExodusII_IO(*_file_mesh);
  // exodusII_io->read(file_name);
  // getMesh().allow_renumbering(true);
  // getMesh().prepare_for_use(/*true*/);
  //
  // // build neighborlist for _fe_mesh elements
  // // getMesh().find_neighbors();
  //
  // // initialize unconventional mesh data for PD mesh
  // // std::cout << "file name = " << file_name << std::endl;
  // std::cout << "WJ : dimension = " << getMesh().mesh_dimension() << std::endl;
  // std::cout << "WJ : number of elements = " << getMesh().n_elem() << std::endl;
  // std::cout << "WJ : number of nodes = " << getMesh().n_nodes() << std::endl;

  // for (MeshBase::node_iterator it = getMesh().nodes_begin(); it != getMesh().nodes_end();
  // ++it)
  // {
  //   Node * node = *it;
  //   std::cout << "node = " << *node << std::endl;
  // }

  // for (MeshBase::element_iterator it = getMesh().elements_begin();
  //      it != getMesh().elements_end();
  //      ++it)
  // {
  //   Elem * elem = *it;
  //   std::cout << "elem subdomain id = " << elem->subdomain_id() << std::endl;
  // }

  // buildNodeSupport();
  // buildInterfacialNodes();
  // duplicateNodes();
  // tearElements();

  // for (MeshBase::element_iterator it = getMesh().elements_begin();
  //      it != getMesh().elements_end();
  //      ++it)
  // {
  //   Elem * elem = *it;
  //   std::cout << "elem = " << *elem << std::endl;
  // }

  MooseMesh::init();

  buildNodeSupport();
  buildInterfacialNodes();
  duplicateNodes();
  tearElements();
  addInterfaceBoundary();
}

void
CohesiveZoneMesh::buildNodeSupport()
{
  for (MeshBase::element_iterator it = getMesh().elements_begin(); it != getMesh().elements_end();
       ++it)
  {
    Elem * elem = *it;
    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
      _node_support[elem->node_id(i)].push_back(elem->id());
  }

  // for (auto & ns : _node_support)
  // {
  //   std::cout << "node id = " << ns.first << std::endl;
  //   for (unsigned int i = 0; i < (ns.second).size(); ++i)
  //     std::cout << "===> elem id = " << (ns.second)[i] << std::endl;
  // }
}

void
CohesiveZoneMesh::buildInterfacialNodes()
{
  const unsigned int n_nodes = getMesh().n_nodes();
  const unsigned int n_elems = getMesh().n_elem();

  std::set<subdomain_id_type> matSet;

  for (MeshBase::node_iterator it = getMesh().nodes_begin(); it != getMesh().nodes_end(); ++it)
  {
    Node * node = *it;
    std::vector<dof_id_type> support = _node_support[node->id()];
    unsigned int suppCount = support.size();

    if (suppCount == 1)
      continue;

    for (unsigned int ie = 0; ie < suppCount; ie++)
    {
      subdomain_id_type imat = getMesh().elem(support[ie])->subdomain_id();
      matSet.insert(imat);
    }

    unsigned int matCount = matSet.size();

    // if interface elements are created everywhere then
    // duplicity = nodal support. Else duplicity = materials.

    unsigned int duplicity = matCount;

    if (matCount != 1) // interfacial node
    {
      _node_duplicity[node->id()] = duplicity;
    }

    matSet.clear(); // clear for the next node

    // check if this node is on an external boundary or not

    // if (find(globdat.boundaryNodes.begin(), globdat.boundaryNodes.end(), index) !=
    //     globdat.boundaryNodes.end())
    // {
    //   np->setIsOnBoundary(true);
    // }
  }
}

void
CohesiveZoneMesh::duplicateNodes()
{
  int idd(0);

  for (auto & dn : _node_duplicity)
  {
    _duplicated_node[dn.first].push_back(dn.first);

    for (int id = 0; id < dn.second - 1; id++)
    {
      idd++;

      Node * new_node = Node::build(getMesh().node(dn.first), getMesh().n_nodes()).release();
      new_node->processor_id() = getMesh().node(dn.first).processor_id();
      getMesh().add_node(new_node);

      std::cout << "added new node " << new_node->id() << std::endl;
      _duplicated_node[dn.first].push_back(new_node->id());
    }
  }

  std::cout << "number of nodes added: " << idd << std::endl;
}

void
CohesiveZoneMesh::tearElements()
{
  for (auto & dn : _duplicated_node)
  {
    std::vector<dof_id_type> support = _node_support[dn.first];
    unsigned int suppCount = support.size();

    subdomain_id_type mat = getMesh().elem(support[0])->subdomain_id();

    for (unsigned int ie = 0; ie < suppCount; ie++)
    {
      subdomain_id_type imat = getMesh().elem(support[ie])->subdomain_id();

      if (imat == mat)
        continue;

      Elem * elem = getMesh().elem(support[ie]);

      for (unsigned int i = 0; i < elem->n_nodes(); ++i)
      {
        if (elem->node_id(i) == dn.first)
        {
          elem->set_node(i) = getMesh().node_ptr((dn.second)[1]);
        }
      }
    }
  }
}

void
CohesiveZoneMesh::addInterfaceBoundary()
{
  BoundaryInfo & boundary_info = getMesh().get_boundary_info();

  for (MeshBase::element_iterator it = getMesh().elements_begin(); it != getMesh().elements_end();
       ++it)
  {
    Elem * elem = *it;
    for (unsigned int i = 0; i < elem->n_sides(); ++i)
    {
      unsigned int j = i + 1;
      if (j == elem->n_nodes())
        j = 0;
      Node * node1 = elem->get_node(i);
      Node * node2 = elem->get_node(j);

      if (_duplicated_node.find(node1->id()) != _duplicated_node.end() &&
          _duplicated_node.find(node2->id()) != _duplicated_node.end())
      {
        std::cout << "Added elem id = " << elem->id() << ", edge = " << i << " to the sideset 100"
                  << std::endl;
        boundary_info.add_side(elem->id(), i, 100);
      }
    }
  }
  boundary_info.sideset_name(100) = "interface";

  // getMesh().find_neighbors();

  for (MeshBase::element_iterator it = getMesh().elements_begin(); it != getMesh().elements_end();
       ++it)
  {
    const Elem * elem = *it;
    for (unsigned int i = 0; i < elem->n_sides(); ++i)
    {
      const Elem * neighbor = elem->neighbor_ptr(i);
      if (neighbor != nullptr)
      {
        // std::cout << "elem id = " << elem->id() << ", neighbor id " << neighbor->id() <<
        // std::endl; std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>ELEM" << std::endl; std::cout <<
        // *elem << std::endl; std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>NEIGHBOR_ELEM" <<
        // std::endl; std::cout << *neighbor << std::endl;
      }
    }
  }

  // getMesh().prepare_for_use(/*skip_renumber =*/true);
}

void
CohesiveZoneMesh::buildMesh()
{
  std::string _file_name = getParam<MeshFileName>("file");

  Moose::perf_log.push("Read Mesh", "Setup");
  if (_is_nemesis)
  {
    // Nemesis_IO only takes a reference to DistributedMesh, so we can't be quite so short here.
    DistributedMesh & pmesh = cast_ref<DistributedMesh &>(getMesh());
    Nemesis_IO(pmesh).read(_file_name);

    getMesh().allow_renumbering(false);

    // Even if we want repartitioning when load balancing later, we'll
    // begin with the default partitioning defined by the Nemesis
    // file.
    bool skip_partitioning_later = getMesh().skip_partitioning();
    getMesh().skip_partitioning(true);
    getMesh().prepare_for_use();
    getMesh().skip_partitioning(skip_partitioning_later);
  }
  else // not reading Nemesis files
  {
    // See if the user has requested reading a solution from the file.  If so, we'll need to read
    // the mesh with the exodus reader instead of using mesh.read().  This will read the mesh on
    // every processor

    if (_app.setFileRestart() && (_file_name.rfind(".exd") < _file_name.size() ||
                                  _file_name.rfind(".e") < _file_name.size()))
    {
      MooseUtils::checkFileReadable(_file_name);

      _exreader = libmesh_make_unique<ExodusII_IO>(getMesh());
      _exreader->read(_file_name);

      getMesh().allow_renumbering(false);
      getMesh().prepare_for_use();
    }
    else
    {
      auto slash_pos = _file_name.find_last_of("/");
      auto path = _file_name.substr(0, slash_pos);
      auto file = _file_name.substr(slash_pos + 1);

      bool restarting;
      // If we are reading a mesh while restarting, then we might have
      // a solution file that relies on that mesh partitioning and/or
      // numbering.  In that case, we need to turn off repartitioning
      // and renumbering, at least at first.
      if (file == "LATEST")
      {
        std::list<std::string> dir_list(1, path);
        std::list<std::string> files = MooseUtils::getFilesInDirs(dir_list);

        // Fill in the name of the LATEST file so we can open it and read it.
        _file_name = MooseUtils::getLatestMeshCheckpointFile(files);
        restarting = true;
      }
      else
        restarting = _file_name.rfind(".cpa") < _file_name.size() ||
                     _file_name.rfind(".cpr") < _file_name.size();

      const bool skip_partitioning_later = restarting && getMesh().skip_partitioning();
      const bool allow_renumbering_later = restarting && getMesh().allow_renumbering();

      if (restarting)
      {
        getMesh().skip_partitioning(true);
        getMesh().allow_renumbering(false);
      }

      MooseUtils::checkFileReadable(_file_name);
      getMesh().read(_file_name);

      if (restarting)
      {
        getMesh().allow_renumbering(allow_renumbering_later);
        getMesh().skip_partitioning(skip_partitioning_later);
      }
    }
  }

  std::cout << "WJ2 : dimension = " << getMesh().mesh_dimension() << std::endl;
  std::cout << "WJ2 : number of elements = " << getMesh().n_elem() << std::endl;
  std::cout << "WJ2 : number of nodes = " << getMesh().n_nodes() << std::endl;

  Moose::perf_log.pop("Read Mesh", "Setup");
}

void
CohesiveZoneMesh::read(const std::string & file_name)
{
  if (dynamic_cast<DistributedMesh *>(&getMesh()) && !_is_nemesis)
    getMesh().read(file_name, /*mesh_data=*/NULL, /*skip_renumber=*/false);
  else
    getMesh().read(file_name, /*mesh_data=*/NULL, /*skip_renumber=*/true);
}
