/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                        Grizzly                               */
/*                                                              */
/*           (c) 2015 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ReadGrainCenterFromFile.h"

template<>
InputParameters validParams<ReadGrainCenterFromFile>()
{
  InputParameters params = validParams<ElementPropertyReadFile>();
  params.addParam<std::string>("grain_center_file_name", "Name of the file with grain center");
  return params;
}

ReadGrainCenterFromFile::ReadGrainCenterFromFile(const InputParameters & parameters) :
    ElementPropertyReadFile(parameters),
    _grain_center_file_name(getParam<std::string>("grain_center_file_name"))
{
  initGrainCenterPoints();
}

void
ReadGrainCenterFromFile::initGrainCenterPoints()
{
  readFromFile();
}

void
ReadGrainCenterFromFile::readFromFile()
{
  _center.resize(_ngrain);
  MooseUtils::checkFileReadable(_grain_center_file_name);

  std::ifstream file_prop;
  file_prop.open(_grain_center_file_name.c_str());

  for ( unsigned int i = 0; i < _ngrain; ++i)
    for ( unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if (!(file_prop >> _center[i](j)))
        mooseError("Error ReadGrainCenterFromFile: Premature end of file");

  file_prop.close();
}
