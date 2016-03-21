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

#ifndef READGRAINCENTERFROMFILE
#define READGRAINCENTERFROMFILE

#include "ElementPropertyReadFile.h"

class ReadGrainCenterFromFile;

template<>
InputParameters validParams<ReadGrainCenterFromFile>();

class ReadGrainCenterFromFile : public ElementPropertyReadFile
{
 public:
  ReadGrainCenterFromFile(const InputParameters & parameters);
  virtual ~ReadGrainCenterFromFile() {}

  virtual void initGrainCenterPoints();
  void readFromFile();

 protected:
  const std::string _grain_center_file_name;
};

#endif
