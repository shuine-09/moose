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

#ifndef RESTARTABLEDATAIO_H
#define RESTARTABLEDATAIO_H

// MOOSE includes
#include "DataIO.h"

// C++ includes
#include <sstream>
#include <string>
#include <list>

// Forward declarations
class RestartableDatas;
class RestartableDataValue;
class FEProblem;

/**
 * Helper class to hold streams for Backup and Restore operations.
 */
class Backup
{
public:
  Backup()
  {
    unsigned int n_threads = libMesh::n_threads();

    _restartable_data.resize(n_threads);

    for (unsigned int i=0; i < n_threads; i++)
      _restartable_data[i] = new std::stringstream;
  }

  ~Backup()
  {
    unsigned int n_threads = libMesh::n_threads();

    for (unsigned int i=0; i < n_threads; i++)
      delete _restartable_data[i];
  }

  std::stringstream _system_data;

  std::vector<std::stringstream*> _restartable_data;
};

template<>
inline void
dataStore(std::ostream & stream, Backup * & backup, void * context)
{
  dataStore(stream, backup->_system_data, context);

  for (unsigned int i=0; i<backup->_restartable_data.size(); i++)
    dataStore(stream, backup->_restartable_data[i], context);
}

template<>
inline void
dataLoad(std::istream & stream, Backup * & backup, void * context)
{
  dataLoad(stream, backup->_system_data, context);

  for (unsigned int i=0; i<backup->_restartable_data.size(); i++)
    dataLoad(stream, backup->_restartable_data[i], context);
}


/**
 * Class for doing restart.
 *
 * It takes care of writing and reading the restart files.
 */
class RestartableDataIO
{
public:
  RestartableDataIO(FEProblem & fe_problem);

  virtual ~RestartableDataIO();

  /**
   * Write out the restartable data.
   */
  void writeRestartableData(std::string base_file_name, const RestartableDatas & restartable_datas, std::set<std::string> & _recoverable_data);

  /**
   * Read restartable data header to verify that we are restarting on the correct number of processors and threads.
   */
  void readRestartableDataHeader(std::string base_file_name);

  /**
   * Read the restartable data.
   */
  void readRestartableData(const RestartableDatas & restartable_datas, const std::set<std::string> & _recoverable_data);

  /**
   * Create a Backup for the current system.
   */
  MooseSharedPointer<Backup> createBackup();

  /**
   * Restore a Backup for the current system.
   */
  void restoreBackup(MooseSharedPointer<Backup> backup, bool for_restart = false);

private:
  /**
   * Serializes the data into the stream object.
   */
  void serializeRestartableData(const std::map<std::string, RestartableDataValue *> & restartable_data, std::ostream & stream);

  /**
   * Deserializes the data from the stream object.
   */
  void deserializeRestartableData(const std::map<std::string, RestartableDataValue *> & restartable_data, std::istream & stream, const std::set<std::string> & recoverable_data);

  /**
   * Serializes the data for the Systems in FEProblem
   */
  void serializeSystems(std::ostream & stream);

  /**
   * Deserializes the data for the Systems in FEProblem
   */
  void deserializeSystems(std::istream & stream);

  /// Reference to a FEProblem being restarted
  FEProblem & _fe_problem;

  /// A vector of file handles, one per thread
  std::vector<std::ifstream *> _in_file_handles;
};

#endif /* RESTARTABLEDATAIO_H */
