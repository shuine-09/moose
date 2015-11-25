/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYUOBASE_H
#define CRYSTALPLASTICITYUOBASE_H

#include "ElementUserObject.h"
#include "RankTwoTensor.h"

class CrystalPlasticityUOBase;

template<>
InputParameters validParams<CrystalPlasticityUOBase>();

/**
 * Crystal plasticity system userobject base class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticityUOBase : public ElementUserObject
{
 public:
  CrystalPlasticityUOBase(const InputParameters & parameters);

  virtual ~CrystalPlasticityUOBase() {}

  void initialize() {}
  void execute() {}
  void finalize() {} 
  void threadJoin(const UserObject &) {}

  /// Returns the slip rate  name  (eg "DislocationGlide")
  virtual std::string crystalPlasticityUOName() const;

  /// Returns the size of data
  virtual unsigned int dataSize() const;

 protected:

  std::string _uo_name;

  unsigned int _uo_data_size;
};

#endif // CRYSTALPLASTICITYUOBASE_H
