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

#ifndef XFEMCRACKTIPMEANSTRESS_H
#define XFEMCRACKTIPMEANSTRESS_H

#include "ElementUserObject.h"
#include "libmesh/string_to_enum.h"
#include "MaterialTensorCalculator.h"

// libMesh
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/elem.h"

class XFEM;

class XFEMCrackTipMeanStress : public ElementUserObject
{
  public:

    XFEMCrackTipMeanStress(const InputParameters & parameters);

    virtual ~XFEMCrackTipMeanStress() {}

    virtual void initialize();
    virtual void execute();
    virtual void threadJoin(const UserObject &y);
    virtual void finalize();

  protected:
    virtual std::vector<Real> getStressTensor();
    unsigned int _crack_front_point_index;

  private:
    MaterialTensorCalculator _material_tensor_calculator;
    const MaterialProperty<SymmTensor> & _tensor;
    std::vector<Real> _stress_tensor;
    Real _radius;
    MooseSharedPointer<XFEM> _xfem;
    std::map<unsigned int, const Elem* > _elem_id_crack_tip;
    std::vector<Point> _crack_front_points;
    std::vector<Real> _weights;
    unsigned int _num_crack_front_points;
    Real _critical_stress;
};

template<>
InputParameters validParams<XFEMCrackTipMeanStress>();

#endif //XFEMCRACKTIPMEANSTRESS_H
