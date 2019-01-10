//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POLYCRYSTALVORONOIBUBBLEIC_H
#define POLYCRYSTALVORONOIBUBBLEIC_H

#include "MultiSmoothCircleIC.h"
#include "MooseRandom.h"
#include "PolycrystalICTools.h"

// Forward Declarationsc
class PolycrystalVoronoiBubbleIC;

template <>
InputParameters validParams<PolycrystalVoronoiBubbleIC>();

/**
 * PolycrystalVoronoiBubbleIC initializes either grain or void values for a
 * voronoi tesselation with voids distributed along the grain boundaries.
 */
class PolycrystalVoronoiBubbleIC : public MultiSmoothCircleIC
{
public:
  PolycrystalVoronoiBubbleIC(const InputParameters & parameters);

  virtual void initialSetup() override;

  static InputParameters actionParameters();

protected:
  const MooseEnum _structure_type;

  const unsigned int _op_num;
  unsigned int _grain_num;
  const unsigned int _op_index;

  const unsigned int _rand_seed;

  const bool _columnar_3D;

  virtual void computeCircleCenters() override;

  virtual Real value(const Point & p) override;
  virtual RealGradient gradient(const Point & p) override;

  virtual Real grainValueCalc(const Point & p);
  virtual void computeGrainCenters();

  std::vector<Point> _centerpoints;
  std::vector<unsigned int> _assigned_op;

  std::vector<Point> _gb_normal;

  Real _R0;
  Real _r0;

  FileName _file_name;

  /// Type for distance and point
  struct DistancePoint
  {
    Real d;
    unsigned int gr;
  };

  /// Sorts the temp_centerpoints into order of magnitude
  struct DistancePointComparator
  {
    bool operator()(const DistancePoint & a, const DistancePoint & b) { return a.d < b.d; }
  } _customLess;
};

#endif // POLYCRYSTALVORONOIBUBBLEIC_H
