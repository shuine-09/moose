# MovingLineSegmentCutSetUserObject

!syntax description /UserObjects/MovingLineSegmentCutSetUserObject

## Description

The `MovingLineSegmentCutSetUserObject` is used to cut the mesh with a set of line segments. The points on those line segments will be moving at a velocity that will be given by `XFEMPhaseTransitionMovingInterfaceVelocity`.

## Example Input File Syntax

!listing modules/xfem/test/tests/moving_interface/phase_transition.i block=UserObjects/moving_line_segments

!syntax parameters /UserObjects/MovingLineSegmentCutSetUserObject

!syntax inputs /UserObjects/MovingLineSegmentCutSetUserObject

!syntax children /UserObjects/MovingLineSegmentCutSetUserObject

!bibtex bibliography
