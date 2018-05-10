# LevelSetBiMaterialProperty

!syntax description /Materials/LevelSetBiMaterialProperty

## Description

This material, `LevelSetBiMaterialProperty` is intended only for use with XFEM. It determines the global material property by switching the two material properties with different base_name based on the level set values.

## Example Input File Syntax

!listing modules/xfem/test/tests/moving_interface/moving_diffusion.i block=Materials/diff_combined

!syntax parameters /Materials/LevelSetBiMaterialProperty

!syntax inputs /Materials/LevelSetBiMaterialProperty

!syntax children /Materials/LevelSetBiMaterialProperty

!bibtex bibliography
