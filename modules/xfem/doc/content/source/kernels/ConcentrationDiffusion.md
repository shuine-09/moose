# ConcentrationDiffusion

!syntax description /Kernels/ConcentrationDiffusion

## Description

The `ConcentrationDiffusion` computes residual/Jacobian contribution for diffusion equation with diffusivity. It has an optional parameter `base_name` that allows the users to define multiple materials systems on the smae block.

## Example Input File Syntax

!listing modules/xfem/test/tests/moving_interface/moving_diffusion.i block=Kernels/diff

!syntax parameters /Kernels/ConcentrationDiffusion

!syntax inputs /Kernels/ConcentrationDiffusion

!syntax children /Kernels/ConcentrationDiffusion

!bibtex bibliography
