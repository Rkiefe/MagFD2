# Finite difference method for the demagnetizing field calculation
Solves the  2D demagnetizing field "Hd" for a given magnetization configuration using the finite difference method.

One script runs Hd considering only the magnetic domain. The other considers the vaccum exterior as well. There is an additional file ending in .c which contains the same algorithm written in C, with about 3x the performance. In the future, the C code will be implemented in Matlab using mex.

## Method:
Main equation:

$$\nabla \cdot \vec{H_d} = -\nabla \cdot \vec{M}$$

Considering $H_d = -\nabla \phi$, results in: $\nabla^2 \phi = \nabla \cdot \vec{M}$, which is a Poisson equation. Now you apply the centered differences to the second order derivatives and iterate over a grid.

### Potential of the magnetic field Hd, example:
![potential](https://github.com/user-attachments/assets/1b468598-1750-4ba6-9834-c1632d11e067)
