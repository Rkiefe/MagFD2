# Finite difference method for the demagnetizing field calculation
Solves the  2D demagnetizing field "Hd" for a given magnetization configuration using the finite difference method.

It considers the vaccum exterior to set the magnetic potential to 0 in "infinity". There is an additional file ending in .c which contains the same algorithm written in C, with about 3x the performance, still in production. In the future, the C code will be implemented in Matlab using mex.

## Method:
Main equation:

$$\nabla \cdot \vec{H_d} = -\nabla \cdot \vec{M}$$

Considering $H_d = -\nabla \phi$, results in: $\nabla^2 \phi = \nabla \cdot \vec{M}$, which is a Poisson equation. Now you apply the centered differences to the second order derivatives and iterate over a grid.
The border condition for the magnetic region is,
$$(\nabla \phi_{int} - \nabla \phi_{ext} )\cdot n = M \cdot n$$
and for the outside: $\phi = 0$

### Potential of the magnetic field Hd, example:
![potential](https://github.com/user-attachments/assets/1b468598-1750-4ba6-9834-c1632d11e067)
