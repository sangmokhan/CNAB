The lid-driven cavity flow is one of the most frequently solved problems in the field of computational fluid dynamics, and it has been serving as a benchmark problem for testing in-house codes developed for incompressible viscous fluid flow.

In this repository, the lid-driven cavity problem is solved using the Crank-Nicolson/Adams-Bashforth scheme for viscous and convection terms, respectively. For spatial discretization, the second-order central difference scheme is used in a uniform staggered grid. The projection method is used to solve the incompressible Navier-Stokes equations, in which the computations of velocity and pressure are decoupled.

**Problem**: Consider two-dimensional unsteady flow in a square domain with Dirichlet boundary conditions on all sides, with three stationary sides and one moving side. Based on the lid velocity $U$ and the cavity height (or width) $H$, Reynolds number can be defined as Re = $\frac{UH}{\nu}$, where $\nu$ is the kinematic viscosity. The velocity of the impulsively started lid is given by a step function $u_{lid} = U$ for $t \geq 0$ and $u_{lid} = 0$ for $t < 0$. 

<p align="center">
<img width="400" alt="Screenshot 2024-07-18 at 2 03 34 PM" src="https://github.com/user-attachments/assets/fe545c64-e016-4406-a8c4-e34d76fc40ff">
</p>
