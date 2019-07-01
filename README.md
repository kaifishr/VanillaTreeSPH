# VanillaTreeSPH
A simple smoothed-particle hydrodynamics code that uses a tree algorithm to find neighbour particles. Using a tree algorithm is not always necessary to compute the direct interactions between particles. For the example below, an approach with linked lists would be sufficient and also speed up the computation. However, a tree structure makes sense if long-range potentials are used. This would be, for example, the gravitational interaction between particles.

Each particle in a SPH simulation does not only have a location, a velocity and a mass like in a normal N-body simulation, but can also have pressure and density. 

Here is a simple example of a collision of two bodies consisting of particles of the same mass. The small body and the large body consist of 100 and 10000 particles, respectively. The initial configuration is shown below:

<p align="center">
    <img src="https://github.com/KaiFabi/VanillaTreeSPH/blob/master/results/res_id_0.png" height="500">
</p>

Here the color of the particles corresponds to the object they belong to. In the following results the particles were colored according to their speed, pressure, and density:

**Velocity**
<p align="center">
    <img src="https://github.com/KaiFabi/VanillaTreeSPH/blob/master/results/sph_vel.gif" height="500">
</p>

**Pressure**
<p align="center">
    <img src="https://github.com/KaiFabi/VanillaTreeSPH/blob/master/results/sph_pressure.gif" height="500">
</p>

**Density**
<p align="center">
    <img src="https://github.com/KaiFabi/VanillaTreeSPH/blob/master/results/sph_rho.gif" height="500">
</p>


Compile and run the program using

`gcc -O -Wall sph.c -o sph -lm && ./sph`
