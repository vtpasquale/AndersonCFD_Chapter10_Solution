# Supersonic flow over a flat plate
Numerical solution by solving the complete Navier-Stokes equations. The problem is outlined in Chapter 10 of Anderson [1]. Spatial derivatives are calculated using finite difference and time marching is implemented using MacCormack method. 

The code is written for learning to solve this particular problem and nothing more. The code is structured differently than the flowcharts proposed by Anderson [1]. No FOR loops are used - this is intentional so the code works efficiently in MATLAB. Instead, array processing techniques are used for all computations. Solutions finish in 10-20 seconds for me, depending on boundary conditions and Courant number. 

The `Primitives.m` class is used to store velocity, pressure, and temperature, and to calculate dependent primitives (e.g., density, energy, viscosity). The independent primitives are selected so they are consistent with those specified at the boundaries. Constant gas properties are stored as class constants, so this is a gas and unit-system specific class implementation.

Anderson [1] Figures 10.10-10.15 can be replicated using `mainReproduceAndersonFigs.m`. Minor deviations from Anderson’s reference solutions are attributed to details that are not documented in the text. Details such as order of accuracy used for boundary derivative computations and Courant number. Additionally, there are typos in the equations for determining time step size, so there may be differences in step size compared to the reference solutions.

[1] Anderson, John David. Computational Fluid Dynamics. Colombia: McGraw-Hill Education, 1995.
