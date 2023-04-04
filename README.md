# 3D-FDTD

### Current Features
- Uniform cubic mesh in cube simulation space
- Any number of current sources with arbitrary direction and position can be placed. 
- Mur absorbing boundary condition implemented
- A finite number of materials can be defined by using real parameters magnetic/electric conductivity, permittivity and permeability
- Model has same space resolution as the simulation grid
- Simulation timestep tied to spacial step to ensure stability and minimise dispersion

### Planned Features
- Support for superconducting materials
- Speed up simulation by use of parallelisation
- Rectangular simulation space with cube elements $N_x \neq N_y$ but $\Delta x = \Delta y$.
- uniform rectangular cuboid elements $N_x \neq N_y$ and $\Delta x \neq \Delta y$
- non-uniform rectangular cuboid elements
