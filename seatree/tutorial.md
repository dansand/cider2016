## HC



1. Start up the hc module. Use the default settings, but load viscosity structure visc.D from the filedialog. Compute mantle flow at several depths, and plot the geoid.
   
   ``` 
   #Dimesnionless depth, viscosity
   0.546225 5e22
   0.8964 1e21
   ...
   ```
   
   ​
   
2. Change the asthenospheric and lower mantle viscosity. Compute the geoid and comment on the changes, for example as expressed in the correlation with the actual geoid.
   
3. By adding another layer, shift the lower mantle viscosity increase from being at 660 km to
   
   ~1000 km. Recompute the geoid.
   
   * significantly reducing the depth of this viscosity jump eventually shifted the geoid around downwelling to negative
   * Also, using a uniform viscosity for the lower and upper mantle, around 10**21, means the geoid goes negative at significant downwellings.
   
4. Change the viscosity jump into a viscosity ramp, from 660 to 1000 km depth, for example.
   
5. Reload the original viscosity structure, for example by restarting SEATREE or by navigating to
   
   the default input data directory $SEATREE/seatree/python/data/hc/viscosity, where $SEATREE
   
   is the installation directory. Change to depth-dependent scaling of velocity to density and
   
   compare the geoid (and the radial tractions at ~100 km) with constant scaling with that where
   
   the upper or the lower 300 km of the mantle are set to zero. How about the lower 500 km of the
   
   mantle zeroed out?
   
6. Using your favorite flow model with free slip surface boundary conditions, plot the velocities at
   
   the surface, at 600 km (layer 26 for the default model), at 750 km (in the lower mantle, layer
   
   23), and just above the CMB (layer 2).
   
7. Change to prescribed plate velocities and replot the velocities, bottom up.
   
8. Compare the geoid for free-slip and prescribed plate computations.
   
9. Save the velocity output for your favorite model in VTK, and use Paraview 3-D visualization to
   
   explore the velocity and density structure.
   
   This covers some of the global dynamics questions discussed in class. However, you might also
   
   consider the next exercises.