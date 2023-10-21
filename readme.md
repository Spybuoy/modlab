need different pvd files for various kinds of aspect ratios of length to breadth
(The code has 3:2 ratio with potentials on each pad being 4,-4,4,-4)

Steps - 
1. Make a folder with name ending with parameters. Ex - E_mag_300_200_10
2. Change values of V_pad in the code, Ex - (4,0,4,0)
3. Change the name of the PVD file to be saved describing the changed parameters
4. Upload the PVD and VTU files to the drive to the folder you've created before

 If possible, try to get slices from the visualization after checking if the solver's working using threshold filter

I've added the dirichlet conditions, having trouble with neumann conditionsj for now