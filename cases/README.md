# Creating Case Files
A *case* consists of a simulation input file (`.json`) and an associated Gmsh mesh file (`.msh`). The JSON file defines the physical model, numerical methods, output options, and the assignment of materials, sources, and boundary conditions to regions in the mesh. The mesh file defines the computational geometry and finite element discretization. The JSON input file specifies:
1. Material properties
2. Boundary conditions
3. Sources
4. Region assignments
5. Angular treatment
6. Spatial discretization
7. Numerical solver parameters
8. Output options

The Gmsh mesh file specifies:
1. Node coordinates
2. Element connectivity
3. Named "Physical Groups" identifying material and source regions and boundaries.

The mesh file, alternatively, contains information about the simulation geometry, including the node placements and element connectivity. Requirements on the mesh include:
1. `gmsh` *must* be used, and the output should be saved in ASCII (`.msh`) format.
2. The mesh must consist *entirely* of quadrilateral (2D) or hexahedral (3D) elements.
3. Boundary regions must be defined using Gmsh Physical Groups:
   - Physical Point for 1D problems,
   - Physical Curve for 2D problems,
   - Physical Surface for 3D problems.
4. Material regions must be defined using Gmsh Physical Groups:
   - Physical Curve for 1D problems,
   - Physical Surface for 2D problems,
   - Physical Volume for 3D problems.
5. In the current implementation, the names of material-source region Physical Groups must exactly match the corresponding region names specified in the JSON file. These names are used to associate mesh regions with material and source definitions.

