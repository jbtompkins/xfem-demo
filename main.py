import numpy as np

# ******************************************************************************
# XFEM Demo Code
# The point of this exercise is to demonstrate various solution methods
# utilizing FEM or XFEM to solve a 1D, linear, steady state heat conduction
# boundary value problem with two discrete material regions. The following
# options are left to the user to supply:
#  - Starting and ending positions for the mesh
#  - Number of initial intervals (elements) on the mesh
#  - Location of the material interface
#  - Left and right boundary condtions:
#    - Dirichlet
#    - Neumann
#  - Thermal conductivity of material regions
#  - Volumetric heat sources in each material region (if desired)
#  - Solution method options:
#    - FEM where the base mesh is cut and std FEM solve occurs
#    - Mesh is not cut; Classical XFEM solution
#    - Mesh is not cut; Phantom Node XFEM solution
# ******************************************************************************

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input/Problem Specification
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mesh Properties
x_min = 0.
x_max = 1.
delta_x = 0.2
init_num_elems = 5
mat_interf_loc = 0.5

# Boundary Condtitions (Options: Neumann, Dirichlet)
left_BC_type = 'Neumann'
left_BC_value = 0.0
right_BC_type = 'Dirichlet'
right_BC_value = 400.0

# Material Properties
mat_A_conductivity = 0.3
mat_B_conductivity = 0.02
mat_A_source = True
mat_A_src_value = 2.0                            # Optional Variable
mat_B_source = True
mat_B_src_value = 3.5                            # Optional Variable

# Solve Type (Options: FEM, XFEM-C, XFEM-PN)
methods = ['FEM', 'XFEM-C']


