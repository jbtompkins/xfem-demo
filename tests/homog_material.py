# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# XFEM Demo Test
# Name: homog_material
# Purpose: Test that has the same material properties on both sides of the
#   material interface. Can be used to verify code using Method of Manufactured
#   Solutions. For now, used to ensure the system is constructed correctly in
#   various methods.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mesh Properties
x_min = 0.
x_max = 1.
init_num_elems = 5
mat_interf_loc = 0.5

# Boundary Conditions (Options: Neumann, Dirichlet)
left_BC_type = 'Neumann'
left_BC_value = 0.0
right_BC_type = 'Dirichlet'
right_BC_value = 400.0

# Material Properties
mat_A_conductivity = 0.3
mat_B_conductivity = 0.3
mat_A_source = True
mat_A_src_value = 2.0                            # Optional Variable
mat_B_source = True
mat_B_src_value = 2.0                            # Optional Variable

# Methods and Solvers Options
methods = ['FEM', 'XFEM-PN']                     # Options: FEM, XFEM-C, XFEM-PN
l_solver = 'Jacobi'                              # Options: Direct, Jacobi
l_tol = 1.0e-6                                   # Only used in iterative linear solve methods
max_iterations = 1.0e4                           # Maximum number of nonconverged linear iterations
l_output = False                                 # Toggle output of the linear solver
