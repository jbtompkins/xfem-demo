# XFEM Demo
This is a demonstration of XFEM solution methods for 1D, linear, steady state
heat transfer boundary value problems with two discrete material regions
written in Python. 

## Features
The problem to be solved is a 1D, linear, steady state heat transfer BVP:
 -k (d^2 T)/(dx)^2 = q

Problem options for users to supply:
 - Starting and ending positions for the mesh
 - Number of initial intervals (elements) on the mesh
 - Location of the material interface
 - Left and right boundary conditions:
   - Dirichlet
   - Neumann
 - Thermal conductivity of material regions
 - Volumetric heat sources in each material region (if desired)

Users also have the choice between methods used to solve this problem:
 - FEM where the base mesh is cut and a standard FEM solve occurs
 - Classical XFEM solution
 - Phantom Node XFEM solution 

## Running the Code
Simply execute the main.py script using the python excutable after editing the
input variables to run a problem:

python main.py

## Input Variables
x_min:          [Real] The minimum location value included on the mesh.
x_max:          [Real] The maximum location value included on the mesh.
init_num_elems: [Int] The number of elements on the mesh before any cuts. 
mat_interf_loc: [Real] The location of the interface between the two materials.

left_BC_type:   [String] Type of boundary condition on the left end of the mesh.
                Options available are Neumann or Dirichlet.
left_BC_value:  [Real] Value of the left boundary condition.
right_BC_type:  [String] Type of boundary condition on the right end of the
                mesh. Options available are Neumann or Dirichlet.
right_BC_value: [Real] Value of the right boundary condition.

mat_A_conductivity: [Real] Thermal conductivity of the material left of the
                    interface.
mat_B_conductivity: [Real] Thermal conductivity of the material right of the
                    interface.
mat_A_source:       [bool] True if there is a source distributed through the
                    material left of the interface.
mat_A_src_value:    [Real](optional) Value of the left material source.
mat_B_source:       [bool] True if there is a source distributed through the
                    material right of the interface.
mat_B_src_value:    [Real](optional) Value of the right material source.

methods:        [[String] List] List of strings detailing the methods to run.
                Options available are FEM, XFEM-C, XFEM-PN.
l_solver:       [String] Type of linear solve to perform. Options available
                are Jacobi.
l_tol:          [Real] Convergence criterion for iterative linear solve methods
                (Jacobi).
max_iterations: [Int] Maximum number of nonconverged linear iterations.
l_output:       [bool] Toggle output of the linear solver.
