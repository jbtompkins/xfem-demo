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
Create a python file with the input variables to describe the desired problem
options, then edit line 6 and line 8 in main.py to point to the directory of
the input file (comment out the line if input is in the same directory as
main.py) and the name of the input file (without .py) respectively.

Then simply execute the main.py script using the python excutable:

python main.py

## Input Variables
x_min:          [Real] The minimum location value included on the mesh.<br/>
x_max:          [Real] The maximum location value included on the mesh.<br/>
init_num_elems: [Int] The number of elements on the mesh before any cuts.<br/>
mat_interf_loc: [Real] The location of the interface between the two materials.

left_BC_type:   [String] Type of boundary condition on the left end of the mesh. Options available are Neumann or Dirichlet.<br/>
left_BC_value:  [Real] Value of the left boundary condition.<br/>
right_BC_type:  [String] Type of boundary condition on the right end of the mesh. Options available are Neumann or Dirichlet.<br/>
right_BC_value: [Real] Value of the right boundary condition.

mat_A_conductivity: [Real] Thermal conductivity of the material left of the interface.<br/>
mat_B_conductivity: [Real] Thermal conductivity of the material right of the interface.<br/>
mat_A_source:       [bool] True if there is a source distributed through the material left of the interface.<br/>
mat_A_src_value:    [Real](optional) Value of the left material source.<br/>
mat_B_source:       [bool] True if there is a source distributed through the material right of the interface.<br/>
mat_B_src_value:    [Real](optional) Value of the right material source.

methods:        [[String] List] List of strings detailing the methods to run. Options available are FEM, XFEM-C, XFEM-PN.<br/>
l_solver:       [String] Type of linear solve to perform. Options available are Jacobi.<br/>
l_tol:          [Real] Convergence criterion for iterative linear solve methods (Jacobi).<br/>
max_iterations: [Int] Maximum number of nonconverged linear iterations.<br/>
l_output:       [bool] Toggle output of the linear solver.<br/>

## Tests
Several tests are found in the 'tests' folder that evaluate or verify the code's
performance in various situations:
 - homog_material
 - middle_interface
 - left_interface
 - right_interface
 - node_interface
