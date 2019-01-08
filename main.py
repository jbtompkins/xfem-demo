import numpy as np
from math import sqrt

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

def input_params():
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Input/Problem Specification
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
  
  # Methods and Solvers Options
  methods = ['FEM', 'XFEM-C']                      # Options: FEM, XFEM-C, XFEM-PN
  l_solver = 'Jacobi'                              # Options: Direct, Jacobi
  l_tol = 1.0e-6                                   # Only used in iterative linear solve methods
  porder = 1                                       # Polynomial degree (must be 1 for now)
  
  # Additional Variables
  eps = 1.0e-8
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Mesh Construction
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Construct base, equally spaced mesh (used in XFEM calculations)
  delta_x = (x_max - x_min)/init_num_elems
  base_mesh = np.arange(x_min, x_max + eps, delta_x)
  
  # Generate cut mesh (used in FEM calculations)
  cut_mesh = np.insert(base_mesh, np.searchsorted(base_mesh, mat_interf_loc),
                       mat_interf_loc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct System to Solve
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Weak Form: 
# int_domain (grad T k grad b) - int_bd_domain (b k grad t n) = 
#     int_domain (b q)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Quadrature Point and Weight Generation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def GLNodeWt(n):
# Generates Gauss-Legendre Quadrature points and weights
  beta = []
  for i in range(1,n):
    beta.append(i/sqrt(4*pow(i,2) - 1))
  J = np.diag(beta,-1) + np.diag(beta,1)
  D, V = np.linalg.eig(J)
  x = np.sort(D)
  w = 2*V[0,:]**2

  return (x,w)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FE Basis Functions and Derivatives
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def feshpln(xv, p):
# Computes the finite element basis functions and their first derivatives for
# any order 'p' using Lagrangian FE functions
  xd = np.linspace(-1, 1, num = p+1)

  shapefunc = np.zeros((len(xv), p+1))
  dhdx = np.zeros((len(xv), p+1))

  # Shape function
  for i in range(0, p+1):
    num = 1.
    den = 1.
    for j in range(0, p+1):
      if j != i:
        num *= (xv - xd[j])
        den *= (xd[i] - xd[j])
    shapefunc[:,i] = num/den
  
  # Derivative of the shape function
  for i in range(0, p+1):
    sum_tot = 0.
    den = 1.
    for j in range(0, p+1):
      if j != i:
        num = 1.
        for k in range(0, p+1):
          if (k == i and k == j):
            num *= (xv - xd[k])
        sum_tot += num
        den *= (xd[i] - xd[j])
    dhdx[:,i] = sum_tot/den

  return (shapefunc, dhdx)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Problem
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input_params()
xq, wq = GLNodeWt(porder + 1)                    # Generate quadrature/weights
b, dbdx = feshpln(xq, porder)                    # Generate shape functs and their first derivatives

