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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assemble material property vectors
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def k_vector_assemble():
  k_elem = mat_A_conductivity * np.ones(len(cut_mesh)-1)
  k_elem[np.searchsorted(base_mesh, mat_interf_loc):] = mat_B_conductivity

  return k_elem

def q_vector_assemble():
  q_elem = np.zeros(len(cut_mesh)-1)
  if mat_A_source == True:
    q_elem[:np.searchsorted(base_mesh, mat_interf_loc)] = mat_A_src_value
  if mat_B_source == True:
    q_elem[np.searchsorted(base_mesh, mat_interf_loc):] = mat_B_src_value

  return q_elem

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct System to Solve
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def construct_fem_system():
# Weak Form: 
# int_domain (k grad T grad b) = int_bd_domain (b k grad T n) + 
#     int_domain (b q)
# [K]{T} = {F} + {B}

  # Construct material vectors
  k_mat_elem = k_vector_assemble()
  q_mat_elem = q_vector_assemble()

  K = np.zeros((num_nodes, num_nodes))         # Stiffness Matrix
  F = np.zeros(num_nodes)                      # Force Vector

  # Assemble Stiffness Matrix and Force Vector
  for elem in range(num_elems):
    K_elem = np.zeros((porder+1,porder+1))
    F_elem = np.zeros(porder+1)
    J_elem = (cut_mesh[elem + 1] - cut_mesh[elem])/2
    for i in range(porder+1):
      for k in range(porder+1):
        F_elem[i] += b[0,i] * q_mat_elem[elem] * J_elem * wq[k]
      for j in range(porder+1):
        for k in range(porder+1):
          K_elem[i,j] += (k_mat_elem[elem] * dbdx[0,i]/J_elem * dbdx[1,j]/J_elem) *\
                         J_elem * wq[k]
    K[elem:elem+porder+1,elem:elem+porder+1] += K_elem
    F[elem:elem+porder+1] += F_elem

  # Apply Boundary Integral
  K, F = BC_apply(K, F)

  print(K)
  print('\n')
  print(F)
  print('\n')
  
  return (K, F)

def construct_xfem_system(xfem_method):
  K = np.zeros((num_nodes+1, num_nodes+1))
  F = np.zeros(num_nodes+1)

  # Determine number of split element
  cut_elnum = np.searchsorted(base_mesh, mat_interf_loc)

  # Forward Assembly of Stiffness Matrix and Force Vector to split element:
  if cut_elnum > 1:
    for felem in range(cut_elnum - 1):
      K_elem = np.zeros((porder+1,porder+1))
      F_elem = np.zeros(porder+1)
      J_elem = (base_mesh[felem + 1] - base_mesh[felem])/2
      for i in range(porder+1):
        if mat_A_source:
          for k in range(porder+1):
            F_elem[i] += b[0,i] * mat_A_src_value * J_elem * wq[k]
        for j in range(porder+1):
          for k in range(porder+1):
            K_elem[i,j] += (mat_A_conductivity * dbdx[0,i]/J_elem * dbdx[1,j]/J_elem) *\
                           J_elem * wq[k]
      K[felem:felem+porder+1,felem:felem+porder+1] += K_elem
      F[felem:felem+porder+1] += F_elem

  # Backward Assembly of Stiffness Matrix and Force Vector to split element:
  if cut_elnum < num_nodes:
    for belem in range(num_nodes-1, cut_elnum+1, -1):
      K_elem = np.zeros((porder+1,porder+1))
      F_elem = np.zeros(porder+1)
      J_elem = (base_mesh[(belem + 1)-2] - base_mesh[(belem)-2])/2
      for i in range(porder+1):
        if mat_B_source:
          for k in range(porder+1):
            F_elem[i] += b[0,i] * mat_B_src_value * J_elem * wq[k]
        for j in range(porder+1):
          for k in range(porder+1):
            K_elem[i,j] += (mat_B_conductivity * dbdx[0,i]/J_elem * dbdx[1,j]/J_elem) *\
                           J_elem * wq[k]
      K[belem:belem+porder+1,belem:belem+porder+1] += K_elem
      F[belem:belem+porder+1] += F_elem

  # Enriched Element Contribution
#  if xfem_method == 'XFEM-C':
#    K[cut_elnum-1:cut_elnum+3,cut_elnum-1:cut_elnum+3], F[cut_elnum-1:cut_elnum+3] = \
#        classic_xfem_cut_elem()
#  elif xfem_method == 'XFEM-PN':
#    K[cut_elnum-1:cut_elnum+3,cut_elnum-1:cut_elnum+3], F[cut_elnum-1:cut_elnum+3] = \
#        pn_xfem_cut_elem()
        
  K[cut_elnum-1:cut_elnum+3,cut_elnum-1:cut_elnum+3], F[cut_elnum-1:cut_elnum+3] = \
      pn_xfem_cut_elem(cut_elnum)

  # Apply Boundary Integral
  K, F = BC_apply(K, F)

  print(K)
  print('\n')
  print(F)
  print('\n')

  return (K, F)

#def classic_xfem_split_elem():
# Classical XFEM Literature Definition:
# u_h(x) = Sum_{j \in NEN} N_j(x) u_j + Sum_{j \in NEN^*} N_j^*(x) psi(x) a_j
#   NEN - Number elemental nodes
#   N_j - FE Basis Functions
#   psi - XFEM Enrichment Function
#   u_j - Storage vector for traditional FEM solution
#   a_j - Storage vector for enriched nodes
#   We use psi(x) = Heaviside(phi(x)) [step enrichments]
#  

def pn_xfem_cut_elem(cut_elnum):
# Phantom Node Literature Definition:
# u^h(x) = u_1^h(x) H(phi(x)) + u_2^h(x) H(-phi(x))
# Phantom Node and Classical XFEM Methods are functionally identical when
#   psi(x) = H(phi(x)), but the implementation is not the same.
  A_k = np.zeros((4,4))
  F_k = np.zeros(4)
  
  for pnode in range(0,len(F_k),2):
    J_elem = abs(mat_interf_loc - base_mesh[(cut_elnum-1)+pnode/2])/2
    for i in range(porder+1):
      for k in range(porder+1):
        if (pnode < 2) and mat_A_source:
          F_k[pnode+i] += b[0,i] * mat_A_src_value * J_elem * wq[k]
        elif (pnode >= 2) and mat_B_source:
          F_k[pnode+i] += b[0,i] * mat_B_src_value * J_elem * wq[k]
      for j in range(porder+1):
        for k in range(porder+1):
          if pnode < 2:
            A_k[pnode+i,pnode+j] += (mat_A_conductivity * dbdx[0,i]/J_elem * dbdx[1,j]/J_elem) *\
                                    J_elem * wq[k]
          elif pnode >= 2:
            A_k[pnode+i,pnode+j] += (mat_B_conductivity * dbdx[0,i]/J_elem * dbdx[1,j]/J_elem) *\
                                    J_elem * wq[k]

  return A_k, F_k

def BC_apply(A, b):
  if left_BC_type == 'Dirichlet':
    for i in range(len(b)):
      b[i] -= (A[i,0] * left_BC_value)
    A = np.delete(np.delete(A, 0 , 0), 0, 1)
    b = np.delete(b, 0)
  elif left_BC_type == 'Neumann':
    b[0] += left_BC_value
  else:
    print("Error!: Applicable Boundary Conditions are Dirichlet or Neumann!")

  if right_BC_type == 'Dirichlet':
    for i in range(len(b)):
      b[i] -= (A[i,-1] * right_BC_value)
    A = np.delete(np.delete(A, -1, 0), -1, 1)
    b = np.delete(b, -1)
  elif right_BC_type == 'Neumann':
    b[-1] += right_BC_value
  else:
    print("Error!: Applicable Boundary Conditions are Dirichlet on Neumann!")

  return (A, b)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear Solver
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Jacobi_l_solve(A, b):
  err = 1.0e12
  iter_count = 0
  
  x = np.zeros_like(b)

  while (iter_count < max_iterations) and (err > l_tol):
    iter_count += 1
    x_new = np.zeros_like(x)

    for i in range(len(x)):
      s1 = np.dot(A[i,:i], x[:i])
      s2 = np.dot(A[i,i + 1:], x[i + 1:])
      x_new[i] = (b[i] - s1 - s2) / A[i,i]
    err = np.linalg.norm(np.dot(A, x_new) - b)
    if l_output:
      print('Iteration Number: ' + str(iter_count) + ' L2 Error Norm: ' +\
            str(err))
    x = x_new

  return x

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
methods = ['FEM', 'XFEM-PN']                      # Options: FEM, XFEM-C, XFEM-PN
l_solver = 'Jacobi'                              # Options: Direct, Jacobi
l_tol = 1.0e-6                                   # Only used in iterative linear solve methods
max_iterations = 1.0e4                           # Maximum number of nonconverged linear iterations
l_output = False                                  # Toggle output of the linear solver
porder = 1                                       # Polynomial degree (must be 1)

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
num_nodes = len(cut_mesh)
num_elems = num_nodes - 1

print(base_mesh)
print(cut_mesh)
print(mat_interf_loc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Problem
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xq, wq = GLNodeWt(porder + 1)                    # Generate quadrature/weights
b, dbdx = feshpln(xq, porder)                    # Generate shape functs and their first derivatives
K = []
F = []
T = []

# Loop over solution methods
for i in range(len(methods)):
  K.append(np.zeros(1))
  F.append(np.zeros(1))
  T.append(np.zeros(1))
  if methods[i] == 'FEM':
    K[i], F[i] = construct_fem_system()
  elif methods[i] == 'XFEM-C' or methods[i] == 'XFEM-PN':
    K[i], F[i] = construct_xfem_system(methods[i])
  else:
    print('Invalid method input: options are FEM, XFEM-C, XFEM-PN')

  T_sub = Jacobi_l_solve(K[i], F[i])

  if left_BC_type == 'Dirichlet':
    T[i] = np.append(np.array(left_BC_value),T_sub)
  else:
    T[i] = T_sub
  if right_BC_type == 'Dirichlet':
    T[i] = np.append(T_sub, np.array(right_BC_value))

  np.savetxt('results_'+methods[i]+'.csv',T[i],fmt='%.8e')
