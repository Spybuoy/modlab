from fenics import *

path = "/mnt/c/Files/modlab/code/"

"""
** Make changes to a copy, not the org**
** Only change values in BLOCK_1 **
** May change equations in BLOCK_2 **
"""

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
---------------------BLOCK_1--------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
# mesh dimensions
mesh_side = 600e-6
cells_x, cells_y, cells_z = 20, 20, 20

# Define the dimensions of the bot and pads
bot_x, bot_y, bot_z = 300e-6, 200e-6, 10e-6
pad_x, pad_y, pad_z = 50e-6, 50e-6, 10e-6

# Define boundary conditions for the cuboid
V_bn_value = 4
V_pad1, V_pad2, V_pad3, V_pad4 = 4, -4, 4, -4

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
---------------------BLOCK_1_END----------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
---------------------BLOCK_2--------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# Define the center of the bot
center = Point(0, 0, 0)

# Define the center of pads
p1_center = Point(bot_x / 4.0, bot_y / 4.0, 0)
p2_center = Point(-bot_x / 4.0, bot_y / 4.0, 0)
p3_center = Point(-bot_x / 4.0, -bot_y / 4.0, 0)
p4_center = Point(bot_x / 4.0, -bot_y / 4.0, 0)
"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
---------------------BLOCK_2_END----------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
---------------------BLOCK_3--------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

outer_mesh = BoxMesh(
    Point(-mesh_side/2.0, -mesh_side/2.0, -mesh_side/2.0),
    Point(mesh_side/2.0, mesh_side/2.0, mesh_side/2.0),
    cells_x, cells_y, cells_z,
)

# Define separate subdomains for pads and cuboid
class Pad1(SubDomain):
    def inside(self, x, on_boundary):
        return (
            (abs(x[0] - p1_center[0]) <= pad_x / 2)
            and (abs(x[1] - p1_center[1]) <= pad_y / 2)
            and (abs(x[2] - p1_center[2]) <= pad_z / 2)
            and (x[0] > 0)
            and (x[1] > 0)
        )


class Pad2(SubDomain):
    def inside(self, x, on_boundary):
        return (
            (abs(x[0] - p2_center[0]) <= pad_x / 2)
            and (abs(x[1] - p2_center[1]) <= pad_y / 2)
            and (abs(x[2] - p2_center[2]) <= pad_z / 2)
            and (x[0] < 0)
            and (x[1] > 0)
        )


class Pad3(SubDomain):
    def inside(self, x, on_boundary):
        return (
            (abs(x[0] - p3_center[0]) <= pad_x / 2)
            and (abs(x[1] - p3_center[1]) <= pad_y / 2)
            and (abs(x[2] - p3_center[2]) <= pad_z / 2)
            and (x[0] < 0)
            and (x[1] < 0)
        )


class Pad4(SubDomain):
    def inside(self, x, on_boundary):
        return (
            (abs(x[0] - p4_center[0]) <= pad_x / 2)
            and (abs(x[1] - p4_center[1]) <= pad_y / 2)
            and (abs(x[2] - p4_center[2]) <= pad_z / 2)
            and (x[0] > 0)
            and (x[1] < 0)
        )


# # Define the subdomain for the cuboid
# class Cuboid(SubDomain):
#     def inside(self, x, on_boundary):
#         return (
#             not any(
#                 [
#                     Pad1().inside(x, on_boundary),
#                     Pad2().inside(x, on_boundary),
#                     Pad3().inside(x, on_boundary),
#                     Pad4().inside(x, on_boundary),
#                 ]
#             )
#         )



# Define function space
V = FunctionSpace(outer_mesh, "P", 1)  

# Define trial and test functions
phi = TrialFunction(V)
q = TestFunction(V)

# Initialize mesh function for the subdomains
markers = MeshFunction("size_t", outer_mesh, 
outer_mesh.topology().dim())
markers.set_all(0)

# Mark the subdomain for the cuboid 
# cuboid = Cuboid()
# cuboid.mark(markers, 1)

pad1 = Pad1()
pad1.mark(markers, 1)
pad2 = Pad2()
pad2.mark(markers, 2)
pad3 = Pad3()
pad3.mark(markers, 3)
pad4 = Pad4()
pad4.mark(markers, 4)


# Define a new measure for integration over the bot
# dx_bot = Measure("dx", domain=outer_mesh, subdomain_data=markers)
dx_pad1 = Measure("dx", domain=outer_mesh, subdomain_data=markers)
dx_pad2 = Measure("dx", domain=outer_mesh, subdomain_data=markers)
dx_pad3 = Measure("dx", domain=outer_mesh, subdomain_data=markers)
dx_pad4 = Measure("dx", domain=outer_mesh, subdomain_data=markers)
# DX = (dx_bot + dx_pad1 + dx_pad2 + dx_pad3 + dx_pad4)
DX = (dx_pad1 + dx_pad2 + dx_pad3 + dx_pad4)

# Defining boundary conditions for the pads
#bc_bot = [DirichletBC(V, Constant(V_bn_value), cuboid)]
bc_p1 = [DirichletBC(V, Constant(V_pad1), pad1)]
bc_p2 = [DirichletBC(V, Constant(V_pad2), pad2)]
bc_p3 = [DirichletBC(V, Constant(V_pad3), pad3)]
bc_p4 = [DirichletBC(V, Constant(V_pad4), pad4)]

# Define weak forms
a_phi = inner(grad(phi), grad(q)) * DX
L_phi = (
    Constant(0) * q * DX
    
)  # Assuming no source term in this case
# Equivalent to ∇Φ^2=0

# Solve for Phi
Phi = Function(V)
BCS = bc_p1 + bc_p2 + bc_p3 + bc_p4
solve(a_phi == L_phi, Phi, BCS)

# Define the electric field function space
V_vec = VectorFunctionSpace(outer_mesh, "P", 1)

# Define the electric field
E = Function(V_vec)

# Define trial and test functions for the EF
E_trial = TrialFunction(V_vec)
v = TestFunction(V_vec)

# Define the weak form for the electric field
a_E = inner(E_trial, v) * DX
L_E = -inner(grad(Phi), v) * DX

# Solve for the electric field
solve(a_E == L_E, E)

# Calculate the magnitude of the electric field
E_magnitude = sqrt(dot(E, E))

# Project the magnitude onto the function space
E_magnitude_func = project(E_magnitude, V)

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
---------------------BLOCK_3_END----------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# Save the magnitude of the electric field
output_file = File(path 
+ "oct19/sort_of_final_draft/E_mag_calc2.pvd")
output_file << E_magnitude_func

print("File Executed")
