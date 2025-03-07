import meshio

mesh = meshio.read("mesh.msh")

print(mesh.cells)

# Cut out everything that isn't a tetrahedron
mesh.cells = [mesh.cells[-1]]

# Just to make sure
print(mesh.cells)

meshio.write("mesh.xdmf", mesh)
