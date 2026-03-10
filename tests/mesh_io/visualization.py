from vtu_output import VtuOutput
#from mesh_io import ReadMesh
from mesh_io import ReadMesh

def main():
    mesh = ReadMesh("example_mesh_1.mdpa")
    VtuOutput("example_mesh_1.mdpa", mesh, [], "", "output")


if __name__ == "__main__":
    main()