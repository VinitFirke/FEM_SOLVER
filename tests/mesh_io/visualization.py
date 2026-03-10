from supporting_files.vtu_output import VtuOutput
#from mesh_io import ReadMesh
from supporting_files.mesh_io import ReadMesh
import os

def main():
    os.makedirs("vtu_files", exist_ok=True)
    mesh = ReadMesh("mdpa_files/example_mesh_1.mdpa")
    VtuOutput("mdpa_files/example_mesh_1.mdpa", mesh, [], "", "vtu_files/output")


if __name__ == "__main__":
    main()