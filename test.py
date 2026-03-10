def ReadMesh(file_name: str):
    with open(file_name,"r") as file_input:
        lines = file_input.readlines()

    for element_line in lines[26 + 1: 36]:
        print(element_line)


if __name__ == "__main__":
    ReadMesh("example_mesh_1.mdpa")