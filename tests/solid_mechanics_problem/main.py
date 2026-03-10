import numpy as np
from mesh_io import ReadMesh
from geometries.geometry import GeometryType
from elements.element import Element
from elements.solid_element_2d import SolidElement2D
from elements.force_nbc_2d import ForceNBC2D
from geometries.system_matrix_builder import SystemMatrixBuilder
from vtu_output import VtuOutput

def main():
    #mesh = ReadMesh("example_mesh_1.mdpa") #original
    mesh = ReadMesh("Task1_QuadrilateralMeshType_2D.mdpa")

    element_properties = {
        "poisson_ratio": 0.29,
        #"youngs_modulus": 2e9 #original
        "youngs_modulus": 30.5e9
    }

    nbc_properties = {
        #"traction": [1.0, 0.0]  #original
        "traction": [0.0, -50.0]
    }

    list_of_elements: list[Element] = []
    for geometry in mesh.geometries_list:
        if geometry.GetGeometryType() == GeometryType.Domain_integral:
            solid_element = SolidElement2D(geometry.geometry_id, geometry, element_properties)
            list_of_elements.append(solid_element)
        elif geometry.GetGeometryType() == GeometryType.Boundary_integral:
            force_nbc = ForceNBC2D(geometry.geometry_id, geometry, nbc_properties)
            list_of_elements.append(force_nbc)

    system_matrix_builder = SystemMatrixBuilder(list_of_elements)
    system_matrix_builder.InitializeDofs()

    #left_boundary_nodes = mesh.GetSubMesh("DISPLACEMENT_left").nodes_list   #original code, change related terms
    left_boundary_nodes_pwh = mesh.GetSubMesh("DISPLACEMENT_DBC").nodes_list
    for node in left_boundary_nodes_pwh:   #this one
        node.values["displacement_x"] = 0.0
        node.values["displacement_y"] = 0.0

    system_matrix_builder.BuildLHSandRHS()
    lhs = system_matrix_builder.global_lhs
    rhs = system_matrix_builder.global_rhs

    u = system_matrix_builder.GetSolution(mesh.nodes_list, ["displacement_x", "displacement_y"])

    system_matrix_builder.global_rhs = rhs - lhs@u
    system_matrix_builder.ApplyDirichletBoundaryConditions(left_boundary_nodes_pwh, ["displacement_x", "displacement_y"]) #this one
    du = np.linalg.solve(lhs, system_matrix_builder.global_rhs)

    system_matrix_builder.UpdateSolution(du, mesh.nodes_list, ["displacement_x", "displacement_y"])

    VtuOutput("Task1_QuadrilateralMeshType_2D.mdpa", mesh, ["displacement_x", "displacement_y"], "displacement", "output_pwh") #this one

if __name__ == "__main__":
    main()
