import numpy as np
from supporting_files.mesh_io import ReadMesh
from elements.heat_conduction_element import HeatConductionElement
from geometries.system_matrix_builder import SystemMatrixBuilder
from supporting_files.vtu_output import VtuOutput
import os
from geometries.geometry import GeometryType

def main():
    os.makedirs("vtu_files", exist_ok=True)
    mesh = ReadMesh("mdpa_files/example_mesh_1.mdpa")

    list_of_elements: 'list[HeatConductionElement]' = []

    #for geometry in mesh.geometries_list:
    #    heat_conduction_element = HeatConductionElement(geometry.geometry_id, geometry, {})
    #    list_of_elements.append(heat_conduction_element)
    for geometry in mesh.geometries_list:
        if geometry.GetGeometryType() != GeometryType.Domain_integral:
            continue
        heat_conduction_element = HeatConductionElement(geometry.geometry_id, geometry, {})
        list_of_elements.append(heat_conduction_element)
    
    system_matrix_builder = SystemMatrixBuilder(list_of_elements)

    system_matrix_builder.InitializeDofs()

    left_boundary_nodes = mesh.GetSubMesh("DISPLACEMENT_left").nodes_list
    for node in left_boundary_nodes:
        node.values["temperature"] = 10.0
    
    right_boundary_nodes = mesh.GetSubMesh("LineLoad2D_right").nodes_list
    for node in right_boundary_nodes:
        node.values["temperature"] = 100.0


    system_matrix_builder.BuildLHSandRHS()
    lhs = system_matrix_builder.global_lhs
    rhs = system_matrix_builder.global_rhs

    # Ku = F
    u = system_matrix_builder.GetSolution(mesh.nodes_list, ["temperature"])

    #R = F - Ku
    # R_i = F_i - K_ij x u_j

    system_matrix_builder.global_rhs = rhs - lhs@u
      
    system_matrix_builder.ApplyDirichletBoundaryConditions(left_boundary_nodes, ['temperature'])
    system_matrix_builder.ApplyDirichletBoundaryConditions(right_boundary_nodes, ["temperature"])
    
    du = np.linalg.solve(lhs, system_matrix_builder.global_rhs)

    system_matrix_builder.UpdateSolution(du, mesh.nodes_list, ["temperature"])

    VtuOutput("mdpa_files/example_mesh_1.mdpa", mesh, ["temperature"], "temperature", "vtu_files/heat_conduction_output")


if __name__ == "__main__":
    main()