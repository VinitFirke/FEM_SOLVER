import KratosMultiphysics as Kratos
import numpy as np
from mesh_io import Mesh

def VtuOutput(file_name: str, mesh: Mesh, dofs_list: 'list[str]', output_dof_name: str, output_file_name: str):
    model = Kratos.Model()
    model_part = model.CreateModelPart("vtu_output")
    Kratos.ModelPartIO(file_name[:-5], Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)

    vtu_output = Kratos.VtuOutput(model_part)
    cexp = Kratos.Expression.NodalExpression(model_part)

    if len(dofs_list) != 0:
        values: 'list[list[float]]' = []
        for node in mesh.nodes_list:
            nodal_values: 'list[float]' = []
            for dof in dofs_list:
                nodal_values.append(node.values[dof])
            values.append(nodal_values)

        Kratos.Expression.CArrayExpressionIO.Read(cexp, np.array(values, dtype=np.float64))
        vtu_output.AddContainerExpression(output_dof_name, cexp)
    vtu_output.PrintOutput(output_file_name)

if __name__ == "__main__":
    VtuOutput(1, 1)