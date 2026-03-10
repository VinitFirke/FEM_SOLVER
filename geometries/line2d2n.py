from geometries.geometry import Geometry, GeometryType
from geometries.node import Node
import numpy as np
from geometries.triangle2d3n import Triangle2D3N


class Line2D2N(Geometry):
    def __init__(self, geometry_id: int, list_of_nodes: list[Node]):
        super().__init__(geometry_id, list_of_nodes)

        if len(list_of_nodes) != 2:
            raise RuntimeError(f"Line geometry requires 2 nodes. Given: {len(list_of_nodes)}")

    def GetGaussIntegrationPoints(self, integration_order: int) -> 'list[tuple[float, list[float]]]':
        gauss_points: 'list[tuple[float, list[float]]]' = []
        if integration_order == 1:
            gauss_point_weight = 2.0
            gauss_point_xi = 0.0
            gauss_points.append((gauss_point_weight, [gauss_point_xi]))
        elif integration_order == 2:
            # gauss point 1
            gauss_point_weight = 1.0
            gauss_point_xi = 1 / (3 ** 0.5)
            gauss_points.append((gauss_point_weight, [gauss_point_xi]))

            # gauss point 2
            gauss_point_weight = 1.0
            gauss_point_xi = -1 / (3 ** 0.5)
            gauss_points.append((gauss_point_weight, [gauss_point_xi]))
        else:
            raise RuntimeError(f"Unsupported integration order = {integration_order}.")

        return gauss_points
        
    
    def GetShapeFunctionValue(self, list_of_local_coordinates: 'list[float]') -> 'list[float]':
        xi = list_of_local_coordinates[0]

        n_0 = (1 + xi) / 2
        n_1 = (1 - xi) / 2

        return [n_0, n_1]

    
    
    def GetJacobian(self, list_of_local_coordinates: 'list[float]') -> np.ndarray:
        raise NotImplementedError("This is not possible yet with the Line2D2N")
    
    def GetJacobianDeterminant(self, list_of_local_coordinates: 'list[float]') -> float:
        node_1 = self.list_of_nodes[0]
        node_2 = self.list_of_nodes[1]

        line_length = ((node_2.x - node_1.x)** 2 + (node_2.y - node_1.y)** 2) ** (0.5)
        return 0.5 * line_length

    
    def GetShapeFunctionDerivatives(self, list_of_local_coordinates: list[float]) -> np.ndarray:
       raise NotImplementedError("This is not possible yet with the Line2D2N")
    
    def GetGeometryType(self) -> GeometryType:
        return GeometryType.Boundary_integral

    def __str__(self):
        msg = "Line2D2N with Nodes: ["
        for node in self.list_of_nodes:
            msg += " " + str(node.node_id)
        msg += "]"
        return msg
    

if __name__ == "__main__":
    node_1 = Node(1,0,0,0)
    node_2 = Node(2,2,2,0)
    node_3 = Node(3,-1,3,0)
    

    ref_jacobian = np.array(
        [
            [2, -1],
            [2,  3]
        ]
    )

    triangle = Triangle2D3N(1,[node_1, node_2, node_3])

    gauss_weight, [xi, eta] = triangle.GetGaussIntegrationPoints(1)[0]

    jacobian = triangle.GetJacobian([xi, eta])

    if np.linalg.norm(ref_jacobian - jacobian, 'fro') == 0.0:
        
        print('The jacobian is correct')
    else:
        print('The jacobian is incorrect')