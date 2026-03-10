from geometries.geometry import Geometry, GeometryType
from geometries.node import Node
import numpy as np

class Triangle2D3N(Geometry):
    def __init__(self, geometry_id: int, list_of_nodes: list[Node]):
        super().__init__(geometry_id, list_of_nodes)

        if len(list_of_nodes) != 3:
            raise RuntimeError(f"Triangle geometry requires 3 nodes. Given: {len(list_of_nodes)}")

    def GetGaussIntegrationPoints(self, integration_order: int) -> 'list[tuple[float, list[float]]]':
        
       if integration_order == 1:
           gauss_point_weight = 0.5
           gauss_point_xi = 1/3
           gauss_point_eta = 1/3

           return [(gauss_point_weight,  [gauss_point_xi, gauss_point_eta])]

       elif integration_order ==2:
           gauss_points: 'list[tuple[float, list[float]]]' = []

           #first gauss point
           gauss_point_weight = 1/6
           gauss_point_xi = 1/6
           gauss_point_eta = 1/6
           gauss_points.append((gauss_point_weight, [gauss_point_xi, gauss_point_eta]))

           #first gauss point
           gauss_point_weight = 1/6
           gauss_point_xi = 2/3
           gauss_point_eta = 1/6
           gauss_points.append((gauss_point_weight, [gauss_point_xi, gauss_point_eta]))

           #first gauss point
           gauss_point_weight = 1/6
           gauss_point_xi = 1/6
           gauss_point_eta = 2/3
           gauss_points.append((gauss_point_weight, [gauss_point_xi, gauss_point_eta]))
           return gauss_points
       else:
           raise RuntimeError(f"Unsupported integration order = {integration_order}")
            
        
         
        #match integration_order:
        #    case 0: print("H")


    def GetShapeFunctionValue(self, list_of_local_coordinates: list[float]) -> list[float]:
        xi = list_of_local_coordinates[0]
        eta = list_of_local_coordinates[1]

        N_0 = 1 - xi - eta
        N_1 = xi
        N_2 = eta

        return [N_0, N_1, N_2]

    
    def GetJacobian(self, list_of_local_coordinates: 'list[float]') -> np.ndarray:
        
        dx_dxi  = self.list_of_nodes[1].x - self.list_of_nodes[0].x
        dx_deta = self.list_of_nodes[2].x - self.list_of_nodes[0].x
        dy_dxi  = self.list_of_nodes[1].y - self.list_of_nodes[0].y
        dy_deta = self.list_of_nodes[2].y - self.list_of_nodes[0].y                    

        J = [
                 [dx_dxi, dx_deta],
                 [dy_dxi, dy_deta]
        ]

        numpy_jacobian = np.array(J)

        return numpy_jacobian
    def GetJacobianDeterminant(self, list_of_local_coordinates: 'list[float]') -> float:
        pass


    def GetShapeFunctionDerivatives(self, list_of_local_coordinates: list[float]) -> np.ndarray:
        dN0_dxi = -1
        dN1_dxi =  1
        dN2_dxi =  0

        dN0_deta = -1
        dN1_deta =  0
        dN2_deta =  1
        
        jacobian = self.GetJacobian(list_of_local_coordinates)
        inv_jacobian = np.linalg.inv(jacobian)

        local_shape_function_gradients = np.array([

            [dN0_dxi,   dN1_dxi,    dN2_dxi],
            [dN0_deta,  dN1_deta,   dN2_deta]
        ])

        global_shape_function_gradients = np.dot(inv_jacobian.T, local_shape_function_gradients)
        return global_shape_function_gradients
    
    def GetGeometryType(self) -> GeometryType:
        return GeometryType.Domain_integral

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

    triangle = Triangle2D3N([node_1, node_2, node_3])

    gauss_weight, [xi, eta] = triangle.GetGaussIntegrationPoints(1)[0]

    jacobian = triangle.GetJacobian([xi, eta])

    if np.linalg.norm(ref_jacobian - jacobian, 'fro') == 0.0:
        
        print('The jacobian is correct')
    else:
        print('The jacobian is incorrect')