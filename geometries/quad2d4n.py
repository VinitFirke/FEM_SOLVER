from geometries.geometry import Geometry, GeometryType
from geometries.node import Node
import numpy as np

class Quad2D4N(Geometry):
    def __init__(self,geometry_id: int, list_of_nodes: 'list[Node]'):
        super().__init__(geometry_id, list_of_nodes)

        if len(list_of_nodes) != 4:
            raise RuntimeError(f'Quad Geometry requires 4 nodes. Given: {len(list_of_nodes)}')
    def GetJacobianDeterminant(self, list_of_local_coordinates: 'list[float]') -> float:
        pass
    def GetGaussIntegrationPoints(self, integration_order: int) -> 'list[tuple[float, list[float]]]':


        if integration_order == 1:
            gauss_point_xi = 0.0
            gauss_point_eta  = 0.0
            gauss_point_weight  = 4.0

            return [(gauss_point_weight,  [gauss_point_xi, gauss_point_eta])]
        

        elif integration_order == 2:
            gauss_points: 'list[tuple[float, list[float]]]' = []

            #point 1
            gauss_point_weight = 1.0
            gauss_point_xi = -np.sqrt(3)/3
            gauss_point_eta = -np.sqrt(3)/3
            gauss_points.append((gauss_point_weight, [gauss_point_xi, gauss_point_eta]))

            #point 2
            gauss_point_weight = 1.0
            gauss_point_xi = -np.sqrt(3)/3
            gauss_point_eta = np.sqrt(3)/3
            gauss_points.append((gauss_point_weight, [gauss_point_xi, gauss_point_eta]))
        
            #point 3
            gauss_point_weight = 1.0
            gauss_point_xi = np.sqrt(3)/3
            gauss_point_eta =-np.sqrt(3)/3
            gauss_points.append((gauss_point_weight, [gauss_point_xi, gauss_point_eta]))

            #point 4
            gauss_point_weight = 1.0
            gauss_point_xi = np.sqrt(3)/3
            gauss_point_eta = np.sqrt(3)/3
            gauss_points.append((gauss_point_weight, [gauss_point_xi, gauss_point_eta]))
            return gauss_points
            
        else:
            raise RuntimeError(f"Unsupported integration order = {integration_order}")
        
    

    def GetShapeFunctionValue(self, list_of_local_coordinates: list[float]) -> list[float]:
        xi = list_of_local_coordinates[0]
        eta = list_of_local_coordinates[1]

        N_0 = (1 - xi)*(1 - eta)/4
        N_1 = (1 + xi)*(1 - eta)/4
        N_2 = (1 + xi)*(1 + eta)/4
        N_3 = (1 - xi)*(1 + eta)/4

        return [N_0, N_1, N_2, N_3]

    
    def GetJacobian(self, list_of_local_coordinates: 'list[float]') -> np.ndarray:
        xi = list_of_local_coordinates[0]
        eta = list_of_local_coordinates[1]
        
        n0 = self.list_of_nodes[0]
        n1 = self.list_of_nodes[1]
        n2 = self.list_of_nodes[2]
        n3 = self.list_of_nodes[3]

        dx_dxi  = -(1 - eta) * n0.x / 4 + (1 - eta) * n1.x / 4 + (1 + eta) * n2.x / 4  - (1 + eta) * n3.x / 4
        dx_deta = -(1 - xi)  * n0.x / 4 - (1 + xi)  * n1.x / 4 + (1 + xi)  * n2.x / 4  + (1 - xi)  * n3.x / 4
        dy_dxi  = -(1 - eta) * n0.y / 4 + (1 - eta) * n1.y / 4 + (1 + eta) * n2.y / 4  - (1 + eta) * n3.y / 4
        dy_deta = -(1 - xi)  * n0.y / 4 - (1 + xi)  * n1.y / 4 + (1 + xi)  * n2.y / 4  + (1 - xi)  * n3.y / 4 


        J = [
                 [dx_dxi, dx_deta],
                 [dy_dxi, dy_deta]
        ]
        numpy_jacobian = np.array(J)
        return numpy_jacobian
        

    def GetShapeFunctionDerivatives(self, list_of_local_coordinates: list[float]) -> np.ndarray:
        xi = list_of_local_coordinates[0]
        eta = list_of_local_coordinates[1]

        # dN0/dxi; N0 = (1 - xi) * (1 - eta) / 4
        dN0_dxi = -0.25 * (1 - eta)

        # dN1/dxi; N1 = (1 + xi) * (1 - eta) / 4
        dN1_dxi = 0.25 * (1 - eta)

        # dN2/dxi; N2 = (1 + xi) * (1 + eta) / 4
        dN2_dxi = 0.25 * (1 + eta)

        # dN2/dxi; N3 = (1 - xi) * (1 + eta) / 4
        dN3_dxi = -0.25 * (1 + eta)

        # dN0/deta; N0 = (1 - xi) * (1 - eta) / 4
        dN0_deta = -0.25 * (1 - xi)

        # dN1/deta; N1 = (1 + xi) * (1 - eta) / 4
        dN1_deta = -0.25 * (1 + xi)

        # dN2/deta; N2 = (1 + xi) * (1 + eta) / 4
        dN2_deta = 0.25 * (1 + xi)

        # dN3/deta; N3 = (1 - xi) * (1 + eta) / 4
        dN3_deta = 0.25 * (1 - xi)

        jacobian = self.GetJacobian(list_of_local_coordinates)

        
        inv_jacobian = np.linalg.inv(jacobian)

        local_shape_function_gradients = np.array([

            [dN0_dxi,   dN1_dxi,    dN2_dxi,  dN3_dxi],
            [dN0_deta,  dN1_deta,   dN2_deta, dN3_deta]
        ])

        global_shape_function_gradients = np.dot(inv_jacobian.T, local_shape_function_gradients)
        return global_shape_function_gradients
    
    def GetGeometryType(self) -> GeometryType:
        return GeometryType.Domain_integral




if __name__ == "__main__":
    node_1 = Node(1, 0, 0, 0)
    node_2 = Node(2, 2, 2, 0)
    node_3 = Node(3, 1, 4, 0)
    node_4 = Node(4,-1, 3, 0)

    ref_area = 6.5
    quad = Quad2D4N(1, [node_1, node_2, node_3, node_4])

    # checking for the integration order 1
    gauss_weight, (xi, eta) = quad.GetGaussIntegrationPoints(1)[0]
    print(quad.GetGaussIntegrationPoints(1)[0])
    jacobian = quad.GetJacobian([xi, eta])
    print(f"Jacobian: {jacobian}, Shape: {np.shape(jacobian)}")
    computed_area = gauss_weight * np.linalg.det(jacobian)

    
    
    computed_area = 0
    print(f"Order 2: {quad.GetGaussIntegrationPoints(2)}")
    for index, (gauss_weight, (xi, eta)) in enumerate(quad.GetGaussIntegrationPoints(2)):
        jacobian = quad.GetJacobian([xi, eta])
        print(f"jacobian {index}: {jacobian}")
        computed_area += gauss_weight * np.linalg.det(jacobian)


    ref_shape_function_derivatives = np.array(
        [
            [
                [-0.107069, 0.349897, 0.0286891, -0.271517],
                [-0.321208, 0.049691, 0.0860674, 0.185449]

            ],
            [
                [0.0539171, 0.124516, 0.286267, -0.4647],
                [-0.263479, -0.031129, 0.178433, 0.116175]
            ],
            [
                [-0.219925, 0.507895, -0.15188, -0.13609],
                [-0.195019, -0.126974, 0.28797, 0.0340225]
            ],
            [
                [-0.037509, 0.252513, 0.139986, -0.354989],
                [-0.112527, -0.242462, 0.419957, -0.0649675]
            ]
        ]
    )


    for index, (gauss_weight, (xi, eta)) in enumerate(quad.GetGaussIntegrationPoints(2)):
        computed_shape_function_derivatives = quad.GetShapeFunctionDerivatives([xi, eta])
        ref_gauss_point_shape_func_derivatives = ref_shape_function_derivatives[index]

        diff_matrix = computed_shape_function_derivatives - ref_gauss_point_shape_func_derivatives
        diff_matrix_norm = np.linalg.norm(diff_matrix)
        print(f"Diff matrix norm for gp {index} = {diff_matrix_norm}")

    print(f"Ref area = {ref_area}, computed area = {computed_area}, jacobian = {jacobian}")


    
