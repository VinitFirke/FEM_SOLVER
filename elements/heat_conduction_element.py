import numpy as np
from elements.element import Element
from geometries.node import Node
from geometries.geometry import Geometry, GeometryType
from typing import Any
from geometries.quad2d4n import Quad2D4N

class HeatConductionElement(Element):
    def __init__(self, element_id: int, geometry: Geometry, parameters: 'dict[str, Any]'):
        super().__init__(element_id, geometry, parameters)
        if geometry.GetGeometryType() != GeometryType.Domain_integral:
            raise RuntimeError(f"{self.__class__.__name__} only supports domain integrals")

    def CalculateLeftHandSideMatrix(self) -> np.ndarray:
        
        number_of_nodes = len(self.geometry.list_of_nodes)
        lhs = np.zeros((number_of_nodes, number_of_nodes))

        #construct the lhs
        #   1. Get the gauss point data
        #   2. Iterate through each gauss point
        #   3. Construct the weak formulation for each gauss point
        #   4. Build the lhs
       
        
        # 2. Iterate through each gauss point
        for gauss_weight, (xi, eta) in self.geometry.GetGaussIntegrationPoints(2):
            jacobian = self.geometry.GetJacobian([xi, eta])
            det_jacobian = np.linalg.det(jacobian)
            dN_dx = self.geometry.GetShapeFunctionDerivatives([xi, eta])

            # loop to compute dNa_dNb
            dNa_dNb = dN_dx.T @ dN_dx

            # loop to compute lhs
            for a in range(number_of_nodes):
                for b in range(a + 1, number_of_nodes):
                    value = gauss_weight * dNa_dNb[a, b] * det_jacobian
                    lhs[a, b] += value #upper half of 4x4 matrix
                    lhs[b, a] += value #lower half of 4x4 matrix
                lhs[a, a] += gauss_weight * det_jacobian * dNa_dNb[a, a] #diagonal of 4x4 matrix 

        return lhs 
                   
    

    def CalculateRightHandSideVector(self) -> np.ndarray:
        rhs = np.zeros(len(self.geometry.list_of_nodes))
        return rhs
    
    def GetDofIdsVector(self) -> list[int]:
        dof_ids: 'list[int]' = []
        for node in self.geometry.list_of_nodes:
            dof_ids.append(node.dof_ids["temperature"])
        return dof_ids
    
    def GetDofNames(self) -> 'list[str]':
        return ["temperature"]
        

"""
import numpy as np
from unittest.mock import MagicMock

if __name__ == "__main__":
    # Mock geometry
    class MockGeometry:
        def __init__(self, nodes):
            self.list_of_nodes = nodes
        
        def GetGeometryType(self):
            return GeometryType.Domain_integral
        
        def GetGaussIntegrationPoints(self, order: int):
            # Return dummy Gauss points and weights
            return [(1.0, (0.0, 0.0)), (1.0, (1.0, 1.0))]  # weight, (xi, eta)
        
        def GetJacobian(self, xi_eta):
            # Return a simple identity matrix for the Jacobian
            return np.eye(2)
        
        def GetShapeFunctionDerivatives(self, xi_eta):
            # Return dummy derivatives of shape functions
            return np.array([[1.0, 0.0], [0.0, 1.0]])
        
    
    # Mock node
    class MockNode:
        def __init__(self, dof_id):
            self.dof_ids = {"temperature": dof_id}
    
    # Create mock nodes
    mock_nodes = [MockNode(dof_id=1), MockNode(dof_id=2)]
    
    # Create mock geometry
    mock_geometry = MockGeometry(mock_nodes)
    
    # Create parameters
    parameters = {"thermal_conductivity": 1.0}
    
    # Instantiate HeatConductionElement
    element = HeatConductionElement(1, mock_geometry, parameters)
    
    # Test CalculateLeftHandSideVector
    lhs = element.CalculateLeftHandSideMatrix()
    print("LHS Matrix:")
    print(lhs)
    
    # Test CalculateRightHandSideVector
    rhs = element.CalculateRightHandSideVector()
    print("RHS Vector:")
    print(rhs)
    
    # Test GetDofIdsVector
    dof_ids = element.GetDofIdsVector()
    print("DOF IDs:")
    print(dof_ids)
    
    # Simple checks
    if lhs.shape == (2, 2):
        print("LHS shape test passed.")
    else:
        print("LHS shape test failed.")
    
    if np.all(rhs == 0):
        print("RHS test passed.")
    else:
        print("RHS test failed.")
    
    if dof_ids == [1, 2]:
        print("DOF IDs test passed.")
    else:
        print("DOF IDs test failed.")
    

"""


if __name__ == "__main__":
    node_1 = Node(1,  0, 0, 0)
    node_2 = Node(2,  2, 2, 0)
    node_3 = Node(3,  1, 4, 0)
    node_4 = Node(4, -1, 3, 0)

    quad = Quad2D4N(1,[node_1, node_2, node_3, node_4])
    element = HeatConductionElement(1, quad, {})

    ref_lhs = np.array([0.0] * 16).reshape((4,4))
    print(f"ref_lhs before:\n {ref_lhs}")
    for gauss_weight, (xi, eta) in quad.GetGaussIntegrationPoints(2):
        jacobian = quad.GetJacobian([xi,eta])
        det_jacobian = np.linalg.det(jacobian)
        dN_dx = quad.GetShapeFunctionDerivatives([xi, eta])

        for a in range(4):
            for b in range(4):
                for i in range(2):
                    ref_lhs[a,b] += gauss_weight * det_jacobian * dN_dx[i,a] * dN_dx[i, b]
        
    element_lhs = element.CalculateLeftHandSideMatrix()
    print(f"ref_lhs after: \n {ref_lhs}")
    print(f"elemment_lhs: \n {element_lhs}")
    if (np.linalg.norm(ref_lhs - element_lhs) > 1e-10):
        raise RuntimeError(f"Reference values and lhs do not match. \n{ref_lhs}\n{element_lhs}")
    else:
        print("Unit test passes")


    