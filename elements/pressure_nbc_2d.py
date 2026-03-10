from elements.element import Element
import typing
import numpy as np
from geometries.geometry import Geometry

class PressureNBC2D(Element):
    def __init__(self, element_id: int, geometry: Geometry, parameters: 'dict[str, typing.Any]' ) -> None:
        super().__init__(element_id, geometry, parameters)
        self.pressure = parameters["pressure"]
        self.normal_direction = np.array([
            geometry.list_of_nodes[1].y - geometry.list_of_nodes[0].y,
            geometry.list_of_nodes[0].x - geometry.list_of_nodes[1].x
        ])
        self.normal_direction = self.normal_direction / np.linalg.norm(self.normal_direction)
        """
        This reflects line loads


        """
    def CalculateLeftHandSideMatrix(self) -> np.ndarray:
        number_of_nodes = len(self.geometry.list_of_nodes)
        np.zeros(self.geometry.list_of_nodes * 2, self.geometry.list_of_nodes * 2)
        return np.zeros((number_of_nodes * 2, number_of_nodes * 2))
    

    def CalculateRightHandSideVector(self) -> np.ndarray:
        number_of_nodes = len(self.geometry.list_of_nodes)
        rhs = np.zeros((number_of_nodes * 2))
        traction = self.normal_direction * self.pressure


        for gauss_weight, (xi, eta) in self.geometry.GetGaussIntegrationPoints(2):
            det_jacobian = self.geometry.GetJacobianDeterminant([xi, eta])
            shape_functions = self.geometry.GetShapeFunctionValue([xi, eta])

            for a in range(number_of_nodes): #N_a
                for i in range(2): #t_i
                    rhs[a * 2 + i] += shape_functions[a] * traction[i] * det_jacobian * gauss_weight
        return rhs

    def GetDofIdsVector(self) -> list[int]:
        dof_ids: 'list[int]' = []
        for node in self.geometry.list_of_nodes:
            dof_ids.append(node.dof_ids["displacement_x"])
            dof_ids.append(node.dof_ids["displacement_y"])
        return dof_ids
    
    def GetDofNames(self) -> 'list[str]':
        return ["displacement_x", "displacement_y"]