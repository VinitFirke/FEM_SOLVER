from elements.element import Element
import typing
import numpy as np
from geometries.geometry import Geometry

class ForceNBC2D(Element):
    def __init__(self, element_id: int, geometry: Geometry, parameters: 'dict[str, typing.Any]' ) -> None:
        super().__init__(element_id, geometry, parameters)
        self.traction = parameters["traction"]
        """
        This reflects line loads


        """
    def CalculateLeftHandSideMatrix(self) -> np.ndarray:
        number_of_nodes = len(self.geometry.list_of_nodes)
        
        return np.zeros((number_of_nodes * 2, number_of_nodes * 2))
    

    def CalculateRightHandSideVector(self) -> np.ndarray:
        number_of_nodes = len(self.geometry.list_of_nodes)
        rhs = np.zeros((number_of_nodes * 2))

        for gauss_weight, local_coordinates in self.geometry.GetGaussIntegrationPoints(2):
            det_jacobian = self.geometry.GetJacobianDeterminant(local_coordinates)
            shape_functions = self.geometry.GetShapeFunctionValue(local_coordinates)

            for a in range(number_of_nodes): #N_a
                for i in range(2): #t_i
                    rhs[a * 2 + i] += shape_functions[a] * self.traction[i] * det_jacobian * gauss_weight
        return rhs

    def GetDofIdsVector(self) -> list[int]:
        dof_ids: 'list[int]' = []
        for node in self.geometry.list_of_nodes:
            dof_ids.append(node.dof_ids["displacement_x"])
            dof_ids.append(node.dof_ids["displacement_y"])
        return dof_ids
    
    def GetDofNames(self) -> 'list[str]':
        return ["displacement_x", "displacement_y"]