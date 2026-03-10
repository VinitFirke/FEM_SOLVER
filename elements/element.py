import abc
from geometries.geometry import Geometry
from typing import Any #scalars, vectors can be anything
import numpy as np

class Element(abc.ABC):
    def __init__(self, element_id: int, geometry: Geometry, parameters: 'dict[str, Any]' ):
        """
        Construction of an element
        Params dictionary can hold both material parameters as well as body force values. 
        (Whatever is requuiredby the eleemtn to fromulate weak form of the PDE.)
        Args:
            element_id
            geometry: Geometry represented by the element parameters (dict[str, Any]): Parameters [Material adn others]
        """
        #Purpose of an element:  Calculate the Stiffness Matrix
        self.element_id = element_id
        self.geometry = geometry
        self.parameters = parameters

    @abc.abstractmethod
    def CalculateLeftHandSideMatrix(self) -> np.ndarray:
        """Calculate the left hand side matrix (stiffness matrix)

        Returns:
            numpy.ndarray: The computed square stiffness matrix
        """    
        pass

    
    @abc.abstractmethod
    def CalculateRightHandSideVector(self) -> np.ndarray:
        """Calculates the right hand side vector (Load Vector)

        Returns:
            numpy.ndarray: The computed right hand side vector
        """
        pass

    @abc.abstractmethod
    def GetDofIdsVector(self) -> 'list[int]':
        """Gives the dof id list which represents the ordering of components in left hand side matrix and right hand side vector.

        Returns:
            list[int]: List of integers representing the order.
        """
        pass

    @abc.abstractmethod
    def GetDofNames() -> 'list[str]':
        """
        Returns list of scalar dofs used in the element formulation
        Returns:
            list[str]: list of scalar dofs used in the element formulation
        """
        pass


    def __str__(self):
        msg = f"-- Element id: {self.element_id}"
        msg += f"-- node_ids: [" + ",".join([str(node.node_id) for node in self.geometry.list_of_nodes])

