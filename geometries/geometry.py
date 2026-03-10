
import abc
from geometries.node import Node
import numpy as np
from enum import Enum

class GeometryType(Enum):
    Boundary_integral = 1
    Domain_integral = 2


class Geometry(abc.ABC):
    def __init__(self, geometry_id: int, list_of_nodes: 'list[Node]') -> None:
        self.list_of_nodes = list_of_nodes
        self.geometry_id = geometry_id

    @abc.abstractmethod
    def GetShapeFunctionValue(self, list_of_local_coordinates: 'list[float]') -> 'list[float]':
        """
        Get the shape function values for given xi, eta...
        Args:
            list_of_local_coordinates(list[float]): local coordinates (xi, eta,...)
        
        Returns:
            list[float]: Values of Shape functions evaluated at given (xi, eta,...)

        """
        pass
    @abc.abstractmethod
    def GetGaussIntegrationPoints(self, integration_order: int) -> 'list[tuple[float, list[float]]]':   #depending on integration order, the no. of points will be different
        
        """
            Gets the gauss inegration point data for specified integration order
            
            Args:
                integration_order (int): integration order

            Returns:
                list[tuple[float, list[float]]]: Returns a list of ( gauss weight, [xi,eta,...]) tuples

        """
        pass
        #return [
        #         (weight1, [xi_1, eta_1])           #why use tuple, intead of list? - Because we know definitely we need weight and coordinates and nothing else
        #         (weight2, [xi_2, eta_2])
        #         (weight3, [xi_3, eta_3])
        #]


    @abc.abstractmethod
    def GetJacobian(self, list_of_local_coordinates: 'list[float]') -> np.ndarray:
        """
        Get the jacobian matrix for given local coordinate system (xi, eta, ...)

        J = [
                 [dx/dxi, dx/deta],
                 [dy/dxi, dy/deta]     jacobain for 2D       
        ]

        Args:
            list_of_local_coordinates(list[float]): local coordinates(xi, eta,...)

        Returns:
            numpy.ndarray: Jacobian matrix of size(2,2) in 2D and (3,3) in 3D

        """
        pass
    @abc.abstractmethod
    def GetJacobianDeterminant(self, list_of_local_coordinates: 'list[float]') -> float:
        """
        Required to calculated B=Neumann Boudnary conditions
        """
        pass

    @abc.abstractmethod
    def GetShapeFunctionDerivatives(self, list_of_local_coordinates: 'list[float]') -> np.ndarray :
        """
        Get shape function derivatives w.r.t global (x,y) cordinates
        
        shape_function_derivates = [
                    [dN0/dx dN1/dx dN2/dx],
                    [dN0/dy dN1/dy dN2/dy]
        ]

        Args:
            list_of_local_coordinates(list[float]): local coordinates(xi,eta,...)
        
        Returns:
            numpy.ndarray: Shape function derivatives (2,3) matrix for 3 noded triangle 2d geometry
        """
        
        
        pass

    
        #output:  shape_function_derivates = [
        #        [dN0/dx, dN1/dx, dN2/dx],
        #        [dN0/dy, dN1/dy, dN2/dy],
        #        .....
        #    ]
    
    @abc.abstractmethod
    def GetGeometryType(self) -> GeometryType:
        """
        Returns Geometry type

        Returns:
                GeometryType: Geometry Type

        """
        pass




      