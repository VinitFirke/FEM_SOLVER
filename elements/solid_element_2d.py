import numpy
from typing import Any
from geometries.node import Node
from elements.element import Element
from geometries.quad2d4n import Quad2D4N
from geometries.geometry import Geometry, GeometryType

class SolidElement2D(Element):
    def __init__(self, element_id: int, geometry: Geometry, parameters: dict[str, Any]):
        super().__init__(element_id, geometry, parameters)
        if geometry.GetGeometryType() != GeometryType.Domain_integral:
            raise RuntimeError("SolidElement2D only supports domain integrals")

        self.nu = parameters["poisson_ratio"]
        self.E = parameters["youngs_modulus"]

    def CalculateLeftHandSideMatrix(self) -> numpy.ndarray:
        number_of_nodes = len(self.geometry.list_of_nodes)
        lhs = numpy.zeros((number_of_nodes * 2, number_of_nodes * 2))

        # The B matrix depends on the gauss point.
        #   dN/dx -> depends on the gauss point
        # 1. Construct C =
        # 2. Construct B =
        # 3. lhs = B^T@C@B

        ### Constructing C matrix
        C = numpy.array([
                [1.0,     self.nu, 0.0              ],
                [self.nu, 1.0,     0.0              ],
                [0.0,     0.0,     (1 - self.nu) / 2]
            ]) * self.E / (1 - self.nu ** 2)

        for gauss_weight, (xi, eta) in self.geometry.GetGaussIntegrationPoints(2):
            jacobian = self.geometry.GetJacobian([xi, eta])
            det_jacobian = numpy.linalg.det(jacobian)
            dN_dx = self.geometry.GetShapeFunctionDerivatives([xi, eta])

            B = numpy.zeros((3, number_of_nodes * 2))

            for i in range(number_of_nodes):
                B[0, i * 2]     = dN_dx[0, i]
                B[1, i * 2 + 1] = dN_dx[1, i]
                B[2, i * 2]     = dN_dx[1, i]
                B[2, i * 2 + 1] = dN_dx[0, i]

            lhs += gauss_weight * det_jacobian * B.T @ C @ B

        return lhs

    def CalculateRightHandSideVector(self) -> numpy.ndarray:
        rhs = numpy.zeros(len(self.geometry.list_of_nodes) * 2)
        return rhs

    def GetDofIdsVector(self) -> 'list[int]':
        dof_ids: 'list[int]' = []
        for node in self.geometry.list_of_nodes:
            dof_ids.append(node.dof_ids["displacement_x"])
            dof_ids.append(node.dof_ids["displacement_y"])
        return dof_ids

    def GetDofNames(self) -> 'list[str]':
        return ["displacement_x", "displacement_y"]

if __name__ == "__main__":
    E = 1e+4
    nu = 0.2

    node_1 = Node(1,  0, 0, 0)
    node_2 = Node(2,  2, 2, 0)
    node_3 = Node(3,  1, 4, 0)
    node_4 = Node(4, -1, 3, 0)

    # using dummy values for stress computation.
    # these does not need to be the solved ones, you
    # can assume these are given to you by a some solution
    # to compute the stresses.
    node_1.values["displacement_x"] = 0.1
    node_1.values["displacement_y"] = 0.2
    node_2.values["displacement_x"] = 0.3
    node_2.values["displacement_y"] = 0.4
    node_3.values["displacement_x"] = 0.1
    node_3.values["displacement_y"] = 0.6
    node_4.values["displacement_x"] = -0.7
    node_4.values["displacement_y"] = 0.8

    geom = Quad2D4N(1, [node_1, node_2, node_3, node_4])
    n = len(geom.list_of_nodes)

    """
        rhs_ref[0] -> node 1 x direction eqn.
        rhs_ref[1] -> node 1 y direction eqn.
        rhs_ref[2] -> node 2 x direction eqn.
        rhs_ref[3] -> node 2 y direction eqn.
        rhs_ref[4] -> node 3 x direction eqn.
        rhs_ref[5] -> node 3 y direction eqn.

    """
    rhs_ref = numpy.zeros((n * 2))

    # we need to compute e_xx, e_xy, e_yy
    for gauss_weight, (xi, eta) in geom.GetGaussIntegrationPoints(2):
        jacobian = geom.GetJacobian([xi, eta])
        det_jacobian = numpy.linalg.det(jacobian)
        dN_dx = geom.GetShapeFunctionDerivatives([xi, eta])

        e_xx = 0.0 # du_dx
        e_yy = 0.0 # dv_dy
        e_xy = 0.0 # 0.5(du_dy + dv_dx)
        for i in range(n):
            e_xx += dN_dx[0, i] * geom.list_of_nodes[i].values["displacement_x"]
            e_yy += dN_dx[1, i] * geom.list_of_nodes[i].values["displacement_y"]
            e_xy += 0.5 * dN_dx[1, i] * geom.list_of_nodes[i].values["displacement_x"]
            e_xy += 0.5 * dN_dx[0, i] * geom.list_of_nodes[i].values["displacement_y"]

        # we cannot use the general sigma = lambda * trace(epsilon) * I + 2 * mu * epsilon
        # where sigma is the stress tensor and epsilon is the strain sensor here because
        # our element is formulated with the plane stress assumption. Therefore, we need
        # to compute the stresses here as well using the same plane stress assumption.
        # hence using the followings to compute the stresses.
        sigma_x = E / (1 - nu ** 2) * (e_xx + nu * e_yy)
        sigma_y = E / (1 - nu ** 2) * (e_yy + nu * e_xx)
        tau_xy = E / (1 + nu) * e_xy

        for i in range(n):
            rhs_ref[i * 2]     += (dN_dx[0, i] * sigma_x + dN_dx[1, i] * tau_xy ) * gauss_weight * det_jacobian
            rhs_ref[i * 2 + 1] += (dN_dx[0, i] * tau_xy  + dN_dx[1, i] * sigma_y) * gauss_weight * det_jacobian

    # now construct the solid element
    elem = SolidElement2D(1, geom, {"poisson_ratio": nu, "youngs_modulus": E})
    # calculate the lhs
    lhs = elem.CalculateLeftHandSideMatrix()
    u = numpy.zeros((2 * n))
    for i in range(n):
        u[i * 2]     = geom.list_of_nodes[i].values["displacement_x"]
        u[i * 2 + 1] = geom.list_of_nodes[i].values["displacement_y"]

    # calculate the residual
    rhs = lhs @ u

    # check
    if (numpy.linalg.norm(rhs - rhs_ref) > 1e-12):
        raise RuntimeError("Error: Solution mismatch")
    else:
        print("It works !")