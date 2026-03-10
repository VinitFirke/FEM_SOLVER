
import numpy as np
from geometries.geometry import Geometry
from geometries.node import Node
from elements.element import Element
from elements.heat_conduction_element import HeatConductionElement
from geometries.triangle2d3n import Triangle2D3N

class SystemMatrixBuilder:
    def __init__(self, list_of_elements: 'list[Element]'):
        self.list_of_elements = list_of_elements

    def InitializeDofs(self):
        current_dof_id = 0
        for element in self.list_of_elements:
            dof_names = element.GetDofNames()  #list of dof: ["temperature"]
            for node in element.geometry.list_of_nodes:
                for dof_name in dof_names:
                    if dof_name not in node.dof_ids.keys():
                        node.dof_ids[dof_name] = current_dof_id
                        node.values[dof_name] = 0.0
                        current_dof_id += 1
        
        self.global_lhs = np.zeros((current_dof_id, current_dof_id))
        self.global_rhs = np.zeros((current_dof_id))
        print(f"global_rhs size: {len(self.global_rhs)}")



    def BuildLHSandRHS(self) -> 'tuple[np.ndarray, np.ndarray]':
        """
        Returns the LHS Matrix and RHS Vector

        Returns:
            tuple[np.ndarray,np.ndarray]: LHS matrix and RHS vector

        """

        for element in self.list_of_elements:
            local_lhs = element.CalculateLeftHandSideMatrix()
            local_rhs = element.CalculateRightHandSideVector()
            dofs_order = element.GetDofIdsVector()
            number_of_dofs = len(dofs_order)
            print(f"dofs_order: {dofs_order}, size: {len(dofs_order)}")
            print(f"local_rhs: {local_rhs}, size: {len(local_rhs)}")
            for a in range(number_of_dofs):
                g_a = dofs_order[a]
                if g_a >= len(self.global_rhs):
                    raise IndexError(f"Invalid index g_a: {g_a}, global_rhs size: {len(self.global_rhs)}")
                for b in range(number_of_dofs):
                    g_b = dofs_order[b]
                    self.global_lhs[g_a,g_b] += local_lhs[a,b]
                self.global_rhs[g_a] += local_rhs[a]
    
    def GetSolution(self, list_of_nodes: 'list[Node]', list_of_dofs: 'list[str]') -> np.ndarray:
        global_u = np.zeros((self.global_rhs.shape[0]))
        for node in list_of_nodes:
            for dof_str in list_of_dofs:
                if dof_str in node.dof_ids and dof_str in node.values:
                    global_u[node.dof_ids[dof_str]] = node.values[dof_str]
                else:
                    print(f"Warning: Node ID {node.node_id} does not have DOF '{dof_str}'")
        return global_u
    
    def UpdateSolution(self, du: np.ndarray, list_of_nodes: 'list[Node]', list_of_dofs: list[str]) -> None:
        for node in list_of_nodes:
            for dof_str in list_of_dofs:
                dof_index = node.dof_ids[dof_str]
                node.values[dof_str] += du[dof_index]



    def ApplyDirichletBoundaryConditions(self, list_of_nodes: 'list[Node]', list_of_dofs: 'list[str]') -> None:
        """
        The list of nodes contains all the nodes which needs to apply the DBC
        The list of dofs contains the dof on each noe which needs DBCs

        self.global_lhs = which is the system matrix
        self.global_rhs  whihc is the residual vector


        1. you need to iterate within the global lhs, and identify the rows and columns 
        which corresponds to dirichlet BC
        a. get dof id from  specific dof(dof_id = Node.dof_ids("temperature))
        b. make the lhs and rhs components accordingly.
        2. You make those components in the lhs, rhs yero except for the diagonal
        """
        # Iterate over all nodes and their corresponding degrees of freedom
        for node in list_of_nodes:
            for dof_str in list_of_dofs:
                # Retrieve the dof_id corresponding to the specified degree of freedom (e.g., "temperature")
                dof_id = node.dof_ids[dof_str]

                # Modify the global LHS matrix and RHS vector
                # Zero out the row in the LHS and RHS for this dof_id
                self.global_rhs[dof_id] = 0.0
                for col in range(self.global_lhs.shape[1]):
                    self.global_lhs[dof_id, col] = 0.0

                # Set the diagonal entry of the LHS to 1.0 for this dof_id
                self.global_lhs[dof_id, dof_id] = 1.0


                
if __name__ == "__main__":
    node_1 = Node(1,0,0,0)
    node_2 = Node(2,2,2,0)
    node_3 = Node(3,-1,3,0)
    node_4 = Node(4,0,1,0)
    
    triangle_1 = Triangle2D3N(1, [node_1, node_2, node_4])
    element_1 = HeatConductionElement(1, triangle_1, {})

    triangle_2 = Triangle2D3N(1,[node_2, node_3, node_4])
    element_2 = HeatConductionElement(1, triangle_2, {})

    system_matrix_builder = SystemMatrixBuilder([element_1, element_2])
    system_matrix_builder.InitializeDofs()

    list_of_all_nodes = [node_1, node_2, node_3, node_4]
    for node in list_of_all_nodes:
        print(f"Temperature dof id of node {node.node_id}- ", node.dof_ids["Temperature"])

    list_of_nodes_to_be_fixed = [node_3, node_4]

    system_matrix_builder.BuildLHSandRHS()
    print("Before applying DBC: \n", system_matrix_builder.global_lhs)
    system_matrix_builder.ApplyDirichletBoundaryConditions(list_of_nodes_to_be_fixed, ["Temperature"])
    print("After applying DBC: \n", system_matrix_builder.global_lhs)