from geometries.node import Node
from geometries.geometry import Geometry
from geometries.quad2d4n import Quad2D4N
from geometries.triangle2d3n import Triangle2D3N
from geometries.line2d2n import Line2D2N
#from __future__ import annotations

class Mesh:
    def __init__(self, nodes_list: list[Node], geometries_list: list[Geometry]):
        self.nodes_list = nodes_list
        self.geometries_list = geometries_list
        self.dict_of_sub_mesh: 'dict[str, Mesh]' = {}
    def AddSubMesh(self, sub_mesh_name: str, mesh: 'Mesh') -> None:
        self.dict_of_sub_mesh[sub_mesh_name] = mesh

    def GetSubMesh(self, sub_mesh_name: str) -> 'Mesh':
        return self.dict_of_sub_mesh[sub_mesh_name]

    def __str__(self)-> str:
        msg = "Mesh info:\n"
        msg += f"Number of nodes: {len(self.nodes_list)}\n"
        msg += f"Number of geometries: {len(self.geometries_list)}"
        for sub_mesh_name, sube_mesh in self.dict_of_sub_mesh.items():
            msg += f"\n---- Sub mesh {sub_mesh_name}\n"
            msg += str(sube_mesh)
        return msg
    
def IsNotEmptyString(value: str) -> bool:
    return len(value) != 0


def IdentifyBlock(lines: 'list[str]', offset_index: int, starting_string: str, ending_string: str) -> tuple[int, int]:
    for i, line in enumerate(lines):
        if line.startswith(starting_string):
            start_index = i
            break
    end_index = lines.index(ending_string)
    return start_index + offset_index, end_index + offset_index



def ReadMesh(filename: str) -> Mesh:
    with open(filename, "r") as file_input:
       lines = file_input.readlines()
    
    node_start_index = lines.index("Begin Nodes\n")   #6            --- in reality at line 7
    node_end_index = lines.index("End Nodes\n")  #23                      --- python indexing start from 0. 
    #print(f"Start: {node_start_index + 1}, End: {node_end_index + 1}")



    list_of_nodes: list[Node] = []  #instantiating empty list to store nodes with id and coords
    #looping over the lines between begin adn end nodes: - 8 -23   {jjust -1 the lines in mpda file and you will get it}
    for node_data in lines[node_start_index + 1: node_end_index]:      
        raw_data = node_data.split(" ")
        data = list(filter(IsNotEmptyString, raw_data))

        
        node_id = int(data[0])
        node_x =  float(data[1])
        node_y =  float(data[2])
        node_z =  float(data[3])
        node = Node(node_id, node_x, node_y, node_z)                     #need to add values "temperatrue" in order to prevent error in main..py
        list_of_nodes.append(node)
    print(f"----- Read Nodes: {len(list_of_nodes)}")
      

    element_start_index, element_end_index = IdentifyBlock(lines[node_end_index + 1: ],
                                                           node_end_index + 1,
                                                           "Begin Elements",
                                                           "End Elements\n")
    #print(f"Ele start: {element_start_index}, End: {element_end_index}")
    
    list_of_geometries: 'list[Geometry]' = []
    for ele_data in lines[element_start_index + 1: element_end_index]:
        raw_data = ele_data.split(" ")
        data = list(filter(IsNotEmptyString, raw_data))
        #print(data)
        element_id = int(data[0])
        nodal_connectivites: 'list[Node]' = []
        for node_id_str in data[2: len(data)-1]:
            node_id = int(node_id_str)
            node_id_not_found = True
            for node in list_of_nodes:
                if node_id == node.node_id:
                    nodal_connectivites.append(node)
                    node_id_not_found = False
                    break
        if node_id_not_found is True:
            raise RuntimeError(f"node id: {node_id} not found")
        
        if len(nodal_connectivites) == 3:
            geometry = Triangle2D3N(element_id, nodal_connectivites)
        
        elif len(nodal_connectivites) == 4:
            geometry = Quad2D4N(element_id, nodal_connectivites)
        else:
            raise RuntimeError(f"Unsupported Geometry")
        list_of_geometries.append(geometry)

    # Identifying the starting index of conditions
    condition_start_index, condition_ending_index = IdentifyBlock(
                        lines[element_end_index + 1:],
                        element_end_index + 1,
                        "Begin Conditions",
                        "End Conditions\n")

    for condition_line in lines[condition_start_index+1: condition_ending_index]:
        raw_data = condition_line.split(" ")
        data = list(filter(IsNotEmptyString, raw_data))

        condition_id = int(data[0])
        nodal_connectivities: 'list[Node]' = []
        for node_id_str in data[2:]:
            node_id = int(node_id_str)
            node_id_not_found = True
            for node in list_of_nodes:
                if node_id == node.node_id:
                    nodal_connectivities.append(node)
                    node_id_not_found = False
                    break

            if node_id_not_found:
                raise RuntimeError(f"Node id {node_id} not found.")

        if len(nodal_connectivities) == 2:
            geometry = Line2D2N(condition_id, nodal_connectivities)
        else:
            raise RuntimeError(f"Unsupported geometry.")

        list_of_geometries.append(geometry)

    print(f"---- Read Geometries: {len(list_of_geometries)}")
    mesh = Mesh(list_of_nodes, list_of_geometries)        

    #taking element_end_index as offset_index 
    start_index = element_end_index

    while(start_index < len(lines)):
        #Algorithm for Bottom submesh nodes
        SubModelPart_start_index, SubModelPart_end_index = IdentifyBlock(lines[start_index:], start_index, "Begin SubModelPart", "End SubModelPart\n")
        #print(f"Submesh begin: {SubModelPart_start_index}, End: {SubModelPart_end_index}")
        submesh_name = lines[SubModelPart_start_index].split(" ")[2]
        #print(f"submesh name: {submesh_name}")
        submesh_node_begin, submesh_node_end = IdentifyBlock(lines[SubModelPart_start_index: SubModelPart_end_index ], SubModelPart_start_index, "    Begin SubModelPartNodes", "    End SubModelPartNodes\n")
        #print(f"Node begin line: {submesh_node_begin}, end line: {submesh_node_end}")
        list_of_submesh_nodes: 'list[Node]' = []
        for submesh_node_str in lines[submesh_node_begin + 1: submesh_node_end]:
            #cleaned_data = [item.strip() for item in data]  ----Not Needed {Wanted to remove \n}
            submesh_node = int(submesh_node_str)
            #print(f"submesh nodes: {submesh_node}")
            node_id_not_found = True
            for node in list_of_nodes:
                if submesh_node == node.node_id:
                    found_node = node
                    node_id_not_found = False
                    break
            if node_id_not_found:
                raise RuntimeError(f"Node id {submesh_node} node found. ")
            list_of_submesh_nodes.append(found_node)
        
        #Algortihm for bottom submesh elements
        submesh_element_begin, submesh_element_end = IdentifyBlock(lines[submesh_node_begin: ], submesh_node_begin, "    Begin SubModelPartElements", "    End SubModelPartElements\n")
        #print(f"Submesh ele begin: {submesh_element_begin}, end: {submesh_element_end}")
        list_of_submesh_geometries: 'list[Geometry]' = []
        for submesh_ele_str in lines[submesh_element_begin + 1: submesh_element_end]:
            submesh_ele = int(submesh_ele_str)
            submesh_ele_not_found = True
            for element in list_of_geometries:
                if submesh_ele == element.geometry_id:
                    found_element = element
                    submesh_ele_not_found = False
                    break
            if submesh_ele_not_found:
                raise RuntimeError(f"Element id {submesh_ele} not found")
            list_of_submesh_geometries.append(found_element)
        
        submesh = Mesh(list_of_submesh_nodes, list_of_submesh_geometries)
        mesh.AddSubMesh(submesh_name, submesh)
        start_index = SubModelPart_end_index + 1
        
    return mesh


if __name__ == "__main__":
    #mesh = ReadMesh("tests/mesh_io/example_mesh_1.mdpa")
    mesh =  ReadMesh("Task1_QuadrilateralMeshType_2D.mdpa")
    print(mesh)


