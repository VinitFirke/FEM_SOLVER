class Node:  
    def __init__(self, node_id: int, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.node_id = node_id
        self.values: 'dict[str, float]' = {}
        self.dof_ids: 'dict[str, int]' = {}


    def __str__(self) -> str:
        msg = f" --Node id: {self.node_id} "
        msg+= f"\n  -- loc x: {self.x}"
        msg+= f"\n  -- loc y: {self.y}"
        msg+= f"\n  -- loc z: {self.z}"

        msg += "\nvalues: "
        for var, value in self.values.items():
            msg += f"\n    -- {var} = {value}"
        #msg+= f"\n  -- values: {self.values}"
        return msg
    
if __name__ == "__main__":
    node_1 = Node(1,0,0,0)
    node_1.values["Temperature"] = 100.0
    node_1.values["Disp_x"] = 2.0
    node_1.dof_ids["Temperature"] = 3

    print(node_1)


                        
