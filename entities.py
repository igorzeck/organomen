# TODO: Change representation to be individual know atoms (include their valence?)
# TODO: Adicionar algum objeto que represente coordenadas 2D para simplficar os cálculos
# 
# Imports
# 
from enum import Enum
import numpy as np
#
# Classes and objects that represent chemical entities
#
# Abstração
el_nome = {
    1: "H",
    6: "C",
    7: "N",
    8: "O"
}

# Connections directions
NORTHWEST = 0
NORTH = 1
NORTHEAST = 2
WEST = 3
EAST = 4
SOUTHWEST = 5
SOUTH = 6
SOUTHEAST = 7

# Coordinates associated with each connection
cd_offsets = [
    (-1, -1),
    (0, -1),
    (1, -1),
    (-1, 0),
    (1, 0),
    (-1, 1),
    (0, 1),
    (1, 1),
]

ConType = Enum("ConType",[
    ("SIMPLE", 1),
    ("DOUBLE", 2),
    ("TRIPLE", 3),
    ("QUADRUPLE", 4)
])


class Conection:
    """
    Representative class of a single connection (except when is to a Hidrogen)

    Attirbutes:
        type: Type of connection (SIMPLE (1) up to QUADRUPLE (4))
        id_ent: Entity id to where the connection points to
        dir: Cardinal direction of the connection
    """
    def __init__(self, type: int, id_to: ConType, dir: int):
        self.type = type
        # Id within the context of the matrix where it came from
        self.id_to = id_to
        self.dir = dir
    #
    # Representation
    #
    def __str__(self):
        return f"Conection(Type: {self.type} To: {self.id_to} Dir: {self.dir})"
    def __repr__(self):
        return self.__str__()

class Entity:
    id = 0
    def __init__(self, z:int, conexs: tuple[Conection]):
        """
        Generalized chemical entity representation

        Attributes
            z: Atomic Number
            *conexs: Connections tuple
        """
        #
        # Identification
        #
        # ID within the context of the originating matrix
        self.id = Entity.id
        Entity.id += 1

        self.conexs = conexs
        #
        # Chemistry
        #
        self.z = z  # Atomic Number
    #
    # Representation
    #
    def __str__(self):
        return f"Entity(\n\tId: {self.id}\n\tZ: {self.z},\n\tConnections: {self.conexs}\n)"
    def __repr__(self):
        return self.__str__()


"""
Representation of a organic chain
"""
class Chain:
    id = 0
    def __init__(self, l_ents: list[Entity], chain_mat: str, name=""):
        #
        # Identification
        #
        self.id = Chain.id
        Chain.id += 1
        
        if not name:
            self.name = f"Chain_{self.id}"

        self.chain_mat = chain_mat

        # Elements
        self.els = l_ents

    # Assessors
    def __getitem__(self, item):
        return self.els[item]
    #
    # Representations
    #
    def parse(self, el_ant: Entity, dir_ant: int):
        # Starts by yield current recursion element
        yield el_ant
        el_atual = None
        tam_conexs = len(el_ant.conexs)
        # If it is equal to the previouse element
        # then that means that it come back to the start!
        while tam_conexs > 1 or el_atual != el_ant:
            for cd_i, _ in enumerate(cd_offsets):
                    curr_con = el_ant.conexs[cd_i] 
                    if curr_con != 0:
                        # Direction that it should ignore
                        # It's the one diametrically opposite to the anterior
                        if dir_ant != -1:
                            if abs(cd_i - dir_ant) == 4:
                                continue
                        if len(el_ant.conexs) > 2:
                            # Recursion
                            pass
                            # self.parse(el_ant, )
                        else:
                            # Iteration
                            for cd_i, _ in enumerate(cd_offsets):
                                curr_con = el_ant.conexs[cd_i] 
                                if curr_con != 0:
                                    # Direction that it should ignore
                                    # It's the one diametrically opposite to the anterior
                                    if abs(cd_i - dir_ant) == 4:
                                        continue
                                    el_atual = self[curr_con.id_to]
                                    dir_ant = cd_i
                                    # Yields element
                                    yield el_atual
        
    def scheme_chain(self):
        scheme_str = ""
        for el in self.parse(self[0],-1):
            print(el)
        for l in self.chain_mat:
            for el in l:
                if isinstance(el, int):
                    if el == 0:
                        scheme_str += "  "
                    else:
                        scheme_str += f" c{el}"
                else:
                    scheme_str += f" {el[1]}"
            scheme_str += "\n"
        # Following along the trail and registring coordinates
        return scheme_str

    def __str__(self):
        return self.scheme_chain()
    def __repr__(self):
        header = f"Chain(Name: {self.name}, Id: {self.id})"
        return header


#
#  File handling
#
"""
File parsing as done fromt Left -> Right and Up -> Bottom
Returns the a chain object.
"""
def open_chain_file(file_name: str) -> Chain:
    els = []
    lines_mat = []
    with open(file_name, "r") as f:
        # Matrix representation of the file
        lines = [line.split() for line in f.readlines()]

        # Gives ids for the representations
        curr_id = 0
        for l in lines:
            _line_mat = []
            for el in l:
                if el.isnumeric():
                    _line_mat.append(int(el))
                else:
                    _line_mat.append((curr_id, el))
                    curr_id += 1
            lines_mat.append(_line_mat)
        # Creates connections
        for i_l, l in enumerate(lines_mat):
            for i_h, el in enumerate(l):
                if isinstance(el, tuple):
                    cardc = _look_around(lines_mat, i_l, i_h)

                    # Instantiate entity object
                    ent = Entity(el, cardc)
                    print(ent)
                    
                    els.append(ent)
    return Chain(els, lines_mat)



"""
Looks "around" - i.e. all 8 directions - a specifc char inside a matrix.
"""
def _look_around(matrix: list[list], *pivot_coor: tuple[int]):
    # Cardinal directions surrounding pivot
    print(f"Looking around {(pivot_coor)}")
    card_dir = []
    for curr_cd_off in cd_offsets:
        # Next target coordinates
        nxt_coor = (pivot_coor[0] + curr_cd_off[0],
                    pivot_coor[1] + curr_cd_off[1])
        # Checks if the next cell is within the matrix bounds
        if _is_out(nxt_coor, len(matrix), len(matrix[0])):
            card_dir.append(0)
            continue
        
        el = matrix[nxt_coor[0]][nxt_coor[1]]

        if el == 0:
            card_dir.append(0)
        else:
            #
            # Appends a connection element
            #
            # Searches where the connection leads to:
            nxt_nxt_coor = (pivot_coor[0] + curr_cd_off[0] * 2,
                    pivot_coor[1] + curr_cd_off[1] * 2)

            nxt_nxt_el = matrix[nxt_nxt_coor[0]][nxt_nxt_coor[1]]
            conx = Conection(ConType.SIMPLE, nxt_nxt_el[0],curr_cd_off)
            
            card_dir.append(conx)

    # print("Cardinal:", card_dir)
    return tuple(card_dir)


# Verifies if coordinate is outisde a matrix bounding box
def _is_out(coor: tuple[int], height, length) -> bool:
    # Checks if is beyond bouding box or before bounding box
    return ((coor[0] >= length or coor[1] >= height)\
        or (coor[0] < 0 or coor[1] < 0))