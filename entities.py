# TODO: Change representation to be individual know atoms (include their valence?)
# 
# Imports
# 
from enum import Enum
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
# TODO: Tirar essses enums...
Cd = Enum("Cd", [
    "NORTHWEST",
    "NORTH",
    "NORTHEAST",
    "WEST",
    "EAST",
    "SOUTHWEST",
    "SOUTH",
    "SOUTHEAST",
])

# Coordinates associated with each connection
Cd_val = [
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
    def __init__(self, type: int, id_el: ConType, dir: Cd):
        self.tipo = type
        self.id_el = id_el
        self.dir = dir
"""
Generalized chemical entity representation
"""
class Entity:
    id = 0
    """
    Attributes
        z: Atomic Number
        *conexs: Connections tuple
    """
    def __init__(self, z:int, *conexs: tuple[Conection]):
        #
        # Identification
        #
        self.id = id
        Entity.id += 1

        self.conexs = conexs
        #
        # Chemistry
        #
        self.z = z  # Atomic Number

"""
Representation of a organic chain
"""
class Chain:
    id = 0
    def __init__(self, l_ids: list[int], name=""):
        #
        # Identification
        #
        self.id = id
        Entity.id += 1
        
        if not name:
            self.name = f"Chain_{self.id}"

        # Elements
        self.els = l_ids

    # Assessors
    def __getitem__(self, item):
        return self.els[item]
    
    #
    # Representations
    #
    def scheme_chain(self):
        return "Chain"

    def __str__(self):
        return self.scheme_chain()
    def __repr__(self):
        return self.scheme_chain()


#
#  File handling
#
"""
File parsing as done fromt Left -> Right and Up -> Bottom
Returns the a chain object.
"""
def open_chain_file(file_name: str) -> Chain:
    els = []
    with open(file_name, "r") as f:
        # Matrix representation of the file
        lines = [line.split() for line in f.readlines()]
        lines_mat = []

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
                    
                    # Instantiate connections object
                    for c in cardc:
                        pass
                    # Instantiate entity object
                    # ent = Entity()
                    
                    els.append(el)


    return Chain(els)



"""
Looks "around" - i.e. all 8 directions - a specifc char inside a matrix.
"""
def _look_around(matrix: list[list], pivot_v: int, pivot_h: int):
    # Cardinal directions surrounding pivot
    print(f"Looking around {(pivot_h, pivot_v)}")
    card_dir = []
    n = len(matrix)

    curr_cd = Cd.NORTHWEST
    i_h_min = max(0, pivot_h - 1)
    i_v_min = max(0, pivot_v - 1)

    i_h = i_h_min
    i_v = i_v_min
    # TODO: Rodar todas as 8 células para poder contar corretamente as direções
    for l in matrix[max(0, pivot_v - 1):min(n, pivot_v + 2)]:
        for el in l[max(0, pivot_h - 1):min(n, pivot_h + 2)]:
            # Ignores pivot
            if i_v != pivot_v or i_h != pivot_h:
                print(f"{curr_cd.name} - coor {(i_h, i_v)} [{el}]")
                if el == 0:
                    card_dir.append(0)
                else:
                    #
                    # Appends a connection element
                    #
                    # Searches where the connection leads to:
                    next_h_inc = Cd_val[curr_cd.value - 1][0]
                    next_v_inc = Cd_val[curr_cd.value - 1][1]
                    print((next_h_inc + i_h, next_v_inc + i_v))

                    next_el = matrix[next_h_inc + i_h][next_v_inc + i_v]
                    conx = Conection(ConType.SIMPLE, int(next_el[0]),Cd(curr_cd.value - 1))
                    
                    card_dir.append(conx)
                # Isso é bem feio...
                if curr_cd.value + 1 < 9:
                    curr_cd = Cd(curr_cd.value + 1)
            
            i_h += 1
        i_h = i_h_min
        i_v += 1
    print("Cardinal:", card_dir)
    return card_dir