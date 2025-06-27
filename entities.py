# TODO: Change representation to be individual know atoms (include their valence?)
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

"""
Representative class of a single conection (except when is to a Hidrogen)
"""
class Conexao:
    def __init__(self, tipo: int, id_el: int):
        self.tipo = tipo
        self.id_el = id_el
"""
Generalized chemical entity representation
"""
class Entity:
    id = 0
    def __init__(self, z:int, *conexs: tuple[Conexao]):
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
    
    # Representations
    def scheme_chain(self):
        return "Chain"
    def __repr__(self):
        return self.scheme_chain


cTest = Entity([0])
print(cTest)