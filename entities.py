# TODO: Entity object with iter overload so you can get its connections
# TODO: Make an field with id instead of elements and maybe negative for connections
from base_structures import Pos, print_field, abrir_arq_chain, CD_OFFS
from pathfinder import *

# - Connection -
class Connection:
    # - Class related -
    def __init__(self, id_from: int, id_to: int, dir: int, type: int):
        self.from_id = id_from
        self.to_id = id_to
        self.dir = dir
        self.type = type
    def __str__(self):
        return f'(From: {self.from_id}, To: {self.to_id}, Dir: {self.dir}, Type: {self.type})'
    def __repr__(self):
        return self.__str__()

# - (Chemical) Entity -
class Entity:
    # - Class Related -
    def __init__(self, ent_id, el: str, cons: list[Connection]):
        # ~ Pos ~
        # self.pos: Pos = pos
        
        # Id
        self.id = ent_id
        self.el = el
        self.cons = cons
    # - Visual -
    def __str__(self):
        return f"({self.id}: {self.el} - {self.cons})"
    def __repr__(self):
        return f"({self.id}: {self.el})"
    # - Features -
    def __iter__(self):
        for i in range(len(self.cons)):
            yield self.cons[i]

# - Classes -
# Poderia guardar os carbonos em um gerador
# Que itera por referência eles
class Chain:
    def __init__(self, builder: str | list):
        """
        Paths are represented by an List of integer IDs
        builder: field array or path string
        """
        if isinstance(builder, str):
            self.load_file(builder)
        elif isinstance(builder, list):
            self.initiate(builder)
        else:
            raise TypeError('Builder is neither an file path (str) or field\'s list!')
    
    # - Class related -
    def initiate(self, field: list[str]):
        # - Field -
        self.field = field
        self.n_row = len(self.field)
        # Pre-supposing simmetry
        self.n_col = len(self.field[0])
        
        # I will keep as a tuple
        self.id_pool: tuple[Pos] = self._make_id_pool()

        # List all of paths
        self.paths: list[int] = []
        # List of indexes for the main paths (name-definer)
        self.main_path: int = []

        # ~ Test ~
        self.chain: list[Entity] = []
        self.edges: list = []

        # Create chain representation of the field
        for id_pos, pos in enumerate(self.id_pool):
            el_str = field[pos.row][pos.col]
            con_info = scout(field, self.id_pool, id_pos, nxt_dir=True, nxt_con_val=True, nxt_str=True)
            # Connections
            cons = []
            c_count = 0
            for raw_con in con_info:
                nxt_id, dir, type, ent_str = raw_con
                cons.append(Connection(id_pos, nxt_id, dir, type))
                if ent_str == 'C':
                    c_count += 1
            self.chain.append(Entity(pos, el_str, cons))
            if el_str == 'C' and c_count == 1:
                self.edges.append(id_pos)
        # - Chemical properties -
        self.name = ""
        self.functional = ""
        self.func_name = ""


    def load_file(self, arq: str):
        field: list[list[str]] = abrir_arq_chain(arq)
        # Maybe is better to return those two isnt it?
        self.initiate(field)
    
    # - Construct -
    def _make_id_pool(self):
        pos_pool = []
        for row in range(self.n_row):
            for col in range(self.n_col):
                if not self.field[row][col].isnumeric():
                    pos_pool.append(Pos(row, col))
        return tuple(pos_pool)
    
    # - Element-wise -
    def add_path(self, path):
        """
        Adds a new path if it is unique
        """
        if not self.paths:
            self.paths.append(path)
            return
        for i_path in self.paths:
            if (len(i_path) == len(path)) and (sorted(i_path) != sorted(path)):
                self.paths.append(path)
                break
    
    # - Visual -
    def print_ids(self):
        mock_field = [[0 for _ in range(self.n_col)] for __ in range(self.n_row)]
        for i in self.id_pool:
            pos = self.id_pool[i]
            mock_field[pos.row][pos.col] = i
        print_field(mock_field)
    def __str__(self):
        summary = f"Name: {self.name}\
                   Functional class: {self.functional}"
        return summary
    def __repr__(self):
        return self.__str__()
    
    # - Functionalities -
    def __len__(self):
        return len(self.paths)
