# TODO: Entity object with iter overload so you can get its connections
# TODO: Make an field with id instead of elements and maybe negative for connections
from base_structures import Pos, print_field, abrir_arq_chain, CD_OFFS
from auxiliary import *

# - Connection -
class Connection:
    # - Class related -
    def __init__(self, id_from: int, id_to: int, dir: int, type: int):
        self.from_id = id_from
        self.to_id = id_to
        self.dir = dir
        self.type = type
    # - Features -
    def __eq__(self, value):
        return (
            (self.from_id == value.from_id) and
            (self.to_id == value.to_id) and
            (self.dir == value.dir) and
            (self.type == value.type)
        )
    def __str__(self):
        return f'(From: {self.from_id}, To: {self.to_id}, Dir: {self.dir}, Type: {self.type})'
    def __repr__(self):
        return self.__str__()

# - (Chemical) Entity -
class Entity:
    # - Class Related -
    def __init__(self, ent_id: int, el: str, cons: list[Connection]):
        # ~ Pos ~
        # self.pos: Pos = pos
        
        # Id
        self.id = ent_id
        self.el = el
        self.cons = cons
    
    def at_dir(self, dir_to: int):
        """
        Finds id (if any) to direction being compared
        
        :param dir_from: Direction of comparison

        Returns id if found or -1 if not
        """
        for con in self:
            if con.dir == dir_to:
                # Found
                return con.to_id
        # Not found
        return -1
    # - Visual -
    def __str__(self):
        return f"({self.id}: {self.el} - {self.cons})"
    def __repr__(self):
        return f"({self.id}: {self.el})"
    # - Features -
    def __eq__(self, value):
        if isinstance(value, int):
            return self.id == value
        elif isinstance(value, str):
            return self.el == value
        else:
            if len(self.cons) == len(value.cons):
                val = (
                    (self.id == value.id) and
                    (self.el == value.el) and
                    (all([self.cons[i] == value.cons[i] for i in range(len(self.cons))]))
                )
                return val
            else:
                return False
    def __iter__(self):
        for i in range(len(self.cons)):
            yield self.cons[i]
    # Maybe accessor function to chain.chain?
    def __len__(self):
        return len(self.cons)

# - Classes -
# Single chain (only one path) representation
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
        # Maybe instead of this, I could get alternatives after...
        # self.paths: list[int] = []
        # List of indexes for the main paths (name-definer)
        self.main_path: list[int] = []

        # ~ Test ~
        self.chain: list[Entity] = []
        self.main_chain: list[Entity] = []
        self.groups: list[list[Entity]] = []

        self.edges: list = []

        # Create chain representation of the field
        # It indexes self.chain by id_pos
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
            self.chain.append(Entity(id_pos, el_str, cons))
            if el_str == 'C' and c_count <= 1:
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
    
    # - Chemistry related -
    def get_main_path_id(self, pool_id):
        """
        Get main path id (starting from iC to fC following main path)
        
        :param self: Chain object
        :param pool_id: The id of the entity within the self (chain) id_pool
        """
        return self.main_path.index(pool_id) + 1

    # - Element-wise -
    # def add_path(self, path):
    #     """
    #     Adds a new path if it is unique
    #     """
    #     if not self.paths:
    #         self.paths.append(path)
    #         return
    #     for i_path in self.paths:
    #         if (len(i_path) == len(path)) and (sorted(i_path) != sorted(path)):
    #             self.paths.append(path)
    #             break
    def get_main_path_size(self):
        return len(self.main_path)
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
        return len(self.main_chain)
    def __getitem__(self, i: int):
         return self.chain[i]
