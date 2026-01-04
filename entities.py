# TODO: Entity object with iter overload so you can get its connections
# TODO: Make an field with id instead of elements and maybe negative for connections
from base_structures import Pos, print_field, open_f, CD_OFFS
from constants import *
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
        if isinstance(value, Connection):
            # Checks equality between values for the connection
            return (
                (self.from_id == value.from_id) and
                (self.to_id == value.to_id) and
                (self.dir == value.dir) and
                (self.type == value.type)
            )
        elif isinstance(value, int):
            # Checks type of connection
            return self.type == value
    def __lt__(self, value):
        if isinstance(value, int):
            return self.type < value
    def __gt__(self, value):
        if isinstance(value, int):
            return self.type > value
    def __le__(self, value):
        if isinstance(value, int):
            return self.type <= value
    def __ge__(self, value):
        if isinstance(value, int):
            return self.type >= value
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

        # Classification
        # (Relative to number of carbons connected to it)
        self.classif = ''
    
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
    # - Chemistry -
    def is_hetero(self) -> bool:
        return self.el in HETEROATOMS
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
                    (self.id == value.id) # and
                    # (self.el == value.el) and
                    # (all([self.cons[i] == value.cons[i] for i in range(len(self.cons))]))
                    # Commented because all those calclations are somewhat redundant!
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
# TODO: > and < should compare size between chains!
class Chain:
    # Make it so its pssible to get pos from id
    def __init__(self, builder: str | list = None):
        """
        Paths are represented by an List of integer IDs
        builder: field array or path string
        """
        if isinstance(builder, str):
            self.load_file(builder)
        elif isinstance(builder, tuple) or isinstance(builder, list):
            self.initiate(builder)
        elif builder is None: # Empty chain 
            pass
        else:
            raise TypeError('Builder is neither an file path (str) or field\'s list!')
    
    # - Class related -
    # TODO: Refactor this whole mess
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

        self.elements = {'H': 0}

        self._make_chain_from_field(field)

        # Classification of the elements (based on number of connected carbons)
        for el in self.chain:
            el.classif = CCLASSIFICATION[self.to_els(el.id, el.el)]
        # - Elements info -
        # Adds as element
        for el in self.chain:
            # Molecular formula
            self._add_element(el)
        
        # - Chemical properties -
        self.name = ""
        self.functional = ""
        self.func_name = ""
        self.mol_formula = ""
        
        for el, num in self.elements.items():
            if num:
                self.mol_formula += f"{el}" +\
                    f"\033[38;5;232m{str(num)}\033[0m"

    def load_file(self, filename: str):
        field: tuple[tuple[str]] = open_f(filename)
        self.initiate(field)

    def load_chain(self, chain: list[Entity]):
        field: tuple[tuple[str]] = self._make_field(chain)
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
    
    def _make_field(self, chain: list[Entity]):
        """
        Make field from a list of entities.
        
        :param self: Chain object
        :param chain: List of entities
        :type chain: list[Entity]
        """
        # - Variables -
        field: list[str] = []
        coords: list[Pos] = []
        els: list[str] = []
        
        curr_pos: Pos = Pos(0, 0)

        # (ID, Entity, Direction, Current Position)
        action_stack: list[tuple] = []
        
        if chain:
            curr_ent_id: int = 0
            curr_entity: Entity = chain[curr_ent_id]
            curr_dir: int = 0
            curr_type: int = 0
            action_stack.append((curr_ent_id, curr_entity, curr_dir, curr_pos, curr_type))

            # First one is the starting atom
            coords.append(curr_pos)
            els.append(curr_entity.el)

        # 1. Give a 2D coordinate for each entity ((0,0) for the first one)
        while action_stack:
            # Descompacts current action info
            curr_ent_id, curr_entity, curr_dir, curr_pos, curr_type = action_stack.pop()
            
            # Get current position (of connection) based on direction
            if curr_dir != 0:
                curr_pos += CD_OFFS[curr_dir]

                # Connection
                coords.append(curr_pos)
                els.append(curr_type)

                # Only if the is curr_type
                curr_pos += CD_OFFS[curr_dir]

                # Atom
                coords.append(curr_pos)
                els.append(curr_entity.el)

            # Adds next step(s) to action stack
            for con in curr_entity:
                if con.dir != -curr_dir:
                    action_stack.append((curr_entity.id,
                                         chain[con.to_id],
                                         con.dir,
                                         curr_pos,
                                         con.type))

        
        # 2. Use coords to fill field out
        # First, need to find maximum and minium in each axis
        min_row = min(coords, key=lambda coord: coord.row).row
        min_col = min(coords, key=lambda coord: coord.col).col

        max_row = max(coords, key=lambda coord: coord.row).row
        max_col = max(coords, key=lambda coord: coord.col).col

        field = [['0' for _ in range(min_col, max_col + 1)] for _ in range(min_row, max_row + 1)]

        for row in range(min_row, max_row + 1):
            for col in range(min_col, max_col + 1):
                curr_pos = Pos(row, col)
                normalized_pos = curr_pos + Pos(abs(min_row), abs(min_col))
                if curr_pos in coords:
                    curr_pos_id = coords.index(curr_pos)
                    field[normalized_pos.row][normalized_pos.col] = str(els[curr_pos_id])
        return field

    def _make_chain_from_field(self,  field: list[str]):
        # Create chain representation of the field
        # It indexes self.chain by id_pos
        _nitrogen_chain = False
        _insasts = set() # for optimization get all insats ahead
        for id_pos, pos in enumerate(self.id_pool):
            el_str = field[pos.row][pos.col]
            con_info = scout(field, self.id_pool, id_pos, nxt_dir=True, nxt_con_val=True, nxt_str=True)
            # Connections
            cons = []
            c_count = 0
            all_count = 0
            for raw_con in con_info:
                nxt_id, dir, type, ent_str = raw_con
                if (type > 1) and (ent_str == 'C'):
                    _insasts.add(id_pos)
                cons.append(Connection(id_pos, nxt_id, dir, type))
                if ent_str == 'C':
                    c_count += 1
                all_count += 1
            self.chain.append(Entity(id_pos, el_str, cons))

            # If is Carbon and with one connection
            # TODO: Handle Nitrogen (Maybe treat them as 'C' if detected)
            # Depends on the situation though!
            if (el_str == 'C' and c_count <= 1) and (not _nitrogen_chain):
                self.edges.append(id_pos)
            # If there is a Nitrogen with ONLY connections to Cs, turns into
            # An edge (amine/amide)
            if (el_str == 'N' and c_count == all_count):
                self.edges = [id_pos]
                _nitrogen_chain = True
    
        # For now, if edgless - and there is more than 1 carbon -
        # Treat insaturations (in _insats) as an edge
        # Otherwise the first position is the "edge"
        # NOTE: This make it possible for it to start and end
        #       whithin a cycle "between" two edges...
        #       though it appear to not cause a problem with pathfinding
        if self.edges == [] and len(self.id_pool) > 0:
            # There will at least be 2 "insaturations" for each double C connections
            if _insasts:
                self.edges.extend(sorted(_insasts))
            else:
                self.edges.append(0)

    def set_main_path(self, main_path: list[int]):
        if main_path:
            self.main_path = main_path
            print("Main chain:")
            print_field(self.field,
                        [self.id_pool[id] for id in self.main_path],
                        show_ids=True)
            self.main_chain = []
        else:
            # No need to raise an error I suppose!
            print("Empty field!")
        
        # Main chain
        for pos_id in self.main_path:
            self.main_chain.append(self.chain[pos_id])
    
    def copy(self):
        """
        Returns copy of itself
        
        :param self: Description
        :return: Chain copy
        :rtype: Chain
        """
        cchain = Chain()
        # TODO: Complete this
        return cchain
    
    # - Chemistry related -
    def get_main_path_id(self, path_id: int):
        if path_id in self.main_path:
            return self.main_path.index(path_id) + 1
        else:
            return -1

    def get_main_path_ids(self, path_ids: list[int]):
        for id in path_ids:
            if id in self.main_path:
                yield self.get_main_path_id(id)
            else:
                yield -1

    # - Element-wise -
    def get_main_path_size(self):
        return len(self.main_path)
    
    def to_els(self, id: int, el: str = '', filter= lambda e: True) -> int:
        """
        Docstring for to_el
        
        :param self: Chain object
        :param id: Id of the atom to 'look around'
        :type id: int
        :param el: Element string to restrict to
        :type el: str
        :param filter: Filter function applied to each Connection
        :type el: function
        :return: Quantity of connections that satisfied the el str and filter function
        :rtype: int
        """
        if el:
            return sum([self.chain[con.to_id].el == el for con in self.chain[id] if filter(con)])
        else:
            return len([con for con in self.chain[id] if filter(con)])
    
    def get_to_els(self, id: int, filter = lambda e: True) -> list[Entity]:
        return [self.chain[con.to_id] for con in self.chain[id] if filter(con.to_id)]

    def _add_element(self, el: Entity):
        el_str = el.el
        # The element itself
        if el_str in self.elements:
            self.elements[el_str] += 1
        else:
            self.elements[el_str] = 1

        # Hydrogen bonds (if any)
        self.elements['H'] += 4 - len(el.cons)

    # - Visual -
    def print_ids(self):
        mock_field = [[0 for _ in range(self.n_col)] for __ in range(self.n_row)]
        for i in self.id_pool:
            pos = self.id_pool[i]
            mock_field[pos.row][pos.col] = i
        print_field(mock_field)
    def __str__(self):
        summary = f"Name: {self.name}\n"+\
        f"Functional class: {self.functional}\n"+\
        f"Molecular formula: {self.mol_formula}"
        return summary
    def __repr__(self):
        return self.__str__()
    
    # - Functionalities -
    def __len__(self):
        return len(self.main_chain)
    def __getitem__(self, i: int):
         return self.chain[i]
