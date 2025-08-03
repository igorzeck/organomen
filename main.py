# Test file for the chain separation algorithm
# -- Definition --
# Stubs are parts of a path

# A path is a collection of carbon starting
# and ending with a one connection carbon
# or with the same 2 carbons on the path

# A chain is made of all unique combinations
# of paths, which in turn define how
# the chain is named

# TODO: Make part division for subnaming
# TODO: make unique hold only unique paths
# TODO: Chain object can be classified as it is read
# TODO: Organize in different files and make the integration with Chain class better
# -- Imports --
from  copy import deepcopy

# -- Variables --
# - Classes -
class Pos:
    """
    2D Vector coordinate
    """
    def __init__(self, row, col):
        self.row = row
        self.col = col
    # Math
    def __add__(self, value):
        return Pos(self.row + value.row, self.col + value.col)
    def __sub__(self, value):
        return Pos(self.row - value.row, self.col - value.col)
    def __mul__(self, value: int):
        # Scalar multiplication
        return Pos(self.row * value, self.col * value)
    def __truediv__(self, value:int):
        return Pos(int(self.row / value), int(self.col / value))
    # Representations
    def __str__(self):
        return f"({self.row}; {self.col})"
    def __repr__(self):
        return self.__str__()
    # Relationals
    def __eq__(self, value):
        return (self.row == value.row and self.col == value.col)


class Chain:
    def __init__(self, field: list):
        """
        Paths are represented by an List of integer IDs
        field: 2D flat representation of the charbon chain
        """
        # - Field -
        self.field = field
        self.n_row = len(self.field)
        self.n_col = len(self.field[0])

        self.id_dict: dict[Pos] = self._make_id_dict()
        # List all of paths
        self.paths = []
        # List of indexes for the main paths (name-definer)
        self.main_path = []

        # - Chemical properties -
        self.name = ""
        self.func_name = ""
    
    def _make_id_dict(self):
        pos_dict = {}
        id_count = 1
        for row in range(self.n_row):
            for col in range(self.n_col):
                if not self.field[row][col].isnumeric():
                    pos_dict[id_count] = Pos(row, col)
                    id_count += 1
        return pos_dict
    
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
        for i in self.id_dict:
            pos = self.id_dict[i]
            mock_field[pos.row][pos.col] = i
        print_field(mock_field)
    def __str__(self):
        summary = f"Name: {self.name}"
        return summary

    def __repr__(self):
        return self.__str__()
# - Constants -
# Connections directions
NORTHWEST = -4
NORTH = -3
NORTHEAST = -2
WEST = -1
EAST = 1
SOUTHWEST = 2
SOUTH = 3
SOUTHEAST = 4

# Coordinates associated with each connection
cd_offsets = {
    -4:Pos(-1, -1),
    -3:Pos(-1, 0),
    -2:Pos(-1, 1),
    -1:Pos(0, -1),
    1:Pos(0, 1),
    2:Pos(1, -1),
    3:Pos(1, 0),
    4:Pos(1, 1),
}

# -- Functions --
# - File handling -
def abrir_arq_chain(nome_arq):
    """
    Handles .chain file oppening

    return: 2D matrix
    """
    with open(nome_arq, "r") as arq:
        field = [linha.split() for linha in arq.readlines()]
        return field


# - Visual -
def print_field(field: list[list], highlights: tuple = [Pos(-1, -1)], zeros_repr='.'):
    """
    Prints the 2D matrix representation

        of a given chain field
    """
    # TODO: Make it so it can have a highlight feature for a specific coordinate
    print("")
    for i_row, row in enumerate(field):
        for i_col, el in enumerate(row):
            # Converts el if needed
            if isinstance(el, int):
                el = str(el)
            curr_pos = Pos(i_row, i_col)
            if el == '0':
                p_el = zeros_repr.ljust(3)
            elif curr_pos in highlights:
                if curr_pos == highlights[0]:
                    # Middle
                    el = 'i' + el
                elif curr_pos == highlights[-1]:
                    # Start
                    el = 'f' + el
                p_el = "\033[34m" + el.ljust(3) + "\033[0m"
            else: 
                p_el = el.ljust(3)
            print(p_el, end = "")
        print("")
    print("")


# - Pathfinder -
# TODO: Change to only detect Carbon to Carbon connections on chain end!
def _edge_c(field, ids):
    """
    Gives the postion of the of an edge Carbon
    field: Flat representation of the chain
    ids: Dictionary with carbon id and its postion
    """
    curr_id = None
    edgeless = True
    for id in ids:
        pos = ids[id]
        el = field[pos.row][pos.col]
        if el == 'C':
            curr_id = id
            if len(_scout(field, ids, id)) == 1:
                if edgeless:
                    edgeless = False
                yield curr_id
    if curr_id and edgeless:
        # If finds no edge, return last curr_pos
        yield curr_id


def _scout(field, ids: dict[Pos], pos_id: int):
    nxt_el_l = []
    for dir in cd_offsets:
        pos = ids[pos_id]

        # Next direction and connectin (to scout)
        nxt_step = cd_offsets[dir]
        nxt_con_pos = pos + nxt_step
        nxt_el_pos = nxt_con_pos + nxt_step
        if _is_inside(nxt_el_pos, len(field), len(field[0])):
            # If is inside, check to see if there is a Carbon in the next path
            nxt_el = field[nxt_el_pos.row][nxt_el_pos.col]
            next_con = field[nxt_con_pos.row][nxt_con_pos.col]
            if next_con.isnumeric() and next_con != '0'\
                and\
                nxt_el == 'C':
                # Kinda dumb... But it wil work...
                selected_id = -1
                for id in ids:
                    _pos = ids[id]
                    if _pos == nxt_el_pos:
                        selected_id = id
                if selected_id < 0:
                    raise ValueError(f"Didn't find position {nxt_el_pos} id's")
                # print("Scout:",dir, nxt_el_pos)
                nxt_el_l.append((dir, selected_id))       
    return nxt_el_l


def _is_inside(dir: Pos, n_row, n_col):
    if (dir.row >= n_row or dir.col >= n_col)\
        or\
        (dir.row < 0 or dir.col < 0):
        return False
    return True


def _is_higher(best: list, contender: list):
    """
    Checks if path is best fit to be the main one.
    """
    if len(contender) > len(best):
        return True
    else:
        return False


# TODO: enroll all variables in a chain class
def per_path(chain: Chain,
             start_pos_id: int,
             path_stub: list,
             best_stub: list,
             origin_dir = 0,
             recur = 0):
    """
    Go through the 2D representation of a carbon chain
    path
    recur: Current recursion (TODO: take it out)
    origin_dir: Previous mirror direction to avoid backtracking
    """
    # Start position of this recursion
    pos_id = start_pos_id
    # Temporary as it doesn't work with recursion!
    while True:
        print_field(chain.field, [chain.id_dict[pos_id]])
        nxt_els = _scout(chain.field, chain.id_dict, pos_id)
        nxt_n_con = len(nxt_els)
        
        # Cycle detection
        if pos_id in path_stub:
            # It will have repeated position
            # But it will make it easier to glance at
            # Cyclic behaviour!
            path_stub.append(pos_id)
            return path_stub
        path_stub.append(pos_id)

        if nxt_n_con == 2 or pos_id == start_pos_id:
            # If it has one possible path, move
            jump = False
            for nxt_dir, nxt_el_pos_id in nxt_els:
                if nxt_dir != origin_dir:
                    origin_dir = -nxt_dir
                    pos_id = nxt_el_pos_id
                    jump = True
                    break
            if jump:
                continue
        if nxt_n_con > 2:
            # Recursion for multiple paths
            for nxt_dir, nxt_el_pos_id in nxt_els:
                if nxt_dir != origin_dir:
                    print(f"({recur + 1}) Going ({nxt_dir}, {nxt_el_pos_id})")
                    print(best_stub)
                    temp_stub = per_path(chain, nxt_el_pos_id, origin_dir=-nxt_dir, recur = recur + 1, path_stub = deepcopy(path_stub), best_stub=best_stub)
                    if _is_higher(best_stub, temp_stub):
                        best_stub = temp_stub
            break
        if nxt_n_con == 1:
            print(f"({recur}) Done")
            print(*path_stub)
            chain.add_path(path_stub)
            if _is_higher(best_stub, path_stub):
                best_stub = path_stub
                chain.main_path = best_stub
            print_field(chain.field, [chain.id_dict[id] for id in path_stub])
            break
    return best_stub


# TODO: Insert chain as a parameter
def per_chain(field):
    """
    Go through the entirety of a chain field representation
    """
    chain = Chain(field)
    for id in _edge_c(chain.field, chain.id_dict):
        print(f"Start from {chain.id_dict[id]}")
        chain.main_path = per_path(chain,
                                id,
                                path_stub=[],
                                best_stub=chain.main_path
                                )
    if chain.main_path:
        print("Longest:")
        print_field(field, [chain.id_dict[id] for id in chain.main_path])
    else:
        print("Empty field!") 
    # return unique_paths, bestest_path   
    return chain  

# -- Naming -- 
# - Main chain naming -
def _name_size_pref(n_main: int):
    # - Name prefix -
    prefixes = [
        'met',
        'et',
        'prop',
        'but',
        'pent',
        'hex',
        'hept',
        'oct',
        'non',
        'dec',
        'undec',
        'duodec'
    ]
    if n_main <= len(prefixes):
        return prefixes[n_main - 1]
    else:
        return None


def _name_con_type(cons: list):
    con_qte = [
        '',
        'adi',
        'atri'
    ]

    infix = ''

    # Indexes within chain
    i2 = [str(i + 1) for i, con in enumerate(cons) if con == '2']
    i3 = [str(i + 1) for i, con in enumerate(cons) if con == '3']
    # n2 = cons.count('2')
    # n3 = cons.count('3')
    n2 = len(i2)
    n3 = len(i3)

    # This way 'enin' can occur naturally
    if n2 != 0:
        if len(con_qte) >= n2:
            infix += con_qte[n2 - 1]
        infix += '-' + ','.join(i2) + '-' + 'en'
    if n3 != 0:
        if len(con_qte) >= n3:
            infix += con_qte[n3 - 1]
        infix += '-' + ','.join(i3) + '-' + 'in'
    if n2 == n3 == 0:
        infix += 'an'
    return infix

def name_chain(field, ids: dict[Pos], all_chains: list, main_chain: list):
    if not main_chain:
        return "Chain is empty!"
    prefix = _name_size_pref(len(main_chain))
    if not prefix:
        return "Chain is too big!"
    
    # Makes chain
    cons = []
    old_pos = None
    for pos in [ids[id] for id in main_chain]:
        if old_pos:
            dif_pos = pos - old_pos
            norm_dif = dif_pos/2
            con_pos = old_pos + norm_dif
            cons.append(field[con_pos.row][con_pos.col])
        old_pos = pos
    infix = _name_con_type(cons)
    
    return  prefix + infix + 'o'


def main():
    chain = "Chains/simple.chain"
    # Pathfinding
    field = abrir_arq_chain(chain)
    c_main = per_chain(field)
    c_main.name = name_chain(field, c_main.id_dict, c_main.paths, c_main.main_path)
    print(c_main)


# -- Start --
if __name__ == "__main__":
    main()