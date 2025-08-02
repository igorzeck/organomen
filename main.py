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
def print_field(field: list[list], highlights: tuple = (Pos(-1, -1)), zeros_repr='.'):
    """
    Prints the 2D matrix representation

        of a given chain field
    """
    # TODO: Make it so it can have a highlight feature for a specific coordinate
    print("")
    for i_row, row in enumerate(field):
        for i_col, el in enumerate(row):
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
def _edge_c(field):
    """
    Gives the postion of the of an edge Carbon
    """
    curr_pos = None
    edgeless = True
    for i_row, row in enumerate(field):
        for i_col, el in enumerate(row):
            if el == 'C':
                curr_pos = Pos(i_row, i_col)
                if len(_scout(field, curr_pos)) == 1:
                    if edgeless:
                        edgeless = False
                    yield curr_pos
    if curr_pos and edgeless:
        # If finds no edge, return last curr_pos
        yield curr_pos


def _scout(field, pos: Pos):
    nxt_el_l = []
    for dir in cd_offsets:
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
                # print("Scout:",dir, nxt_el_pos)
                nxt_el_l.append((dir, nxt_el_pos))       
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
def per_path(field: list[list],
             start_pos: Pos,
             path_stub: list,
             best_stub: list,
             unique_paths: list,
             origin_dir = 0,
             recur = 0):
    """
    Go through the 2D representation of a carbon chain
    path
    recur: Current recursion (TODO: take it out)
    origin_dir: Previous mirror direction to avoid backtracking
    """
    # Start position of this recursion
    pos = start_pos
    # Temporary as it doesn't work with recursion!
    while True:
        print_field(field, [pos])
        nxt_els = _scout(field, pos)
        nxt_n_con = len(nxt_els)
        
        # Cycle detection
        if pos in path_stub:
            # It will have repeated position
            # But it will make it easier to glance at
            # Cyclic behaviour!
            path_stub.append(pos)
            return path_stub
        path_stub.append(pos)

        if nxt_n_con == 2 or pos == start_pos:
            # If it has one possible path, move
            jump = False
            for nxt_dir, nxt_el_pos in nxt_els:
                if nxt_dir != origin_dir:
                    origin_dir = -nxt_dir
                    pos = nxt_el_pos
                    jump = True
                    break
            if jump:
                continue
        if nxt_n_con > 2:
            # Recursion for multiple paths
            for nxt_dir, nxt_el_pos in nxt_els:
                if nxt_dir != origin_dir:
                    print(f"({recur + 1}) Going ({nxt_dir}, {nxt_el_pos})")
                    temp_stub = per_path(field, nxt_el_pos, origin_dir=-nxt_dir, recur = recur + 1, path_stub = deepcopy(path_stub), best_stub=best_stub, unique_paths=unique_paths)
                    if _is_higher(best_stub, temp_stub):
                        best_stub = temp_stub
            break
        if nxt_n_con == 1:
            print(f"({recur}) Done")
            print(*path_stub)
            unique_paths.append(path_stub)
            if _is_higher(best_stub, path_stub):
                best_stub = path_stub
            print_field(field, path_stub)
            break
    return best_stub


def per_chain(field):
    """
    Go through the entirety of a chain field representation
    """
    bestest_path = []
    unique_paths = []
    for pos in _edge_c(field):
        print(f"Start from {pos}")
        bestest_path = per_path(field,
                                pos,
                                path_stub=[],
                                best_stub=bestest_path,
                                unique_paths=unique_paths
                                )
    if bestest_path:
        print("Longest:")
        print_field(field, bestest_path)
    else:
        print("Empty field!") 
    return unique_paths, bestest_path   

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

def name_chain(field, all_chains: list, main_chain: list):
    if not main_chain:
        return "Chain is empty!"
    prefix = _name_size_pref(len(main_chain))
    if not prefix:
        return "Chain is too big!"
    
    # Makes chain
    cons = []
    old_pos = None
    for pos in main_chain:
        if old_pos:
            dif_pos = pos - old_pos
            norm_dif = dif_pos/2
            con_pos = old_pos + norm_dif
            cons.append(field[con_pos.row][con_pos.col])
        old_pos = pos
    infix = _name_con_type(cons)
    
    return  prefix + infix + 'o'


def main():
    path = "Chains/simple.chain"
    # Pathfinding
    field = abrir_arq_chain(path)
    unique, best = per_chain(field)
    # Naming
    print(name_chain(field, unique, best))


# -- Start --
if __name__ == "__main__":
    main()