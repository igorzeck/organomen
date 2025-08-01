# Test file for the chain separation algorithm
# TODO: make the field a global variable???
# -- Imports --
from  copy import deepcopy

# -- Variables --
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
    def __mul__(self, value: int):
        # Scalar multiplication
        return Pos(self.row * value, self.col * value)
    # Representations
    def __str__(self):
        return f"({self.row}, {self.col})"
    def __repr__(self):
        return self.__str__
    # Relationals
    def __eq__(self, value):
        return (self.row == value.row and self.col == value.col)


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
    -3:Pos(0, -1),
    -2:Pos(1, -1),
    -1:Pos(-1, 0),
    1:Pos(1, 0),
    2:Pos(-1, 1),
    3:Pos(0, 1),
    4:Pos(1, 1),
}


# -- Functions -- 
def abrir_arq_chain(nome_arq):
    """
    Handles .chain file oppening

    return: 2D matrix
    """
    with open(nome_arq, "r") as arq:
        field = [linha.split() for linha in arq.readlines()]
        return field


def print_field(field: list[list], highlight: Pos = (-1, -1)):
    """
    Prints the 2D matrix representation

        of a given chain field
    """
    # TODO: Make it so it can have a highlight feature for a specific coordinate
    for i_row, row in enumerate(field):
        for i_col, el in enumerate(row):
            p_el = "\033[34m" + el + "\033[0m" if highlight == Pos(i_row, i_col) else el
            print(p_el, end = " ")
        print("")

def _edge_c(field):
    """
    Gives the postion of the of an edge Carbon
    """
    for i_row, row in enumerate(field):
        for i_col, el in enumerate(row):
            curr_pos = Pos(i_row, i_col)
            if el == 'C' and len(_scout(field, curr_pos)) == 1:
                # print(curr_pos)
                # mock_field = deepcopy(field)
                # mock_field[curr_pos.row][curr_pos.col] = '\033[34mX\033[0m'
                # rint_field(field, curr_pos)
                yield Pos(i_row, i_col)


def _scout(field, pos: Pos):
    nxt_el_l = []
    for dir in cd_offsets:
        # Next direction and connectin (to scout)
        nxt_step = cd_offsets[dir]
        nxt_con_pos = nxt_step + pos
        nxt_el_pos = nxt_con_pos + nxt_step
        if _is_inside(nxt_el_pos, len(field), len(field[0])):
            # If is inside, check to see if there is a Carbon in the next path
            nxt_el = field[nxt_el_pos.row][nxt_el_pos.col]
            next_con = field[nxt_con_pos.row][nxt_con_pos.col]
            if next_con.isnumeric() and next_con != '0'\
                and\
                nxt_el == 'C':
                nxt_el_l.append((dir, nxt_el_pos))       
    return nxt_el_l


def _is_inside(dir: Pos, n_row, n_col):
    if (dir.row >= n_row or dir.col >= n_col)\
        or\
        (dir.row < 0 or dir.col < 0):
        return False
    return True

def per_path(field: list[list], start_pos, recur = 0):
    """
    Go through the 2D representation of a carbon chain
    path
    recur: Current recursion
    """
    # Start position of this recrusion
    pos = start_pos
    origin_dir = 0
    # Temporary as it doesn't work with recursion!
    while True:
        print_field(field, pos)
        nxt_els = _scout(field, pos)
        nxt_n_con = len(nxt_els)
        if nxt_n_con == 2 or pos == start_pos:
            # If it has one possible path, move
            jump = False
            for nxt_dir, nxt_el_pos in nxt_els:
                if nxt_dir != origin_dir:
                    print(origin_dir)
                    origin_dir = -nxt_dir
                    pos = nxt_el_pos
                    jump = True
                    break
            if jump:
                continue
        if nxt_n_con > 2:
            # Recursion for multiple paths
            pass
        if nxt_n_con == 1:
            break



def per_chain(field):
    """
    Go through the entirety of a chain field representation
    """
    for pos in _edge_c(field):
        per_path(field, pos)

def main():
    path = "Chains/simple.chain"
    # Campo
    field = abrir_arq_chain(path)
    per_chain(field)
    list(_edge_c(field))


# -- Start --
if __name__ == "__main__":
    main()