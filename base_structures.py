# -- Classes --
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

# -- Functions --
# - File handling -
def abrir_arq(nome_arq: str) -> tuple[tuple[str]]:
    # Checks extension
    _ext = nome_arq.partition('.')[-1].lower()
    if _ext not in exts:
        raise TypeError(f'File extension is not valid! Supported extensions:\n{exts}')
    # Maybe a default dict with functions?
    match _ext:
        case 'field':
            return _abrir_arq_field(nome_arq)
        case _:
            raise TypeError(f'\'{_ext}\' is not a supported extension!')

def _abrir_arq_field(nome_arq: str) -> tuple[tuple[str]]:
    """
    Handles .field file opening
    
    :param nome_arq: The name of the .field to open
    :type nome_arq: str
    :return: A 2D matrix
    :rtype: tuple[tuple[str]]
    """
    with open(nome_arq, "r") as arq:
        field = tuple([tuple([el.upper() for el in linha.split()]) for linha in arq.readlines() if linha[0] != COMMENT_WILDCARD])
        # Upper case everything
        return field


# - Visual -
# TODO: Still using field I see...
# TODO: Change it to show on chain path instead
def print_field(field: list[list],
                highlights: tuple = [Pos(-1, -1)],
                highlight_color_only = False,
                show_ids = False,
                zeros_repr='.'):
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
                el = el if not show_ids else str(highlights.index(curr_pos))
                if not highlight_color_only:
                    if len(highlights) > 1:
                        if curr_pos == highlights[0] == highlights[-1]:
                            el = 'o' + el
                        elif curr_pos == highlights[0]:
                            # Start
                            el = 'i' + el
                        elif curr_pos == highlights[-1]:
                            # End
                            el = 'f' + el
                p_el = "\033[34m" + el.ljust(3) + "\033[0m"
            else: 
                p_el = el.ljust(3)
            print(p_el, end = "")
        print("")
    print("")

# -- Constants --
# - File related -
COMMENT_WILDCARD = '#'
exts = ['field']
# Enun would be just lovely, but alas...
# Connections directions
NORTHWEST = -4
NORTH = -3
NORTHEAST = -2
WEST = -1
EAST = 1
SOUTHWEST = 2
SOUTH = 3
SOUTHEAST = 4

# Connection type
SIMPLE = 1
DOUBLE = 2
TRIPLE = 3
QUADRUPLE = 4


# Atoms in all caps for ease of comparison
HALIDES = {
    'F':'Fluorine',
    'CL':'Clhorine',
    'BR':'Bromine',
    'I':'Iodine'
}

# Merges dictionaries ( Python >= 3.9.0 )
HETEROATOMS = HALIDES | {'O':'Oxygen'}

# Coordinates associated with each connection
# -4 -3 -2
# -1  C  1
#  2  3  4
CD_OFFS = {
    -4:Pos(-1, -1),
    -3:Pos(-1, 0),
    -2:Pos(-1, 1),
    -1:Pos(0, -1),
    1:Pos(0, 1),
    2:Pos(1, -1),
    3:Pos(1, 0),
    4:Pos(1, 1),
}

# Failsafe for iterative routines
MAX_ITER = 21