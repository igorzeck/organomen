from constants import *
import warnings

# -- Classes --
class Pos:
    # TODO: Upgrade to 3D
    """
    2D Vector coordinate
    """
    def __init__(self, row: int = 0, col: int = 0):
        self.row = row
        self.col = col
    # Math
    def __add__(self, value):
        return Pos(self.row + value.row, self.col + value.col)
    def __sub__(self, value):
        return Pos(self.row - value.row, self.col - value.col)
    def __mul__(self, value: int):
        # Scalar multiplications
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


class Pos3D:
        # TODO: Upgrade to 3D
    """
    3D Vector coordinate - Initializes to origin
    """
    # NOTE: For now is 2D too!
    # TODO: Make it 3D!
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y
    # Math
    def __add__(self, value):
        return Pos3D(self.x + value.x, self.y + value.y)
    def __sub__(self, value):
        return Pos3D(self.x - value.x, self.y - value.y)
    def __mul__(self, value: int):
        # Scalar multiplications
        return Pos3D(self.x * value, self.y * value)
    def __truediv__(self, value):
        if isinstance(value, (int, float)):
            return Pos3D(self.x / value, self.y / value)
        elif isinstance(value, Pos3D):
            return Pos3D(self.x / value.x, self.y / value.y)
    def __pow__(self, value):
        return Pos3D(self.x ** value, self.y ** value)
    # Representations
    def __str__(self):
        return f"({self.x}; {self.y})"
    def __repr__(self):
        return self.__str__()
    def __iter__(self):
        for p in (self.x, self.y):
            yield p
    # Relationals
    def __eq__(self, value):
        # Float comparison though, should be implemented!
        return (self.x == value.x and self.y == value.y)
    

# -- Functions --
# - File handling -
def open_f(filepath: str):
    """
    Handle oppening of valid filepaths
    
    :param filepath: string with filepath
    :type filepath: str
    :return: returns field representation of the file as well as type of file
    """
    # - Check if file exists -
    file: Path = Path(filepath)

    if not file.exists():
        raise FileExistsError(f"File '{file.name}' don't exist!")
    
    # - Checks extension -
    _ext = file.suffix

    if _ext not in EXTS:
        raise TypeError(f"File extension '{_ext}' is not valid! Supported extensions:\n{EXTS}")

    # - Recent files -
    recent_files_path: Path = CONF_PATH / 'recent_files'
    
    if recent_files_path.is_file():
        with open(recent_files_path, "r+") as f:
            # See if content is the same (For now keeps one recent file)
            if filepath not in f.readlines():
                f.seek(0)
                f.write(filepath)
                f.truncate() # If size is not specified it is the IO position!
  
    else:
        warnings.warn(f"'{recent_files_path}' not found!")
    
    # Maybe a default dict with functions?
    match _ext:
        case '.field':
            return _open_f_field(filepath), 'field'
        case '.smi':
            # Reading as a single string here, perhaps?
            return _open_smile(filepath), 'smile'
        case _:
            # Impossible to get here, to be honest
            raise TypeError(f'\'{_ext}\' is not a supported extension!')


def _open_f_field(filepath: Path) -> tuple[tuple[str]]:
    with open(filepath, "r") as arq:
        field = tuple([tuple([el.upper() for el in linha.split()]) for linha in arq.readlines() if linha[0] != COMMENT_WILDCARD])
        # Upper case everything
        return field


def _open_smile(filepath: Path):
    # - Opening of the file - 
    smile_str = ''
    with open(filepath, 'r') as f:
        # Assumes 1 single line
        smile_str = f.readline()
    
    if not smile_str:
        warnings.warn('Empty file')
    return smile_str
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
                el = el if not show_ids else str(highlights.index(curr_pos) + 1)
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
# TODO: For 3D structure I'll have to do away with all this nonsense
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