from base_structures import Pos, print_field, abrir_arq_chain

# - Connection -
# TODO: Não é melhor por referência logo?
class Connection:
    # - Class related -
    def __init__(self, from_id: int, to_id: int, dir: int, type: int):
        self.from_id = from_id
        self.to_id = to_id
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
        self.id = ent_id
        self.el = el
        self.cons = cons
    # - Visual -
    def __str__(self):
        return f"({self.id}: {self.el} - {self.cons})"
    def __repr__(self):
        return f"({self.id}: {self.el})"

# - Classes -
class Chain:
    def __init__(self, builder: str | list):
        """
        Paths are represented by an List of integer IDs
        builder: field array or path string
        """
        # 
        if isinstance(builder, str):
            self.load_file(builder)
        elif isinstance(builder, list):
            self.initiate(builder)
        else:
            raise TypeError('Builder is neither an file path (str) or field\'s list!')
    
    # - Class related -
    def initiate(self, field: list[str], entities: tuple[Entity] = None):
        # - Field -
        self.field = field
        self.n_row = len(self.field)
        self.n_col = len(self.field[0])

        self.id_dict: dict[Pos] = self._make_id_dict()
        # List all of paths
        self.paths = []
        # List of indexes for the main paths (name-definer)
        self.main_path = []

        # ~ Test ~
        self.chain = entities

        # - Chemical properties -
        self.name = ""
        
        self.functional = ""
        self.func_name = ""


    def load_file(self, arq: str):
        field: list[list[str]] = abrir_arq_chain(arq)
        chain = []
        # Create chain representation of the field
        ent_id = 0
        for row in field:
            for el in row:
                if not el.isnumeric():
                    
                    chain.append(Entity(
                        ent_id,
                        el,

                    ))
                ent_id += 1
        self.initiate(field, tuple(chain))
    
    # - Construct -
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
        summary = f"Name: {self.name}\
                   Functional class: {self.functional}"
        return summary
    def __repr__(self):
        return self.__str__()
    
    # - Functionalities -
    def __len__(self):
        return len(self.paths)