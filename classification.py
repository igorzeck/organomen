# CLassification is done in 4 steps:
# 1. Atom by atom
# 2. Group by group
# 3. Host by Host (plus its groups)
# 4. Entire chain
# -- Imports --
from base_structures import Pos
from pathfinder import scout, _get_host, run_subpath
from entities import Chain, Entity

# TODO: Here might be better to pass pos as list of POS instead of ids
# TODO: Multifunction procedure might be harder than it appears
# TODO: Need to be careful with counting! You start counting from one edge onwards!
# TODO: Treat cycles and ramifications like function

# -- Naming --
# - Name prefix -
# - Constants -
SIMPLE = '1'
DOUBLE = '2'

HETEROATOMS = [
    'F',
    'O',
    'N'
]

# Functional look-up table with composite key index
functional_sub = {
    'Hydrocarbon':'',
    'Alcohol':'Alcohol'
}

# Relative to N of Cs in main chain
PREFIXES = [
    '',
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
    'dodec',
    'tridec',
    'tetradec',
    'pentadec',
    'hexadec',
    'heptadec',
    'octadec',
    'nonadec',
    'icos'
]

INFIXES = [
    'UNDER',
    'an',
    'en',
    'yn'
]

# Relative to chain functional class
SUFIXES = {
    'Hydrocarbon':'e',
    'Alcohol':'ol',
    'Aldehyde':'al',
    'Acid':'oico',
    'Keton':'ona',
    'Ether':'ico',
    # Radicals
    'Radical':'hyl',
    'Unsupported!':'NONE',
}

CON_Q = [
        'UNDER',
        '',
        'di',
        'tri'
    ]
# For now infixes are not included
# - Class -
# Necessary to find a way to integrate this
class Classifier:
    """
    Classifier object for an atom
    """
    # - Class related -
    def __init__(self, chain: Chain):
        self.chain: Chain = chain

        # Already verified entities pool
        self.ents_pool: set[int] = set()

        # First subgroup is main chain itself
        self.subgroups: list[tuple[Entity]] = [tuple(chain.main_chain)]

        self.classif: dict = {}  # TODO: initalize with some functionals keys (but empty sets)

        # Function call
        self._define_subgroups()
    
    def _append_classif(self, _class: str, subg_id: int):
        if _class in self.classif:
            self.classif[_class].append(subg_id)
        else:
            self.classif[_class] = [subg_id]

    # - Chemical logic -
    def _define_subgroups(self):
        for ent in self.chain:
            if (ent.id not in self.ents_pool):
                if (ent not in self.chain.main_chain):
                    # Maybe I should treat the main chain as a subgroup!
                    _subgroup = run_subpath(self.chain, ent)
                    
                    for ent in _subgroup:
                        self.ents_pool.add(ent.id)
                    
                    self.subgroups.append(_subgroup)
                    print(_subgroup)
                    
        # ID 0 is the main chain!
        for subg_id in range(len(self.subgroups)):
            self.root_question(subg_id)

        print("Classifications:", self.classif)

    # - Classification logic - 
    # (tree-like function cascade)
    # I think is better to have an if forest tbh...
    def root_question(self, subg_id: int):
        # 0. Is a cycle?
        if (self.subgroups[subg_id][0] == self.subgroups[subg_id][-1]):
            self._append_classif('Cycle', subg_id)
        # 1. Is it an hteroatom?
        if any([ent.el in HETEROATOMS for ent in self.subgroups[subg_id]]):
            # print(f"{self.subgroups[subg_id]} has an heteroatom!")
            # If it is
            self.is_oxy(subg_id)
        else:
            # print(f"{self.subgroups[subg_id]} does NOT have an heteroatom!")
            # If it isn't, is it the main chain?
            if subg_id == 0:
                # If so, then is a standard hydrocarbon chain
                self._append_classif('Hydrocarbon', subg_id)
            else:
                # If not, classificate it
                self.which_radical(subg_id)
                self._append_classif('Radical', subg_id)
    
    def which_radical(self, subg_id: int):
        pass

    def is_oxy(self, subg_id: int):
        if any([ent.el == 'O' for ent in self.subgroups[subg_id]]):
            self.is_oxy_2(subg_id)
    
    def is_oxy_2(self, subg_id: int):
        for ent in self.subgroups[subg_id]:
            if ent == 'O':
                _host = _get_host(self.chain.main_chain, ent)
                if ent != _host:
                    print(_host, end="\n\n")
                    if any([((con.type == SIMPLE))\
                            for con in ent.cons]):
                        # print('Is an alcohol!')
                        # self.classif[subg_id] = [_host.id, 'Alcohol']
                        self._append_classif('Alcohol', subg_id)
                    if any([(con.to_id == _host.id)\
                             and\
                            (con.type == DOUBLE)\
                            for con in ent.cons]):
                        # print('Is an aldehyde!')
                        self._append_classif('Aldehyde', subg_id)


# - Function -
# - Main chain naming -
def _name_size_pref(n_main: int):
    if n_main <= len(PREFIXES):
        return PREFIXES[n_main]
    else:
        return None


def _name_con_type(cons: list):
    # For now always show position, even when redundant!
    con_qte = [
        '',
        'adi',
        'atri'
    ]

    infix = ''

    # Indexes within chain
    i2 = [str(i + 1) for i, con in enumerate(cons) if con == '2']
    i3 = [str(i + 1) for i, con in enumerate(cons) if con == '3']
    n2 = len(i2)
    n3 = len(i3)

    # This way 'enin' can occur naturally
    # TODO: Correct it so it's only between carbon to carbon connections
    if n2 != 0:
        if len(con_qte) >= n2:
            infix += con_qte[n2 - 1]
        infix += '-' + ','.join(i2) + '-' + INFIXES[2]
    if n3 != 0:
        if len(con_qte) >= n3:
            infix += con_qte[n3 - 1]
        infix += '-' + ','.join(i3) + '-' + INFIXES[3]
    if n2 == n3 == 0:
        infix += INFIXES[1]
    return infix


def _name_functional(functional: dict[list[int]]):
    final_str = ''
    for funct in functional:
        # Checks size
        funct_info = functional[funct]
        n_gps = len(funct_info)
        if n_gps >= 2:
            func_pos = [str(_pos) for _pos in funct_info]
            n_str = ''
            if n_gps < len(CON_Q):
                n_str += CON_Q[n_gps]
            final_str += '-' + ','.join(func_pos) + '-' + n_str
        final_str += SUFIXES[funct]
    return final_str

def _name_radical(classific: Classifier):
    final_str = ''
    # Apparently always specify the position of the radical!
    if 'Radical' in classific.classif:
        gp = 'Radical'

        # For now only count atom in the radical subgroup!
        # Separate elements into different types
        radical_types: dict = {}

        # Fills radical_types (based on size)
        for i_radical in classific.classif[gp]:
            # Note that i_subg is also the id for the host of the group!
            # -1 to ignore 'host'
            n_els = len(classific.subgroups[i_radical]) - 1
            if n_els in radical_types:
                radical_types[n_els].append(i_radical)
            else:
                radical_types[n_els] = [i_radical]
       
        # Sorts radical_types ids in their names alphabetical order
        func_names = list(map(lambda k: PREFIXES[k] if k < len(PREFIXES) else f"UNDEF_LEN[{k}]", radical_types.keys()))
        # Zips names with their values
        pair_names = list(zip(func_names, radical_types.keys()))
        # Orders based on first element (a.k.a its name)
        pair_names.sort(key=lambda el: el[0])
        # Gets last element (id)
        # Note that the multiplier prefixes do not count on the alphabetical order!

        print("Pares:", pair_names)

        for r_prefix, i_radical in pair_names:
            # Use radical types to name the prefixes
            # Position
            rt_len = len(radical_types[i_radical])

            print(f"Type: {i_radical} | n: {rt_len}")
            
            # The zero index is the 'host'
            hosts_ids = [_id for _id in [classific.subgroups[_i_subg][0].id for _i_subg in radical_types[i_radical]]]
            func_pos = [str(classific.chain.get_main_path_id(_pos)) for _pos in hosts_ids]
            
            n_str = ''
            if rt_len < len(CON_Q):
                n_str += CON_Q[rt_len]
            if final_str:
                final_str += '-'
            final_str += ','.join(func_pos) + '-' + n_str

            # - Functional name -
            if 'Cycle' in classific.classif:
                if i_radical in classific.classif['Cycle']:
                    final_str += 'cycle'
            final_str += r_prefix
            final_str += SUFIXES[gp]
    
    return final_str
            

# TODO: Make id go over the entire chain, not only main one
def _classif_atom(field: list[Pos], ids, pos_id):
    pos = ids[pos_id]
    el = field[pos.row][pos.col]
    info = scout(field, ids, pos_id)
    els = [field[ids[id[0]].row][ids[id[0]].col] for id in info]

    if el != 'O' and el !='C':
        return 'Unsupported!', pos_id
    
    # - Radicals -
    if (el == 'C'):
        if els.count('C') > 2:
            return 'Radical', pos_id

    # - Oxygen -
    # 1. How many Cs
    if els.count('C') > 1:
        return 'Ether', pos_id

    # Gets host info
    host_info = []
    host_id = -1

    # TODO: Fix host id being different from real host if it is in a subgroup
    # Treating subgroup as a function would be nice (including ignoring its elements for classification)
    for id in info[0]:
        nxt_info = scout(field, ids, id)
        for nxt_id in nxt_info:
            nxt_el = field[ids[nxt_id[0]].row][ids[nxt_id[0]].col]
            if nxt_el == 'C':
                host_info = scout(field, ids, id, nxt_con_val=True)
                host_id = id
                break
        if host_id != -1:
            break
    
    host_els = [field[ids[id].row][ids[id].col] for id, _ in host_info]

    # 2. Is it a edge C?
    c_edge = host_els.count('C') == 1

    # 3. Does the host connect to other Os?
    c_alone = host_els.count('O') == 1

    # 4. How many double connections with Os?
    c_n_double = sum([con == '2' and field[ids[id].row][ids[id].col] == 'O' for id, con in host_info])
    
    # To be honest I'm not happy with this, but it works
    func_name = 'Unknown'
    if c_n_double == 1:
        if c_edge:
            if c_alone:
                func_name = 'Aldehyde'
            else:
                func_name = 'Acid'
        elif c_alone:
            func_name = 'Keton'
    elif c_n_double == 0:
        func_name = "Alcohol"

    return func_name, host_id


def _get_class(chain: Chain):
    functionals: dict = {}

    field = chain.field
    
    atoms_dict: dict = {}
    for pos_id, pos in enumerate(chain.id_pool):
        el = field[pos.row][pos.col]
        if el in atoms_dict:
            atoms_dict[el].append(pos_id)
        else:
            atoms_dict[el] = [pos_id]
    
    print(atoms_dict)

    # Decides in its functional class
    if any(el == 'C' for el in atoms_dict):
        if all(el == 'C' for el in atoms_dict):
            # functionals['Hydrocarbon'] = {chain.main_path[0]}
            return {'Hydrocarbon': {chain.main_path[0]}}
        else:
            atoms_dict.pop('C')
    else:
        return "not organic"

    # Heteroatoms
    for pos_id, pos in enumerate(chain.id_pool):
        el = field[pos.row][pos.col]
        if el in atoms_dict:
           # Note that it uses the host ID
           _class, _host_id = _classif_atom(field, chain.id_pool, pos_id)
           _host_id_path = chain.get_main_path_id(_host_id)
           if _class in functionals:
               functionals[_class].add(_host_id_path)
           else:
               functionals[_class] = {_host_id_path}

    return functionals


def class_chain(chain: Chain):
    if not chain.main_path:
        return "Chain is empty!"
    # Clasificates chain
    chain.functional = _get_class(chain)
    _classfier = Classifier(chain)

    prefix = _name_radical(_classfier)

    # If is a cycle (main chain)
    if 'Cycle' in _classfier.classif:
        if 0 in _classfier.classif['Cycle']:
            prefix += 'cycle'

    # Adds a '-' if first word is a consonant
    # For now only look into unique ids (sets) as its possible to have cycles!
    main_prefix = _name_size_pref(len(set(chain.main_path)))

    # Needs to be better understood when to use the hyphen!
    # if not main_prefix[0] in ['a', 'e', 'i', 'o', 'u']:
    #     prefix += '-'
    
    prefix += main_prefix

    if not prefix:
        return "Chain is too big!"

    # Makes chain
    cons = []
    old_pos = None
    for pos in [chain.id_pool[id] for id in chain.main_path]:
        if old_pos:
            # This would be perfect in a function!
            dif_pos = pos - old_pos
            norm_dif = dif_pos/2
            con_pos = old_pos + norm_dif
            cons.append(chain.field[con_pos.row][con_pos.col])
        old_pos = pos
    
    # TODO: Use connections already stored in chain object?
    infix = _name_con_type(cons)

    # sufix = _name_func_class(_classfier)
    sufix = _name_functional(chain.functional)

    # TODO: Hopefully temporary
    chain.name = prefix + infix + sufix
    return prefix + infix + sufix
