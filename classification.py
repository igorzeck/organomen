# CLassification is done in 4 steps:
# 1. Atom by atom
# 2. Group by group
# 3. Host by Host (plus its groups)
# 4. Entire chain
# -- Imports --
from base_structures import Pos, HETEROATOMS, HALIDES, SIMPLE, DOUBLE, TRIPLE, QUADRUPLE  # ... No comments...
from pathfinder import scout, _get_host, run_subpath
from entities import Chain, Entity

# TODO: Here might be better to pass pos as list of POS instead of ids
# TODO: Multifunction procedure might be harder than it appears
# TODO: Need to be careful with counting! You start counting from one edge onwards!
# TODO: Treat cycles and ramifications like function

# -- Naming --
# - Name prefix -
# - Constants -

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

AFIXES = {
    'Ether':'oxi',
    'Radical':'hyl',
}

# Relative to chain functional class
# Maybe use a default dict?
SUFFIXES = {
    'Hydrocarbon':'e',
    'Alcohol':'ol',
    'Aldehyde':'al',
    'Acid':'oic',
    'Keton':'ona',
    'Ether':'e',
    'Halides':'e',  # Kinda hacky, not gonna lie...
    'UNRESOLVED':'NONE',
}

CON_Q = [
        '', # Possible for hydrocarbons
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

        # Define higher class of the compound
        # For example, cyclane and cyclene would have the same higher class (Hydrocarbon)
        # Fenol would be an Alcohol and etc.
        self.hclass: str = "UNCLASSIFIED"

        # Already verified entities pool
        self.ents_pool: set[int] = set()

        # First subgroup is main chain itself
        self.subgroups: list[tuple[Entity]] = [tuple(chain.main_chain)]

        self.classif: dict = {}  # TODO: initalize with some functionals keys (but empty sets)

        self.classif_by_host: dict[int] = {}  # Host: [functions]

        self.host_by_classif: dict[str] = {} # Functions: [host]
    
        # Function call
        self._root_question(0)

        self._define_subgroups()
    
    def _append_classif(self, _class: str, subg_id: int):
        if _class in self.classif:
            self.classif[_class].append(subg_id)
        else:
            self.classif[_class] = [subg_id]
        
        if subg_id != 0:
            _host_id = self.subgroups[subg_id][0].id

            if _class in self.host_by_classif:
                self.host_by_classif[_class].append(_host_id)
            else:
                self.host_by_classif[_class] = [_host_id]
            
            if _host_id in self.classif_by_host:
                self.classif_by_host[_host_id].append(_class)
            else:
                self.classif_by_host[_host_id] = [_class]

    def get_hclass(self) -> str:
        return self.hclass

    def get_classif(self, id: int) -> str:
        for classif in self.classif:
            if id in self.classif[classif]:
                return classif
        return 'UNRESOLVED'
    
    # - Chemical logic -
    def _define_subgroups(self):
        # TODO: Cleanup
        # Starts from a point on the main_chains
        start_id = 1 if self.get_classif(0) == "Cycle" else 0
        for main_ent in self.chain.main_chain[start_id:]:
            for nxt_con in main_ent.cons:
                nxt_ent = self.chain.chain[nxt_con.to_id]
                if (nxt_ent not in self.chain.main_chain):
                    _subgroup = run_subpath(self.chain, nxt_ent)
                    
                    for ent in _subgroup:
                        self.ents_pool.add(ent.id)
                    
                    self.subgroups.append(_subgroup)

                    print("Captured subgroup: ")
                    print(_subgroup)
                    
        # ID 0 is the main chain!
        for subg_id in range(1, len(self.subgroups)):
            self._root_question(subg_id)

        print("Classifications:", self.classif)

    # - Classification logic - 
    # (tree-like function cascade)
    # I think is better to have an if forest tbh...
    def _root_question(self, subg_id: int):
        # 0. Is it a cycle?
        if (self.subgroups[subg_id][0] == self.subgroups[subg_id][-1]) and\
            (len(self.subgroups[subg_id]) > 1):
            self._append_classif('Cycle', subg_id)
        # 1. Has an hteroatom?
        if any([ent.is_hetero() for ent in self.subgroups[subg_id]]):
            # print(f"{self.subgroups[subg_id]} has an heteroatom!")
            # If it does
            self._is_oxy(subg_id)
            self._is_halide(subg_id)
        else:
            # print(f"{self.subgroups[subg_id]} does NOT have an heteroatom!")
            # If it doesn't, is it the main chain?
            if subg_id == 0:
                # Tries to resolve main chain class
                self._resolve_hclass()
            else:
                # If not, treat as a radical
                self._append_classif('Radical', subg_id)
    
    # This way runs two times though...
    # - Oxygen -
    def _is_oxy(self, subg_id: int):
        if any([ent.el == 'O' for ent in self.subgroups[subg_id]]):
            self._is_oxy_2(subg_id)
    
    def _is_oxy_2(self, subg_id: int):
        for ent in self.subgroups[subg_id]:
            if ent == 'O':
                _host = _get_host(self.chain.main_chain, ent)
                if ent != _host:
                    if any([((con.type == SIMPLE))\
                            for con in ent.cons]) and\
                            self.chain.to_el(ent.id,'C') == 1:
                        # print('Is an alcohol!')
                        # self.classif[subg_id] = [_host.id, 'Alcohol']
                        self._append_classif('Alcohol', subg_id)
                    if any([(con.to_id == _host.id)\
                             and\
                            (con.type == DOUBLE)\
                            for con in ent.cons]):
                        # print('Is an aldehyde!')
                        self._append_classif('Aldehyde', subg_id)
                    if (self.chain.to_el(ent.id,'C') == 2):
                       # For now flags it if there is 2 connection to carbons
                       self._append_classif('Ether', subg_id)
    
    # - Halides -
    def _is_halide(self, subg_id: int):
        for ent in self.subgroups[subg_id]:
            if ent.el in HALIDES:
                _host = _get_host(self.chain.main_chain, ent)
                if ent != _host:
                    print(_host, end="\n\n")
                    if any([((con == SIMPLE))\
                            for con in ent.cons]):
                        # Then, it's an Halide
                        # Appends each element individually
                        self._append_classif(ent.el, subg_id)

    # - Main chain -
    # def _define_host_classif(self):
    #     self._append_classif('UNRESOLVED', 0)
    #     # Host dict
    #     # THIS "list" (dict really) comprehension woks! Wow
    #     host_dict: dict[str] = {_subg[0].id: [] for _subg in self.subgroups}

    #     for _subg_id, _subg in enumerate(self.subgroups):
    #         _host_id = _subg[0].id
    #         host_dict[_host_id].append(self.get_classif(_subg_id))
    

    def _resolve_hclass(self) -> str:
        for _host in self.classif_by_host:
            _host_classif = self.classif_by_host[_host]
            has_alcohol = 'Alcohol' in _host_classif
            has_ether = 'Ether' in _host_classif
            has_keton = 'Keton' in _host_classif
            if has_alcohol:
                self.hclass = 'Alcohol'
            if has_ether and has_keton:
                self.hclass = 'Ester'
            if has_ether:
                self.hclass = 'Ether'
            if has_keton:
                self.hclass = 'Keton'
        self.hclass = 'Hydrocarbon'

# - Function -
# - Helpers -
# Get connection infix
def con_infix(ids: int = [], hide_ids = True):
    n = len(ids)
    if n < len(CON_Q):
        if hide_ids or n == 0:
            return CON_Q[n]
        else:
            return '-' + ','.join(map(str, ids)) + '-'
    else:
        return 'UNSUPPORTED'

# - Main chain naming -
def _name_size_pref(n_main: int) -> str:
    if n_main <= len(PREFIXES):
        return PREFIXES[n_main]
    else:
        return 'NONE'


def _name_con_type(cons: list) -> str:
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


def _name_radical(classific: Classifier) -> str:
    final_str = ''

    # TODO: each class in a function and auxiliary functions should be made!
    
    # - Non hydrocarbon radicals -
    # Sorts to make it easy on naming compounding
    for gp in sorted(classific.classif):
        # For now jumps "Radical"
        if gp == 'Radical':
            continue

        # If it is a Haladie
        if gp in HALIDES:
            # Name of the halide
            for i_radical in classific.classif[gp]:
                # For now looks into the final element
                el = classific.subgroups[i_radical][-1].el
                n_str = HALIDES[el].lower()
                if final_str:
                    final_str += '-'
            
                # - Position (of the host) -
                func_pos = []

                _host = classific.subgroups[i_radical][0]
                func_pos.append(str(classific.chain.get_main_path_id(_host.id)))
            
            # - Functional name -            
            final_str += ','.join(func_pos) + '-' + n_str

    # - Hydrocarbon radicals -
    # Apparently always specify the position of the radical!
    if 'Radical' in classific.classif:
        # TODO: function for this logic!
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
            
            # TODO: function for this logic!
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
            final_str += AFIXES[gp]

    if 'Ether' in classific.classif:
        gp = 'Ether'

        radical_types: dict = {}

        # Fills radical_types (based on size)
        for i_radical in classific.classif[gp]:
            # Note that i_subg is also the id for the host of the group!
            # -2 to ignore 'host' and oxygen
            n_els = len(classific.subgroups[i_radical]) - 2
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
            
            n_str = ''
            if rt_len < len(CON_Q):
                n_str += CON_Q[rt_len]
            # - Functional name -
            if 'Cycle' in classific.classif:
                if i_radical in classific.classif['Cycle']:
                    final_str += 'cycle'
            final_str += r_prefix
            final_str += AFIXES[gp]

    return final_str


def _name_suffix(classific: Classifier) -> str:
    final_str = ''
    main_classif = classific.get_hclass()
    # - Functional positions (ids) -
    _main_ids = []
    if main_classif in classific.host_by_classif:
        # Counts position
        _ids = classific.host_by_classif[main_classif]
        _main_ids = list(classific.chain.get_main_path_ids(_ids))
    final_str += con_infix(_main_ids, hide_ids=False)
    # - Functional suffix -
    final_str += SUFFIXES[classific.get_hclass()]
    return final_str


def class_chain(chain: Chain):
    if not chain.main_path:
        return "Chain is empty!"

    prefix = infix = suffix = ''
    # Clasificates chain
    # TODO: Make it less redundant... Two ways to classify huh???
    # chain.functional = _get_class(chain)
    _classfier = Classifier(chain)

    prefix = _name_radical(_classfier)

    # If is a cycle (main chain)
    if 'Cycle' in _classfier.classif:
        if 0 in _classfier.classif['Cycle']:
            prefix += 'cycle'

    # Adds a '-' if first word is a consonant
    # For now only look into unique ids (sets) as its possible to have cycles!
    
    # For now jumps over if an ether
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

    suffix += _name_suffix(_classfier)
    
    chain.name = prefix + infix + suffix
    return prefix + infix + suffix
