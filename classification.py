# CLassification is done in 4 steps:
# 1. Atom by atom
# 2. Group by group
# 3. Host by Host (plus its groups)
# 4. Entire chain
# -- Imports --
from base_structures import Pos, print_field
from constants import *
from pathfinder import scout, _get_host, iterate_subpaths
from entities import Chain, Entity
from auxiliary import pair_subgs
from collections import defaultdict

# TODO: Here might be better to pass pos as list of POS instead of ids
# TODO: Multifunction procedure might be harder than it appears
# TODO: With one atom it shouldn't show position of the radicals 

# -- Classification --
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
        self.hclass: str = {}

        # First subgroup is main chain itself
        self.subgroups: list[tuple[Entity]] = [tuple(chain.main_chain)]

        self.classif: dict = {}  # TODO: initalize with some functionals keys (but empty sets)

        self.classif_by_host: dict[int] = {}  # Host: [functions]

        self.host_by_classif: dict[str] = {} # Functions: [host]

        # Classificates main chain
        self._root_question(0)
    
        # Function call
        self._define_subgroups()
    
    def _append_classif(self, _class: str, subg_id: int):
        if _class in self.classif:
            self.classif[_class].append(subg_id)
        else:
            self.classif[_class] = [subg_id]
        
        # if subg_id != 0:
        _host_id = self.subgroups[subg_id][0].id

        if _class in self.host_by_classif:
            self.host_by_classif[_class].append(_host_id)
        else:
            self.host_by_classif[_class] = [_host_id]
        
        if _host_id in self.classif_by_host:
            self.classif_by_host[_host_id].append(_class)
        else:
            self.classif_by_host[_host_id] = [_class]

    # TODO: Make it less awkward..
    def get_hclass(self) -> str:
        if self.hclass:
            return list(self.hclass.keys())[0]
        else:
            return ""
    
    def get_trad_hclass(self, plural: bool = False) -> str:
        """
        Returns translated higher classification (for print purpose)
        
        :param self: Classifier object
        :return: Translated higher class
        :rtype: str
        """
        return list(self.hclass.values())[0] + (PLURAL_STR if plural else '')

    def get_classif(self, id: int) -> str:
        for classif in self.classif:
            if id in self.classif[classif]:
                return classif
        return 'UNRESOLVED'  # Maybe adding an immutable constant?

    def get_trad_classif(self, id: int) -> str:
        for classif in self.classif:
            if id in self.classif[classif]:
                return CLASSIFICATION[classif]
        return UNRESOLVED
    
    # - Chemical logic -
    def _define_subgroups(self):
        # TODO: Cleanup
        # Starts from a point on the main_chains
        cyclical = self.get_classif(0) == 'Cycle'
        self.subgroups.extend(iterate_subpaths(self.chain, cyclical))
                    
        # ID 0 is the main chain!
        for subg_id in range(1, len(self.subgroups)):
            self._root_question(subg_id)

        # Resolves main class
        self.hclass = self._resolve_hclass()
        if self.get_classif(0) == 'UNRESOLVED':
            self._append_classif('Chain', subg_id)
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
            self._is_nitrogenic(subg_id)
            self._is_oxy(subg_id)
            self._is_halide(subg_id)
        else:
            # print(f"{self.subgroups[subg_id]} does NOT have an heteroatom!")
            # If it doesn't, is it the main chain?
            if subg_id != 0:
                # If not, treat as a radical
                self._append_classif('Radical', subg_id)
    
    # This way runs two times though...
    # TODO: Merge it all in a single function!
    # - Nitrogenic -
    def _is_nitrogenic(self, subg_id: int):
        if any([ent.el == 'N' for ent in self.subgroups[subg_id]]):
            for ent in self.subgroups[subg_id]:
                if ent == 'N':
                    if len(self.subgroups[subg_id]) > 1:  # If not only one single nitrogen
                        if (self.chain.to_el(ent.id,'O') == 2) and\
                                (self.chain.to_el(ent.id,'C') == 1):
                            self._append_classif('Nitro', subg_id)
                        else:
                            self._append_classif('Radical', subg_id)
                    elif subg_id == 0: # Only the main chain can get here
                        # TODO: Add "Amide"
                        self._append_classif('Amine', subg_id)

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
    def _resolve_hclass(self) -> str:
        for _host in self.classif_by_host:
            _host_classif = self.classif_by_host[_host]
            is_amine = 'Amine' in _host_classif
            is_amide = 'Amide' in _host_classif
            has_alcohol = 'Alcohol' in _host_classif
            has_ether = 'Ether' in _host_classif
            has_keton = 'Keton' in _host_classif

            if is_amine or is_amide:
                return {'Nitrogenic':CLASSIFICATION['Nitrogenic']}
            if has_alcohol:
                return {'Alcohol':CLASSIFICATION['Alcohol']}
            if has_ether and has_keton:
                return {'Ester':CLASSIFICATION['Ester']}
            if has_ether:
                return {'Ether':CLASSIFICATION['Ether']}
            if has_keton:
                return {'Keton':CLASSIFICATION['Keton']}
        return {'Hydrocarbon':CLASSIFICATION['Hydrocarbon']}

# -- Naming --
# - Function -
# - Helpers -
# Get connections infix
def mult_prefix(ids: int = [], hide_ids = True, trailing_hyphen = True):
    n = len(ids)
    _trailing = '-' if trailing_hyphen else ''
    _mult_prefix_str = ''
    if n < len(MULT_PREFS):
        if (not hide_ids) and (n > 0):
            _mult_prefix_str = _trailing + ','.join(map(str, ids)) + '-'
        if n < len(MULT_PREFS):
            _mult_prefix_str += MULT_PREFS[n]
        else:
            _mult_prefix_str += "NONE"
    else:
        _mult_prefix_str = 'UNSUPPORTED'
    return _mult_prefix_str

get_prefix = lambda n: PREFIXES[n] if n < len(PREFIXES) else UNDEF

# - Main chain naming -
def _name_size_pref(n_main: int) -> str:
    if n_main <= len(PREFIXES):
        return PREFIXES[n_main]
    else:
        return UNDEF


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

def __pair_func_pos(host_ids: int, subgs: list[int], opt_suffix:str =''):
    """
    Pair function name and its position on the main path id.
    
    :param host_ids: Ids in the main path
    :type host_ids: int
    :param subgs: List of subgroups global ids
    :type subgs: list[int]
    :param suffix: Optional parameter to add to the function name
    :type suffix: str
    """
    # Pair the main path host ids with their subgroups
    zip_subg = zip(host_ids, subgs)
    radical_types = defaultdict(list)  # Defaults to an empty list

    # Fills radical_types (based on size)
    for _host_id, _subg in zip_subg:
        # Note that i_subg is also the id for the host of the group!
        # len - 1 to ignore 'host'
        n_els = len(_subg) - 1
        n_type = get_prefix(n_els) + opt_suffix
        radical_types[n_type].append(_host_id)
    
    # Orders radical_type (type, ids) tuple based on alphabetical order
    # NOTE: -thyl is the same for every radical, so don't matter on the ordering
    return dict(sorted(radical_types.items(), key=lambda k: k[0]))


def __name_pairs(gp_pairs: dict, hide_ids = False, trailing_first_hyphen = False) -> str:
    final_str = ''
    trailing_hyphen = trailing_first_hyphen
    # Name of the halide
    # Note they are already sorted in alphabetical order
    # and multiplier prefixes don't count towards aplhabetical ordering
    # Adds to string with multiplier
    for _gp_name in gp_pairs:
        _gp_ids = gp_pairs[_gp_name]
        final_str += mult_prefix(_gp_ids, hide_ids, trailing_hyphen=trailing_hyphen) + _gp_name.lower()
        trailing_hyphen = True
    return final_str
    

def _name_substitutive(classifier: Classifier, hide_ids = False) -> str:
    final_str = ''
    radical_str = ''
    halides_str = ''
    ether_str = ''

    # - Non hydrocarbon radicals -
    # Should they be dictionaries or just tuples?
    gp_halides = {}
    gp_radicals = {}
    for gp in sorted(classifier.classif):
        _gp_name = ''
        _gp_host = []
        
        for _host in classifier.host_by_classif[gp]:
            _gp_host.append(classifier.chain.get_main_path_id(_host))

        if gp in HALIDES:
            _gp_name = HALIDES[gp]

            gp_halides[_gp_name] = _gp_host
        if gp == 'Radical':
            # Separate into groups
            subgs = [classifier.subgroups[subg_id] for subg_id in classifier.classif[gp]]

            gp_radicals = __pair_func_pos(_gp_host, subgs, AFIXES[gp])
        if gp == 'Ether':
            # Name it directly for now
            print(classifier.subgroups)
            print(classifier.classif[gp])

            # -2 to avoid counting the oxygen AND the host
            ether_str += get_prefix(len(classifier.subgroups[-1]) - 2) + AFIXES[gp]


    halides_str = __name_pairs(gp_halides, hide_ids=hide_ids)
    final_str += halides_str

    radical_str = __name_pairs(gp_radicals, hide_ids=hide_ids, trailing_first_hyphen=final_str != '')
    final_str += radical_str

    final_str += ether_str

    return final_str


def _name_suffix(classific: Classifier, hide_ids = False) -> str:
    final_str = ''
    main_classif = classific.get_hclass()
    # - Functional positions (ids) -
    _main_ids = []
    if main_classif in classific.host_by_classif:
        # Counts position
        _ids = classific.host_by_classif[main_classif]
        _main_ids = list(classific.chain.get_main_path_ids(_ids))
    final_str += mult_prefix(_main_ids, hide_ids)
    # - Functional suffix -
    final_str += SUFFIXES[classific.get_hclass()]
    return final_str


def class_chain(chain: Chain):
    if not chain.main_path:
        return "Chain is empty!"

    prefix = qte_prefix = infix = suffix = ''
    # Clasificates chain
    # TODO: Make it less redundant... Two ways to classify huh???
    # chain.functional = _get_class(chain)
    _classfier = Classifier(chain)
    
    _is_single_atom = len(chain) == 1
    _is_nitrogenic = _classfier.get_hclass() == 'Nitrogenic'
    _is_oxygenic = _classfier.get_hclass() == 'Ester' or _classfier.get_hclass() == 'Ether'
    
    _hide_id = _is_single_atom or _is_nitrogenic

    prefix = _name_substitutive(_classfier, hide_ids=_hide_id)

    # If is a cycle (main chain)
    if 'Cycle' in _classfier.classif:
        if 0 in _classfier.classif['Cycle']:
            prefix += AFIXES['Cycle']

    # Adds a '-' if first word is a consonant
    # For now only look into unique ids (sets) as its possible to have cycles!
    
    if not _is_nitrogenic:
        qte_prefix = _name_size_pref(len(set(chain.main_path)))
    else:
        qte_prefix = AFIXES[_classfier.get_classif(0)]

    # Needs to be better understood when to use the hyphen!
    # if not main_prefix[0] in ['a', 'e', 'i', 'o', 'u']:
    #     prefix += '-'
    
    prefix += qte_prefix

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
    if not _is_nitrogenic:
        infix = _name_con_type(cons)

    suffix += _name_suffix(_classfier, _is_oxygenic)
    
    chain.name = prefix + infix + suffix
    chain.func_name = "NOT IMPLEMENTED!"
    # TODO: Cyclan, Cyclin and whatnot?
    # I don't like it. I think is best to separete it!
    if REVERSE_STR:
        chain.functional = _classfier.get_trad_classif(0) + CON_STR + _classfier.get_trad_hclass(plural=True)
    else:    
        chain.functional = _classfier.get_trad_hclass(plural=True) + CON_STR + _classfier.get_trad_classif(0)
