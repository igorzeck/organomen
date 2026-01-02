# CLassification is done in 4 steps:
# 1. Atom by atom
# 2. Group by group
# 3. Host by Host (plus its groups)
# 4. Entire chain
# -- Imports --
from base_structures import Pos, print_field
from constants import *
from pathfinder import scout, _get_host, iterate_subpaths
from entities import Chain, Entity, Connection
from collections import defaultdict

# TODO: Here might be better to pass pos as list of POS instead of ids
# TODO: Multifunction procedure might be harder than it appears
# TODO: With one atom it shouldn't show position of the radicals 

# -- Classification --
# - Class -
# Necessary to find a way to integrate this
# Should detect info on connection too!
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
        self._hclass = {}
        self.subfunction = {}
        # NOTE: Cyclic only for the main chain
        self.main_cyclical: bool = False

        # First subgroup is main chain itself
        self.subgroups: list[tuple[Entity]] = [tuple(chain.main_chain)]

        # self.insaturations: defaultdict[list[Connection]] = defaultdict(list)  # { Type: [Connections between elements] }
        self.insaturations: list[Connection] = []

        self.classif: dict = {}  # TODO: initalize with some functionals keys (but empty sets)

        self.classif_by_host: dict[int] = {}  # Host: [functions]

        self.host_by_classif: dict[str] = {} # Functions: [host]

        # - Function call -
        # Classificates main chain
        self._root_question(0)
    
        self._define_subgroups()

        self._get_instaturations()
        # Collapses classification
        # For now ignores id 1
        self._collapse_classif()
        # Resolves main class
        self._hclass = self._resolve_hclass()
        if self.get_classif(0) == 'UNRESOLVED':
            self._insert_classif('Chain', 0)
        print("Classifications:", self.classif)

        # Define subfunctions (alkane, alkadiene, ...)
        # NOTE: it depends on the higher classification
        # self.subfunction = self._resolve_subfunction()
        # print("(CLASSIFIER)Subfunction:", self.subfunction)
    
    def _insert_classif(self, _class: str, subg_id: int = -1):
        if subg_id < 0:
            subg_id = len(self.subgroups) - 1
        
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

    # TODO: Make it less awkward...
    def get_hclass(self) -> str:
        if self._hclass:
            return list(self._hclass.keys())[0]
        else:
            return ""
    
    def get_trad_hclass(self, plural: bool = False) -> str:
        """
        Returns translated higher classification (for print purpose)
        
        :param self: Classifier object
        :return: Translated higher class
        :rtype: str
        """
        return list(self._hclass.values())[0] + (PLURAL_STR if plural else '')

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
    
    def get_host_subgs_ids(self, host_id: int) -> tuple[int]:
        # So it don't count the main_chain!
        # Maybe take it out the subgroups...
        return [subg_id for subg_id, subg in enumerate(self.subgroups[1:]) if subg[0] == host_id]
    
    def get_main_insat(self):
        return [insat for insat in self.insaturations if insat.from_id in self.chain.main_path]
    
    # - Chemical logic -
    def _define_subgroups(self):
        # TODO: Cleanup
        # Starts from a point on the main_chains
        cyclical = self.get_classif(0) == 'Cycle'
        self.subgroups.extend(iterate_subpaths(self.chain, cyclical))
                    
        # ID 0 is the main chain!
        for subg_id in range(1, len(self.subgroups)):
            self._root_question(subg_id)

    def _get_instaturations(self):
        # NOTE: It finds instaturations globally not only on the main chain
        self.insaturations = []
        for el in self.chain:
            self.insaturations.extend([con for con in el.cons if con > 1])
        # for el in self.chain:
        #     for con in el:
        #         if con > 1:
        #             self.insaturations[con.type].append(con)

    # - Classification logic - 
    # (tree-like function cascade)
    # I think is better to have an if forest tbh...
    def _root_question(self, subg_id: int):
        # 0. Is it a cycle?
        if (self.subgroups[subg_id][0] == self.subgroups[subg_id][-1]) and\
            (len(self.subgroups[subg_id]) > 1):
            self.main_cyclical = True
            # Tries to see if it's a special cycle
            self._insert_classif(self._resolve_cycle( self.subgroups[subg_id]), subg_id)
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
                self._insert_classif('Radical', subg_id)
    
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
                            self._insert_classif('Nitro', subg_id)
                        else:
                            self._insert_classif('Radical', subg_id)
                    elif subg_id == 0: # Only the main chain can get here
                        # Name in alphabetical order
                        self._insert_classif('Amine', subg_id)

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
                        self._insert_classif('Alcohol', subg_id)
                    if any([(con.to_id == _host.id)\
                             and\
                            (con.type == DOUBLE)\
                            for con in ent.cons]):
                        # print('Is an aldehyde!')
                        self._insert_classif('Aldehyde', subg_id)
                    if (self.chain.to_el(ent.id,'C') == 2):
                       # For now flags it if there is 2 connection to carbons
                       self._insert_classif('Ether', subg_id)
    
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
                        # TODO: Add only halide as a classification maybe?
                        self._insert_classif(ent.el, subg_id)

    # - Main chain -
    def _resolve_cycle(self, subg: list[Entity]):
        # TODO: Should classificate by number of benzenes!
        #       Maybe a function to find and count them?
        _classif = 'Cycle'
        n2 = 0
        for el in subg[1:]:
            if el.cons.count(2) == 1:
                n2 += 1
        if n2 / 2 == 3:  # Double counting!
            _classif = 'Benzene'
        return _classif

    def _collapse_classif(self):
        # Carbonic Acid
        for host_id, cbh in self.classif_by_host.items():
            if ('Aldehyde' in cbh) and\
                ('Alcohol' in cbh):
                # For now maintains the original classifications too
                # Gets relevant subgroups
                subgs_ids = self.get_host_subgs_ids(host_id)
                # Order in this type (collapsed) of group don't matters
                # As long as the host still is the one in postion 0
                new_subg = []
                for subg_id in subgs_ids:
                    host_subg = self.subgroups[subg_id]
                    for ent in host_subg:
                        if ent.id not in new_subg:
                            new_subg.append(ent)
                # new_subg = [self.subgroups[new_subg_id] for new_subg_id in new_subg_ids]
                self.subgroups.append(tuple(new_subg))
                self._insert_classif('Carbonic Acid')

    # TODO: Maybe just return a tuple pair?
    def _resolve_hclass(self) -> dict:
        for _host in self.classif_by_host:
            _host_classif = self.classif_by_host[_host]
            # TODO: Add "Amide"
            # I think a match would be better...
            is_acid = 'Carbonic Acid' in _host_classif
            is_amine = 'Amine' in _host_classif
            is_amide = 'Amide' in _host_classif
            has_alcohol = 'Alcohol' in _host_classif
            has_ether = 'Ether' in _host_classif
            has_keton = 'Keton' in _host_classif
            
            if is_acid:
                return{'Carbonic Acid':CLASSIFICATION['Carbonic Acid']}
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

    # Maybe specify on a conf file how they should be decided?
    # def _resolve_subfunction(self) -> dict:
    #     subfunc_prefix = ''
    #     subfunc_infix = ''
    #     subfunc_suffix = ''
    #     subfunc = SUBCLASSIFICATION[self.get_hclass()]

    #     if self.get_hclass() in SUBCLASSIFICATION:
    #         # For now reads the first if a list
    #         if isinstance(subfunc, list):
    #             subfunc_prefix = subfunc['prefix'][0]
    #         else:
    #             subfunc_prefix = subfunc['prefix']
    #     else:
    #         subfunc_prefix = UNRESOLVED
        
    #     match self.get_hclass():
    #         case 'Hydrocarbon':
    #             subfunc_infix = _name_con_type(self.insaturations, self.chain.get_main_path_id, show_ids=False)
        
    #     subfunc_suffix = SUFFIXES[self.get_hclass()]
    #     final_str = subfunc_prefix + subfunc_infix + subfunc_suffix
    #     return final_str
            

# -- Naming --
# - Function -
# - Helpers -
# Get connections infix
def mult_prefix(ids: list[int] = [],
                show_ids = False,
                trailing_hyphen = True):
    n = len(ids)
    _trailing = '-' if trailing_hyphen else ''
    
    # NOTE: For now assumes all multiple connection have this 'a' added (hardcoded until figured out)
    _mult_prefix_str = 'a' if n > 1 else ''
    
    if n < len(MULT_PREFS):
        if (show_ids) and (n > 0):
            _mult_prefix_str += _trailing + ','.join(map(str, ids)) + '-'
        _mult_prefix_str += MULT_PREFS[n]
    else:
        _mult_prefix_str = UNDEF
    return _mult_prefix_str

get_prefix = lambda n: PREFIXES[n] if n < len(PREFIXES) else UNDEF

# - Main chain naming -
def _name_size_pref(n_main: int) -> str:
    if n_main <= len(PREFIXES):
        return PREFIXES[n_main]
    else:
        return UNDEF


def _name_con_type(cons: list[Connection],
                   on_main: Chain.get_main_path_id,
                   show_ids: bool = True,
                   hide_mult: bool = False) -> str:
    """
    Name connection type (an, en, in, etc.) using the id from the main_path of the chain
    
    :param cons: List of connections to name
    :type cons: list[Connection]
    :param on_main: Function, bound to a Chain instance, to get main path id of the connection
    :type on_main: Chain.get_main_path_id
    :param show_ids: Flag to decide if main path ids are shown or not (default is True) - only valid for main path
    :type show_ids: bool
    :param hide_mult: Flag to decide if main path multiplier prefixes are shown
    :type hide_mult: bool
    :return: Infix string
    :rtype: str
    """
    infix = ''
    cons_dict = defaultdict(list)
    # - Differenciate into unique classes - 
    # TODO: Make it leaner, is a little rough tbh
    to_id_pool: list = []
    for con in cons:
        from_id = con.from_id
        to_id = con.to_id

        if show_ids:
            from_id_main = on_main(con.from_id)
            if from_id_main < 0:
                raise ValueError("From ID outside main path with shoud_ids set to True!")
        else:
            from_id_main = con.from_id
        # Ideally only valid values would get to this point (controlling by maximum connections)
        if con.type > 1 and con.type < len(INFIXES):
            if from_id in to_id_pool:
                # To avoid possible duplicates
                # The "from" [0] on the connection should be different
                # from any "to" [1] already in cons_dict
                continue
            cons_dict[con.type].append(from_id_main)
            to_id_pool.append(to_id)
    # Sorts keys in ascending order
    cons_dict = dict(sorted(cons_dict.items(), key = lambda e: INFIXES[e[0]]))
    #  - Make infix - 
    # Positions
    for type in cons_dict:
        # Adds ids
        # A little bit unoptmized tbh, but leaner
        cons_from = [con_from for con_from in cons_dict[type]]
        if not hide_mult:
            infix += mult_prefix(cons_from, show_ids = show_ids)
    # Addition of type name
    for type in cons_dict:
        if type < len(INFIXES) - 1:
            infix += INFIXES[type]
        else:
            infix += UNDEF
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


def _name_suffix(classific: Classifier, hide_ids = False, hide_mult = False) -> str:
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

    prefix = infix = suffix = ''
    # Clasificates chain
    # TODO: Make it less redundant... Two ways to classify huh???
    # chain.functional = _get_class(chain)
    _classfier = Classifier(chain)
    
    _is_cyclic = _classfier.main_cyclical
    _is_aromatic = 'Benzene' in _classfier.classif  # For now only this!
    _is_single_atom = len(chain.main_chain) == 1
    _is_nitrogenic = _classfier.get_hclass() == 'Nitrogenic'  # TODO: rename it, as it's now, it includes Nitros!
    _is_oxygenic = _classfier.get_hclass() == 'Ester' or _classfier.get_hclass() == 'Ether'
    # Strictly speaking, if it can only occur at end of the chain it shoukd be numbered!
    # TODO: Make it now if should be numbered another way!
    _is_acid = _classfier.get_hclass() == 'Carbonic Acid'

    # TODO: Reafctor to or HIDE or SHOW
    _hide_prefix_id = _is_single_atom
    _infixed = not _is_nitrogenic
    _show_infix_id = not _is_aromatic
    _hide_suffix_id = _is_oxygenic or _is_acid or _is_aromatic
    _hide_mult = _is_aromatic
    _is_afixed = _is_cyclic or _is_nitrogenic

    prefix = _name_substitutive(_classfier, hide_ids=_hide_prefix_id)

    # If is a cycle (main chain)
    # Maybe add null values for afixes of other classes
    # Afixables
    if _is_afixed:
        prefix += AFIXES[_classfier.get_classif(0)]
    else:
        prefix += _name_size_pref(len(set(chain.main_path)))

        # TODO: better understand when to use the hyphen!

        if not prefix:
            return "Chain is too big!"

        # Makes chain
        cons = []
        # This could be strealined to use connection class!
        old_pos = None
        for pos in [chain.id_pool[id] for id in chain.main_path]:
            if old_pos:
                # This would be perfect in a function!
                dif_pos = pos - old_pos
                norm_dif = dif_pos/2
                con_pos = old_pos + norm_dif
                cons.append(chain.field[con_pos.row][con_pos.col])
            old_pos = pos
    
    if _infixed:
        infix = _name_con_type(_classfier.get_main_insat(),
                               _classfier.chain.get_main_path_id,
                               show_ids=_show_infix_id,
                               hide_mult=_hide_mult)

    suffix += _name_suffix(_classfier, hide_ids=_hide_suffix_id, hide_mult=_hide_mult)
    
    chain.name = prefix + infix + suffix

    # Make this better!
    if _is_acid:
        if REVERSE_STR:
            chain.name = AFIXES['Acid'] + chain.name
        else:
            chain.name += AFIXES['Acid']
    
    chain.func_name = "NOT IMPLEMENTED!"
    # TODO: Cyclan, Cyclin and whatnot?
    # I don't like it. I think is best to separete it!
    # TODO: needs to add 'acid'to chain!
    if REVERSE_STR:
        chain.functional = _classfier.get_trad_classif(0) + CON_STR + _classfier.get_trad_hclass(plural=True)
    else:    
        chain.functional = _classfier.get_trad_hclass(plural=True) + CON_STR + _classfier.get_trad_classif(0)
