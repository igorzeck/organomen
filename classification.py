from base_structures import Pos
from pathfinding import scout
from entities import Chain

# TODO: Here might be better to pass pos as list of POS instead of ids

# -- Naming --
# - Name prefix -
# - Constants -
HETEROATOMS = [
    'F',
    'O',
    'N'
]

# Relative to N of Cs in main chain
PREFIXES = [
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
# Relative to chain functional class
SUFIXES = [
    'o', # Hydrocarbon
    'ol', # Alcohol
    'al', # Aldehyde
    'oico', # Acid
    'ona', # Keton
]
# For now infixes are not included
# - Function -
# - Main chain naming -
def _name_size_pref(n_main: int):
    if n_main <= len(PREFIXES):
        return PREFIXES[n_main - 1]
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
    # TODO: Change it so it's only between carbon to carbon connections
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


# Will be more complicated for cyclic types...
def _name_functional(field, ids: dict[Pos], main_chain: list):
    # atoms_set = set()
    # atoms_list = []
    atoms_dict = {}
    for id in main_chain:
        pos = ids[id]
        el = field[pos.row][pos.col]
        if el in atoms_dict:
            atoms_dict[el].append(id)
        else:
            atoms_dict[el] = [id]
        # atoms_set.add(el)
        # atoms_list.append((el, pos))
    # Checks heteroathoms
    if len(atoms_dict) == 1:
        return 'o'
    # Separated function?
    # For now by hand
    # TODO: Separate into functions as its possible to have multiple functional groups at the same time!
    if 'O' in atoms_dict:
        # Checks postion of Oxygen(s)
        # Alcohol
        for pos_id in atoms_dict['O']:
            cons_list = scout(field, ids, pos_id, nxt_dir=False, con_value=True)
            if len(cons_list) >= 2:
                # Assumes ceton for now
                pass  # Ceton code...
            else:  # Sees connection type
                nxt_id, con_value = cons_list[0]
                nxt_pos = ids[nxt_id]
                print(field[nxt_pos.row][nxt_pos.col], con_value)
                if con_value == '1':
                    # can be acid too!
                    # TODO: Needs to check next pos to confirm
                    return 'ol'
                elif con_value == '2':
                    return 'al'
    return 'o'


# TODO: Make id go over the entire chain, not only main one
def _classif_atom(field: list[Pos], ids, pos_id):
    pos = ids[pos_id]
    el = field[pos.row][pos.col]
    info = scout(field, ids, pos_id)
    els = [field[ids[id[0]].row][ids[id[0]].col] for id in info]
    
    
    # Oxygen
    if el != 'O':
        return 'Unsuported!'
    
    # 1. How many Cs
    if els.count('C') > 1:
        return 'Ether'

    # Gets host info
    host_info = []
    breakout = False
    for id in info[0]:
        nxt_info = scout(field, ids, id)
        for nxt_id in nxt_info:
            nxt_el = field[ids[nxt_id[0]].row][ids[nxt_id[0]].col]
            if nxt_el == 'C':
                host_info = scout(field, ids, id, con_value=True)
                breakout = True
                break
        if breakout:
            break
    
    host_els = [field[ids[id].row][ids[id].col] for id, _ in host_info]
    
    # 2. Is it a edge C?
    c_edge = host_els.count('C') == 1

    # 3. Does the host connect to other Os?
    c_alone = host_els.count('O') == 1

    # 4. How many double connections with Os?
    c_n_double = sum([con == '2' and field[ids[id].row][ids[id].col] == 'O' for id, con in host_info])
    
    # To be honest I'm not happy with this, but it works
    if c_n_double == 1:
        if c_edge:
            if c_alone:
                return 'Aldehyde'
            else:
                return 'Acid'
        if c_alone:
            return 'Ceton'
    elif c_n_double == 0:
        return "Alchol"
    # 4. 
    return 'Unknown'


def _get_class(chain: Chain):
    functionals = []

    field = chain.field
    ids = chain.id_dict
    main_path = chain.main_path
    
    atoms_dict = {}
    for id in ids:
        pos = ids[id]
        el = field[pos.row][pos.col]
        if el in atoms_dict:
            atoms_dict[el].append(id)
        else:
            atoms_dict[el] = [id]
    print(atoms_dict)
    # Decides in its functional class
    if any(el == 'C' for el in atoms_dict):
        if all(el == 'C' for el in atoms_dict):
            functionals.append("hydrocarbon")
            return functionals
        else:
            atoms_dict.pop('C')
    else:
        functionals.append("not organic")
        return functionals
    
    # Heteroatoms
    for id in ids:
        pos = ids[id]
        el = field[pos.row][pos.col]
        if el in atoms_dict:
           functionals.append(_classif_atom(field, ids, id))

    return functionals


def class_chain(chain: Chain):
    field = chain.field
    ids = chain.id_dict
    main_path = chain.main_path

    if not main_path:
        return "Chain is empty!"
    prefix = _name_size_pref(len(main_path))
    if not prefix:
        return "Chain is too big!"
    
    # Clasificates chain
    chain.functional = _get_class(chain)

    # Makes chain
    cons = []
    old_pos = None
    for pos in [ids[id] for id in main_path]:
        if old_pos:
            # This would be perfect in a function!
            dif_pos = pos - old_pos
            norm_dif = dif_pos/2
            con_pos = old_pos + norm_dif
            cons.append(field[con_pos.row][con_pos.col])
        old_pos = pos
    infix = _name_con_type(cons)
    sufix = _name_functional(field, ids, main_path)
    
    # TODO: Hopefully temporary
    return prefix + infix + sufix
