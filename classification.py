from base_structures import Pos
from pathfinding import scout
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
            cons_list = scout(field, ids, pos_id, direction=False, con_value=True)
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
    

def name_chain(field, ids: dict[Pos], all_chains: list, main_chain: list):
    if not main_chain:
        return "Chain is empty!"
    prefix = _name_size_pref(len(main_chain))
    if not prefix:
        return "Chain is too big!"
    
    # Makes chain
    cons = []
    old_pos = None
    for pos in [ids[id] for id in main_chain]:
        if old_pos:
            # This would be perfect in a function!
            dif_pos = pos - old_pos
            norm_dif = dif_pos/2
            con_pos = old_pos + norm_dif
            cons.append(field[con_pos.row][con_pos.col])
        old_pos = pos
    infix = _name_con_type(cons)
    sufix = _name_functional(field, ids, main_chain)
    
    # TODO: Hopefully temporary
    return prefix + infix + sufix
