from copy import deepcopy
from base_structures import *
from entities import Chain

# - Pathfinder -
# Valid atom: Not Hydrogen basically
# TODO: Change to only detect Carbon to Carbon connections on chain end!
def _edge_c(field, ids):
    """
    Gives the postion of the of an edge Carbon
    field: Flat representation of the chain
    ids: Dictionary with carbon id and its postion
    """
    curr_id = None
    edgeless = True
    for id in ids:
        pos = ids[id]
        el = field[pos.row][pos.col]
        if el == 'C':
            curr_id = id
            if len(scout(field, ids, id)) == 1:
                if edgeless:
                    edgeless = False
                yield curr_id
    if curr_id and edgeless:
        # If finds no edge, return last curr_pos
        yield curr_id


def _is_inside(dir: Pos, n_row, n_col):
    if (dir.row >= n_row or dir.col >= n_col)\
        or\
        (dir.row < 0 or dir.col < 0):
        return False
    return True


# Should first selected based on: heteroatom
# Then: Connections
# Then: size
# In this order!
def _is_higher(best: list, contender: list):
    """
    Checks if path is best fit to be the main one.
    """
    if len(contender) > len(best):
        return True
    else:
        return False


# TODO: Version to accept positions AND id positions
def scout(field, ids: dict[Pos], pos_id: int, nxt_dir = False, con_value = False):
    """
    Scouts the vicinity of an given position
    field: Chain field
    ids: Chain dict containing a table of ids to postion in the field
    pos_id: Pos id to investigate around
    direction: If True returns direction
    con_value: If true returns connection value
    returns: list of a tuple (package) with up to the following elements (dir, next_id, connection)
    """
    nxt_el_l = []
    for dir in cd_offsets:
        pos = ids[pos_id]
        # Next direction and connectin (to scout)
        nxt_step = cd_offsets[dir]
        nxt_con_pos = pos + nxt_step
        nxt_el_pos = nxt_con_pos + nxt_step
        if _is_inside(nxt_el_pos, len(field), len(field[0])):
            # If is inside, check to see if there is a valid atom in the next path
            nxt_el = field[nxt_el_pos.row][nxt_el_pos.col]
            next_con = field[nxt_con_pos.row][nxt_con_pos.col]
            if next_con.isnumeric() and next_con != '0'\
                and\
                not nxt_el.isnumeric():
                # Kinda dumb... But it wil work...
                selected_id = -1
                for id in ids:
                    _pos = ids[id]
                    if _pos == nxt_el_pos:
                        selected_id = id
                if selected_id < 0:
                    raise ValueError(f"Didn't find position {nxt_el_pos} id's")
                # print("Scout:",dir, nxt_el_pos)
                # TODO: Reorder this bit of code so selected id is the first
                _package = []
                if nxt_dir:
                    _package.append(dir)
                _package.append(selected_id)
                if con_value:
                    _package.append(next_con)
                nxt_el_l.append(_package)       
    return nxt_el_l


# TODO: enroll all variables in a chain class
def per_path(chain: Chain,
             start_pos_id: int,
             path_stub: list,
             best_stub: list,
             origin_dir = 0,
             recur = 0):
    """
    Go through the 2D representation of a carbon chain
    path
    recur: Current recursion (TODO: take it out)
    origin_dir: Previous mirror direction to avoid backtracking
    """
    # Start position of this recursion
    pos_id = start_pos_id
    # Temporary as it doesn't work with recursion!
    while True:
        print_field(chain.field, [chain.id_dict[pos_id]])
        nxt_els = scout(chain.field, chain.id_dict, pos_id, nxt_dir=True)
        nxt_n_con = len(nxt_els)
        
        # Cycle detection
        if pos_id in path_stub:
            # It will have repeated position
            # But it will make it easier to glance at
            # Cyclic behaviour!
            path_stub.append(pos_id)
            return path_stub
        path_stub.append(pos_id)

        if nxt_n_con == 2 or pos_id == start_pos_id:
            # If it has one possible path, move
            jump = False
            for nxt_dir, nxt_el_pos_id in nxt_els:
                if nxt_dir != origin_dir:
                    origin_dir = -nxt_dir
                    pos_id = nxt_el_pos_id
                    jump = True
                    break
            if jump:
                continue
        if nxt_n_con > 2:
            # Recursion for multiple paths
            for nxt_dir, nxt_el_pos_id in nxt_els:
                if nxt_dir != origin_dir:
                    print(f"({recur + 1}) Going ({nxt_dir}, {nxt_el_pos_id})")
                    print(best_stub)
                    temp_stub = per_path(chain, nxt_el_pos_id, origin_dir=-nxt_dir, recur = recur + 1, path_stub = deepcopy(path_stub), best_stub=best_stub)
                    if _is_higher(best_stub, temp_stub):
                        best_stub = temp_stub
            break
        if nxt_n_con == 1:
            print(f"({recur}) Done")
            print(*path_stub)
            chain.add_path(path_stub)
            if _is_higher(best_stub, path_stub):
                best_stub = path_stub
                chain.main_path = best_stub
            print_field(chain.field, [chain.id_dict[id] for id in path_stub])
            break
    return best_stub


# TODO: Insert chain as a parameter
def per_chain(field):
    """
    Go through the entirety of a chain field representation
    """
    chain = Chain(field)
    for id in _edge_c(chain.field, chain.id_dict):
        print(f"Start from {chain.id_dict[id]}")
        chain.main_path = per_path(chain,
                                id,
                                path_stub=[],
                                best_stub=chain.main_path
                                )
    if chain.main_path:
        print("Longest:")
        print_field(field, [chain.id_dict[id] for id in chain.main_path])
    else:
        print("Empty field!") 
    # return unique_paths, bestest_path   
    return chain  
