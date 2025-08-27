from copy import deepcopy
from base_structures import *
from auxiliary import *
from entities import Chain, Entity, Connection

# - Runners -
# Should first selected based on: heteroatom
# Then: Connections
# Then: size
# In this order!
def _get_host(chain: Chain, ent: Entity):
    """
    Return the host - atom inside main path - of the entity
    If it doesn't have an host, return itself
    """
    return run_subpath(chain, ent)


def run_subpath(chain: Chain, ent: Entity, subpath: list[int], con: Connection = None):
    """
    Explore all connections directions
    Stops if End of Path EOP, that is:
    - Finds itself (in a loop)
    - If finds an edge
    - If finds an main_chain atom (return it)
    subpath is changed by reference
    """
    nxt_ent = ent
    failsafe = 0
    while True:
        subpath.append(ent)
        # Found main path
        if nxt_ent in chain.main_chain:
            print("FIM")
            return subpath
            # return curr_ent
        # First iteration
        if not con:
            for _con in ent:
                subpath = run_subpath(chain, chain.chain[_con.to_id], subpath, _con)
            return subpath
        else:
            nxt_ent = chain.chain[con.to_id]
        # Edge-case
        if len(nxt_ent) == 1:
            return subpath
        # Recursion
        elif len(nxt_ent) > 2:
            for _con in nxt_ent.cons:
                if _con != con:
                    subpath = run_subpath(chain, nxt_ent, subpath, _con)
        # Cyclic-case
        if nxt_ent.id in subpath:
            return subpath
        if failsafe > 100:
            raise IndexError("Too big of a chain: Too big subapath, malformed chain!")



def _is_higher(best: list, contender: list):
    """
    Checks if path is best fit to be the main one.
    """
    if len(contender) > len(best):
        return True
    else:
        return False


def run_path(chain: Chain,
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
        print_field(chain.field, [chain.id_pool[pos_id]])
        nxt_els = chain.chain[pos_id]
        nxt_n_con = len(nxt_els.cons)
        print(nxt_els, pos_id)
        # Cycle detection
        if pos_id in path_stub:
            # It will have repeated position
            # But it will make it easier to glance at
            # Cyclic behaviour!
            path_stub.append(pos_id)
            return path_stub
        path_stub.append(pos_id)

        # Lone heteroatom detection
        # Kinda hacky tbh
        pos = chain.id_pool[pos_id]
        el = chain.field[pos.row][pos.col]
        if el != 'C':
            if len(scout(chain.field, chain.id_pool,pos_id)) == 1:
                return []

        if (nxt_n_con == 2 and pos_id not in chain.edges) or pos_id == start_pos_id:
            # If it has one possible path, move
            jump = False
            for _con in nxt_els:
                nxt_dir = _con.dir
                nxt_el_pos_id = _con.to_id
                if nxt_dir != origin_dir:
                    origin_dir = -nxt_dir
                    pos_id = nxt_el_pos_id
                    jump = True
                    break
            if jump:
                continue
        if nxt_n_con > 2 and pos_id not in chain.edges:
            # Recursion for multiple paths
            for _con in nxt_els:
                nxt_dir = _con.dir
                nxt_el_pos_id = _con.to_id
                if nxt_dir != origin_dir:
                    print(f"({recur + 1}) Going ({nxt_dir}, {nxt_el_pos_id})")
                    print(best_stub)
                    temp_stub = run_path(chain, nxt_el_pos_id, origin_dir=-nxt_dir, recur = recur + 1, path_stub = deepcopy(path_stub), best_stub=best_stub)
                    if _is_higher(best_stub, temp_stub):
                        best_stub = temp_stub
            break
        if pos_id in chain.edges:
            print(f"({recur}) Done")
            print(*path_stub)
            # chain.add_path(path_stub)
            if _is_higher(best_stub, path_stub):
                best_stub = path_stub
                chain.main_path = best_stub
            print_field(chain.field, [chain.id_pool[id] for id in path_stub])
            break
    return best_stub


# Obsolet
def get_sub_groups(chain: Chain):
    # Get substitutive groups
    Tto_highlight: Pos = []
    groups: list[list[Entity]] = []
    for ent in chain.main_chain:
        # Iterates through all connections outside main_chain
        for con in ent.cons:
            if con.to_id not in chain.main_chain:
                # Follows lead
                groups.append(follow_path(chain.chain, con.to_id, con.dir))
                print(groups)

        if any([con.to_id not in chain.main_path for con in ent.cons]):
            Tto_highlight.append(chain.id_pool[ent.id])
            # print(ent)
    print(Tto_highlight)
    print_field(chain.field, Tto_highlight)


def follow_path(chain_path: list[Entity], to_id: int, dir: int):
    ent = chain_path[to_id]
    pos_id = to_id
    while len(ent) > 1 and pos_id not in chain_path:
        print(ent)
        return ent


# TODO: Insert chain as a parameter
# TODO: It needs to start and end at an edge!
def run_chain(field):
    """
    Go through the entirety of a chain field representation
    """
    chain = Chain(field)
    for pos_id in chain.edges:
        print(f"Start from {chain.id_pool[pos_id]}")
        chain.main_path = run_path(chain,
                                pos_id,
                                path_stub=[],
                                best_stub=chain.main_path
                                )
    if chain.main_path:
        print("Longest:")
        print_field(field, [chain.id_pool[id] for id in chain.main_path])
    else:
        print("Empty field!") 
    # ~ Test ~
    for pos_id in chain.main_path:
        chain.main_chain.append(chain.chain[pos_id])

    # chain.groups = get_sub_groups(chain)
    get_sub_groups(chain)

    return chain  