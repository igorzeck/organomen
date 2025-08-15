from copy import deepcopy
from base_structures import *
from pathfinder import *
from entities import Chain

# - Runners -
# Valid atom: Not Hydrogen basically
def _edge_c(field: list, pos_pool: tuple[int]):
    """
    Gives the postion of the of an edge Carbon
    field: Flat representation of the chain
    ids: Tuple with carbon id and its postion
    """
    curr_id = None
    edgeless = True
    for id_pos, pos in enumerate(pos_pool):
        el = field[pos.row][pos.col]
        if el == 'C':
            curr_id = id_pos
            el_info = scout(field, pos_pool, id_pos, nxt_id=False, nxt_str=True)
            if el_info.count('C') == 1:
                if edgeless:
                    edgeless = False
                yield curr_id
    if curr_id and edgeless:
        # If finds no edge, return last curr_pos
        yield curr_id


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


def run_path_test(chain: Chain,
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
                    temp_stub = run_path_test(chain, nxt_el_pos_id, origin_dir=-nxt_dir, recur = recur + 1, path_stub = deepcopy(path_stub), best_stub=best_stub)
                    if _is_higher(best_stub, temp_stub):
                        best_stub = temp_stub
            break
        if pos_id in chain.edges:
            print(f"({recur}) Done")
            print(*path_stub)
            chain.add_path(path_stub)
            if _is_higher(best_stub, path_stub):
                best_stub = path_stub
                chain.main_path = best_stub
            print_field(chain.field, [chain.id_pool[id] for id in path_stub])
            break
    return best_stub



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
        nxt_els = scout(chain.field, chain.id_pool, pos_id, nxt_dir=True)
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
                    temp_stub = run_path(chain, nxt_el_pos_id, origin_dir=-nxt_dir, recur = recur + 1, path_stub = deepcopy(path_stub), best_stub=best_stub)
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
            print_field(chain.field, [chain.id_pool[id] for id in path_stub])
            break
    return best_stub


# TODO: Insert chain as a parameter
# TODO: It needs to start and end at an edge!
def run_chain(field):
    """
    Go through the entirety of a chain field representation
    """
    chain = Chain(field)
    for pos_id in chain.edges:
        print(f"Start from {chain.id_pool[pos_id]}")
        chain.main_path = run_path_test(chain,
                                pos_id,
                                path_stub=[],
                                best_stub=chain.main_path
                                )
    if chain.main_path:
        print("Longest:")
        print_field(field, [chain.id_pool[id] for id in chain.main_path])
    else:
        print("Empty field!") 
    # return unique_paths, bestest_path   
    return chain  