from copy import deepcopy
from base_structures import *
from auxiliary import *
from entities import Chain, Entity, Connection

# - Runners -
# Should first selected based on: heteroatom
# Then: Connections
# Then: size
# In this order!
def _get_host(main_chain: list[Entity], ent: Entity):
    """
    Return the host - atom inside main path - connected to the entity
    If it doesn't have an host, return itself
    """
    if ent in main_chain:
        return ent
    for el in main_chain:
        if any([el.id == con.to_id for con in ent.cons]):
            return el
    # Returns itself if no host in main path (e.g. if "host" on subgroup)
    return ent


# TODO: Go to all ids not already in subpath list
# If instead of recursion I use a stack of ids (FIFO)?
def run_subpath(chain: Chain, ent: Entity):
    """
    Explore all connections directions
    With stack of ids
    Stops if End of Path EOP, that is:
    - Finds itself (in a loop)
    - If finds an edge
    - If finds an main_chain atom (return it)
    
    :param chain: Chain to be expored
    :type chain: Chain
    :param ent: Entity where the subpath starts
    :type ent: Entity
    """
    failsafe = 0
    ents_pool: list[Entity] = [ent]
    queue: list[Entity] = []  # Stack with ids by order
    while True:
        # Updates queue
        for con in ent:
            nxt_ent = chain[con.to_id]
            if nxt_ent not in ents_pool:
                if nxt_ent in chain.main_chain:
                    # For now assumes its the host and put it on 0 index!
                    # TODO: handle multiple hosts
                    ents_pool.insert(0,nxt_ent)
                    continue
                else:
                    ents_pool.append(nxt_ent)
                    queue.append(nxt_ent)
        if not queue:
            return tuple(ents_pool)
        # Edge-case
        # Goes FIFO through queue
        ent = queue.pop()

        failsafe += 1
        if failsafe > 100:
            raise IndexError("Too big of a chain: Too big subapath, malformed chain!")


# TODO: For same size subpaths, the ones with the most ramifications should be the main one
# TODO: The lower position number should fallback to alphabetical order when equivalent positions
def _is_higher(best: list, contender: list, chain: list[Entity]):
    """
    Checks if path is best fit to be the main one.
    """
    # TODO: order this to check on proper order
    # - Checks for empty values -
    if not contender:
        return False
    
    # - Checks for early exit -
    if not best:
        return True
    # - Check variables -
    contender_hetero = False
    best_hetero = False

    # O ideal é contar o número de insaturações!
    contender_insat = 0
    best_insat = 0
    
    # For the most groups and the closest to start
    contender_group = []
    best_group = []

    # - Checks -
    for id in contender:
        curr_el = chain[id]
        # 1. Heteroatoms
        if curr_el.is_hetero():
            contender_hetero = True
        # 2. Insaturation
        if any([con != SIMPLE for con in curr_el.cons if chain[con.to_id] == 'C']):
            # As it's a linear path I can add 1 whenever it finds a instaturation regardless of how many!
            contender_insat += 1
        # 3. Closest group (careful with cyclical logic)
        # Looks to more than 3 connections to carbons
        # Note that any group should flag this, not onlu substitutive groups
        if len(curr_el) > 2:
            contender_group.append(contender.index(id))

    for id in best:
        curr_el = chain[id]
        # 1. Heteroatoms
        if curr_el.is_hetero():
            best_hetero = True
        # 2. Insaturation
        if any([con != SIMPLE for con in curr_el.cons if chain[con.to_id] == 'C']):
            best_insat += 1
        # 3. Closest group
        if len(curr_el) > 2:
            best_group.append(best.index(id))
    # - Comparisons -
    # 1. Heteroatoms (Obsolete?)
    if contender_hetero:
        if not best_hetero:
            return True
    elif best_hetero:
        return False
    
    # 2. Insaturation
    # Does nothing if they both are equal!
    if contender_insat > best_insat:
        return True
    elif contender_insat < best_insat:
        return False

    # 3. Group
    if contender_group and best_group:
        if len(contender_group) > len(best_group):
            # Contender has more groups 
            return True
        elif len(contender_group) == len(best_group):
            # TODO: Checks to see if, in alphabetical order, the lowets is closest
            # Checks the one with the closest group
            if min(contender_group) < min(best_group):
                return True
        else:
            # Current best has more groups
            False
    
    # 4. Size comparison
    if len(contender) > len(best):
        return True
    else:
        return False
    # If both are equal don't mtter which is returned


# TODO: Make it follows all rules for defining main chain!
# TODO: Make it not start from a heteroatom
def run_path_recursive(chain: Chain,
             start_pos_id: int,
             path_stub: list,
             best_stub: list,
             origin_dir = 0,
             recur = 0):
    """
    Go through the 2D representation of a carbon chain
    path
    recur: Current recursion
    origin_dir: Previous mirror direction to avoid backtracking
    """
    # Start position of this recursion
    pos_id = start_pos_id
    
    # Temporary as it doesn't work with recursion!
    while True:
        print_field(chain.field, [chain.id_pool[pos_id]])
        nxt_els = chain.chain[pos_id]
        nxt_n_con = len(nxt_els.cons)
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

        # if (nxt_n_con == 2 and pos_id not in chain.edges) or pos_id == start_pos_id:
        if (nxt_n_con <= 2):
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
                    print(f"({recur + 1}) Going (Dir: {nxt_dir}, ID: {nxt_el_pos_id})")
                    temp_stub = run_path_recursive(chain, nxt_el_pos_id, origin_dir=-nxt_dir, recur = recur + 1, path_stub = deepcopy(path_stub), best_stub=best_stub)
                    if _is_higher(best_stub, temp_stub, chain.chain):
                        best_stub = temp_stub
                    print_field(chain.field, [chain.id_pool[id] for id in path_stub])

            break
        if pos_id in chain.edges:
            print(f"({recur}) Done")
            if _is_higher(best_stub, path_stub, chain.chain):
                best_stub = path_stub
                chain.main_path = best_stub
            print_field(chain.field, [chain.id_pool[id] for id in path_stub])
            break
    return best_stub


# TODO: Oxygen should go to main chain except if on edge! (So does Nitrogen)
# Maybe they should BE the main_chain (registered as edge if connection to multiple Cs?)
def run_path_iterative(chain: Chain):
    best_path: list[int] = []

    for start_id in chain.edges:
        curr_el: Entity = chain.chain[start_id]
        curr_path: list[int] = [start_id]
        action_stack: list[tuple] = []

        print(f"Starting at edge {start_id}")
        print_field(
            chain.field,
            highlights=[chain.id_pool[start_id]]
        )

        # End of Path flag
        eop = False

        # Appends initial direction on action stack
        for con in curr_el:
            action_stack.append((len(curr_path), con.to_id))

        # Follows action_stack
        while action_stack:
            # Goes to the next unless is an Nitrogen edge
            if not curr_el.el == 'N':
                curr_path_id, curr_el_id = action_stack.pop()
                curr_el = chain.chain[curr_el_id]
                print_field(
                    chain.field,
                    highlights=[chain.id_pool[curr_el_id]]
            )
            else:
                action_stack.clear()

            # - End of Path check -
            # At edge or already on the path (cyclical)
            eop = (curr_el in chain.edges) or\
                    (curr_el_id in curr_path)
            
            if curr_path_id < len(curr_path):
                curr_path = curr_path[:curr_path_id]
            
            # Such that heteroatoms are not on the main chain
            if curr_el.el not in HETEROATOMS:
                curr_path.append(curr_el_id)
            else:
                eop = True

            if eop:
                # Checks if curr_path is a better main path than best_path
                if _is_higher(best_path, curr_path, chain):
                    # New best_path
                    best_path = curr_path
                print("Full path:")
                print_field(
                    chain.field,
                    highlights=[chain.id_pool[id] for id in curr_path]
                )
            else:
                for con in curr_el:
                    if con.to_id not in curr_path:
                        action_stack.append((len(curr_path), con.to_id))


    return best_path


def get_sub_groups(chain: Chain):
    # Get substitutive groups
    to_highlight: Pos = []
    groups: list[list[Entity]] = []
    for ent in chain.main_chain:
        # Iterates through all connections outside main_chain
        for con in ent:
            if con.to_id not in chain.main_path:
                # Follows lead
                groups.append(follow_subpath(chain.chain, chain.main_path, [], con.dir, con.to_id))
                print("Captured group:\n",groups)

        if any([con.to_id not in chain.main_path for con in ent.cons]):
            to_highlight.append(chain.id_pool[ent.id])
    print_field(chain.field, to_highlight, highlight_color_only=True)


# TODO: Generalize as the main pathfinder?
def follow_subpath(chain: list[Entity], main_path: list[Entity], group: list[Entity], origin_dir: int, curr_id: int, cont: int = 0):
    """
    Follow subpath recursively (at each junction) until
    all entities outside of main_path are in group
    
    :param chain: Chain with all entities
    :type chain: list[Entity]
    :param main_path: Main path (acts as mask to keep on group)
    :type main_path: list[Entity]
    :param group: Group of unique entities on this. Needs to be 
    inserted initally as a empty list for safe-copy reasons
    :type group: list[Entity]
    :param dir: Origin direction (forbidden on first iteration)
    :type dir: int
    :return: Returns final group
    :rtype: list[Entity]
    """
    print("Following subpath!")
    ent = chain[curr_id]

    dir = -origin_dir

    _cycled = ent in group
    _back_at_main = ent in main_path

    while not _cycled and not _back_at_main:
        # Failsafe
        cont += 1
        if cont >= MAX_ITER:
            break

        group.append(ent)

        # - Advance -
        if len(ent) == 2:
            # If there is two directions, go to the not-origin direction
            for con in ent:
                if con.dir != dir:
                    to_id = ent.at_dir(con.dir)
                    dir = -con.dir
        elif len(ent) > 2:
            # Recursion for the conjunction
            for con in ent:
                if con.dir != dir:
                    new_ent_id = ent.at_dir(con.dir)
                    group = follow_subpath(chain, main_path, group, con.dir, new_ent_id, cont)
            break
        elif len(ent) <= 1:
            # End of path
            break

        # - Additional finishing conditions -
        # Id not in chain (-1)
        if to_id >= 0:
            ent = chain[to_id]
        else:
            break

        _cycled = ent in group
        _back_at_main = ent in main_path


    return group


# TODO: Insert chain as a parameter
# TODO: It needs to start and end at an edge!s
def run_chain(field):
    """
    Go through the entirety of a chain field representation
    """
    chain = Chain(field)
    # for pos_id in chain.edges:
    #     print(f"Start from {chain.id_pool[pos_id]}")
    #     # With path function
    #     chain.main_path = run_path_recursive(chain,
    #                             pos_id,
    #                             path_stub=[],
    #                             best_stub=chain.main_path
    #                             )
    chain.main_path = run_path_iterative(chain)
    if chain.main_path:
        print("Main chain:")
        print_field(field, [chain.id_pool[id] for id in chain.main_path])
    else:
        print("Empty field!") 
    # Elements info
    for pos_id in chain.main_path:
        chain.main_chain.append(chain.chain[pos_id])

    # get_sub_groups(chain)

    return chain
