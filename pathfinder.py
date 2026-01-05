from copy import deepcopy
from base_structures import *
from auxiliary import *
from entities import Chain, Entity, Connection

# - Runners -
# TODO: Correct cyclical naming
# So, if I only look through edges (if any) how do i correctly
# define the main chain of a cyclical compound with an edge from
# a radical? Something to do like: If detects cyclical behaviour
# For comparision sake treat as a main chain in its own right, perhaps?

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


def iterate_subpaths(chain: Chain,
                     main_cyclical:bool = False) -> list[tuple[Entity]]:
    """
    Iterates and captures all subgroups connected to chain.main_chain
    
    :param chain: Chain object with main_chain correctly set
    :type chain: Chain
    :param cyclical: If it should iterate as a cyclical chain
    :type cyclical: bool
    :return: Subgroups (each group as a tuple) list of the captured subgroups
    :rtype: list[tuple[Entity]]
    """
    subgroups = []
    cyclicals = []
    for main_ent in chain.main_chain[(1 if main_cyclical else 0):]:
        for nxt_con in main_ent.cons:
            nxt_ent = chain.chain[nxt_con.to_id]
            if (nxt_ent not in chain.main_chain):
                _subgroup = run_subpath(chain, nxt_ent)
                
                subgroups.append(_subgroup)

                print(f"Captured subgroup ({len(subgroups)}):")
                print(_subgroup)
                print_field(chain.field,
                            [chain.id_pool[el.id] for el in _subgroup],
                            highlight_color_only=True)
    return subgroups


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

    # - Cyclical -
    max_cycle_size = 0
    max_cycle_junct_id = -1
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
                    queue.append(nxt_ent)
                    ents_pool.append(nxt_ent)
            else:
                # TODO: take into account the possibility of multiple cycles
                #       having the same junctions (e.g. Antracene)
                # It may flag multiple cycle points, in this case, the farthest apart is the canonical
                nxt_ent_pool_id = ents_pool.index(nxt_ent)
                curr_cycle_size = len(ents_pool) - nxt_ent_pool_id
                if curr_cycle_size > 3:
                    if curr_cycle_size > max_cycle_size:
                        max_cycle_size = curr_cycle_size
                        max_cycle_junct_id = nxt_ent_pool_id
        if not queue:
            # Appends entry point for the group
            if max_cycle_junct_id >= 0:
                # NOTE: Assumes that if it sees the same
                #       id is within a cycle if there is at least on carbon in between
                ents_pool.append(ents_pool[max_cycle_junct_id])
            # Cheks if last id was nxt_ent
            return tuple(ents_pool)
        # Edge-case
        # Goes FIFO through queue
        ent = queue.pop()

        failsafe += 1
        if failsafe > 100:
            raise IndexError("Too big of a chain: Too big subapath, malformed chain!")

def has_hetero(main_ents: Entity | list[Entity] | tuple[Entity]) -> bool:
    if isinstance(main_ents, Entity):
        return main_ents.is_hetero()
    elif isinstance(main_ents, list | tuple):
        return any([ent.is_hetero() for ent in main_ents])
    else:
        return False

# TODO: Maybe presaving best chain and using multiple cores
# TODO: For same size subpaths, the ones with the most ramifications should be the main one
# TODO: The lower position number should fallback to alphabetical order when equivalent positions
def _is_higher(best: list[int], contender: list[int], chain: Chain):
    """
    Checks if path is best fit to be the main one.
    
    Criteria:
    
    - Cycle (ring) for mixed chain is considered the main chain.

    """
    # TODO: order this to check on proper order
    # Should start counting from the closest to the functional group (trumping radicals)!
    # - Checks for empty values -
    if not contender:
        return False
    
    # - Checks for early exit -
    if not best:
        return True

    # - Check if it's an amine -
    # Kinda buklsome though...
    # for el_id in contender:
    #     el = chain.chain[el_id]
    #     if el == 'N':
    #         # Checks connection
    #         pass
    # - Check variables -
    contender_cyclical = False
    best_cyclical = False

    contender_hetero = False
    best_hetero = False

    contender_insat = 0
    best_insat = 0
    
    # For the most groups and the closest to start
    contender_group = []
    best_group = []

    # - Checks -
    if len(contender) > 1:
        if contender[0] == contender[-1]:
            contender_cyclical = True

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
        # Note that any group should flag this, not only substitutive groups
        if chain.to_els(curr_el.id, filter=lambda con: chain.get_main_path_id(con.to_id) > 0):
            contender_group.append(contender.index(id))

    if len(best) > 1:
        if best[0] == best[-1]:
            best_cyclical = True

    for id in best:
        curr_el = chain[id]
        # 1. Heteroatoms
        if curr_el.is_hetero():
            best_hetero = True
        # 2. Insaturation
        if any([con != SIMPLE for con in curr_el.cons if chain[con.to_id] == 'C']):
            best_insat += 1
        # 3. Closest group
        if chain.to_els(curr_el.id, filter=lambda con: chain.get_main_path_id(con.to_id) > 0):
            best_group.append(best.index(id))
    
    # - Comparisons -
    # 1. Cyclicality
    if contender_cyclical ^ best_cyclical:
        # XOR operator above
        return contender_cyclical

    # 2. Heteroatoms
    if contender_hetero ^ best_hetero:
        return contender_hetero
    
    # 3. Insaturation
    # Does nothing if they both are equal!
    if contender_insat > best_insat:
        return True
    elif contender_insat < best_insat:
        return False

    # 4. Group
    if contender_group and best_group:
        # 4.1. Sees which one has the most groups
        if len(contender_group) > len(best_group):
            # Contender has more groups 
            return True
        # 4.2. Sees nature of groups and follow alphabetical order
        # I think is better to numerate after....
        elif len(contender_group) == len(best_group):
            # Ordering pair function
            ordering_func = lambda el: (el[1], el[0])
            
            contender_mock_chain = chain
            # For now changes real object (not recomended!!!)
            # Though, appears to be a copy...?
            contender_mock_chain.set_main_path(contender)

            contender_subgs = iterate_subpaths(contender_mock_chain, contender_cyclical)

            # Separate elements into different types
            radical_types: dict = {}

            # Fills radical_types (based on size)
            for subg in contender_subgs:
                n_els = len(subg) - 1
                # Appends the position on the main chain
                m_pos = contender_mock_chain.get_main_path_id(subg[0])
                hetero_subg = has_hetero(subg)
                if n_els in radical_types:
                    radical_types[n_els].append((m_pos, hetero_subg))
                else:
                    radical_types[n_els] = [(m_pos, hetero_subg)]

            # Sorts radical_types ids in their names alphabetical order
            # Note that multiplier prefixes (di, tri, etc.) do not count
            # towards the alphabetical order (p. 70) 

            contender_func_names = list(map(lambda k: PREFIXES[k] if k < len(PREFIXES) else f"UNDEF_LEN[{k}]", radical_types.keys()))
            # Zips names with their values
            contender_pair_names = list(zip(contender_func_names, radical_types.values()))
            # Ordering:
            #   1. Based on if it's an heteroatom or not
            #   2. Based on first element (a.k.a its name)
            contender_pair_names.sort(key=ordering_func)

            best_mock_chain = chain
            # Techically is a copy!
            best_mock_chain.set_main_path(best)

            best_subgs = iterate_subpaths(best_mock_chain, best_cyclical)

            # Separate elements into different types
            radical_types: dict = {}

            # Fills radical_types (based on size)
            # TODO: Finish this logic!
            for subg in best_subgs:
                n_els = len(subg) - 1
                # Appends the position on the main chain
                m_pos = best_mock_chain.get_main_path_id(subg[0])
                hetero_subg = has_hetero(subg)
                if n_els in radical_types:
                    radical_types[n_els].append((m_pos, hetero_subg))
                else:
                    radical_types[n_els] = [(m_pos, hetero_subg)]
            # Sorts radical_types ids in their names alphabetical order
            best_func_names = list(map(lambda k: PREFIXES[k] if k < len(PREFIXES) else f"UNDEF_LEN[{k}]", radical_types.keys()))
            # Zips names with their values
            best_pair_names = list(zip(best_func_names, radical_types.values()))
            # Orders based on first element (a.k.a its name)
            best_pair_names.sort(key=ordering_func)
            # 4.2.1. Sees which one has the closest group to the start of the path
            #        (Or the one with most functional groups)
            contender_first_pair = contender_pair_names[0]
            best_first_pair = best_pair_names[0]
            # TODO: I need to know If a priorize longer radicals
            # 4.2.2. COmpares which one has the most heteroatoms
            # 4.2.3. Compares their main chain position for their "best" sorted elements
            if contender_first_pair[1][0][0] < best_first_pair[1][0][0]:
                return True
        else:
            # Current best has more groups
            False
    
    # 5. Size comparison
    if len(contender) > len(best):
        return True
    else:
        return False
    # If both are equal don't mtter which is returned
    # But return a False would be faster!

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

        # Checks to see if it's an ntriogen edge
        if curr_el == "N":
            best_path = curr_path
        else:
            # Appends initial direction on action stack
            for con in curr_el:
                action_stack.append((len(curr_path), con.to_id, con.dir))

        # Follows action_stack
        while action_stack:
            curr_start = 0
            # Goes to the next unless is an Nitrogen edge
            # If it's an nitrogen edge, it MUST be a amine
            # (That's the only way for an heteroatom to be an edge here)
            # A little out there code though
            curr_path_id, curr_el_id, origin_dir = action_stack.pop()
            curr_el = chain.chain[curr_el_id]
            print_field(
                chain.field,
                highlights=[chain.id_pool[curr_el_id]])
                
            # - End of Path check -
            # At edge
            # Detects if it's an edge (Primary) or at the start of path (Cyclical)
            # NOTE: Somewhat of a patchwork. Tbh I already know if is a cylce at this
            #       point, so maybe a less janky code could have been made here
            eop = (curr_el.classif == PRIM or curr_el.id == curr_path[0])
            
            if curr_path_id < len(curr_path):
                curr_path = curr_path[:curr_path_id]

            # Cyclical behaviour
            if curr_el_id in curr_path:
                eop = True
                # If so, cut out the cyclcal chain
                curr_start = curr_path.index(curr_el_id)
            
            # Such that heteroatoms are not on the main chain
            if curr_el.el not in HETEROATOMS:
                curr_path.append(curr_el_id)
            else:
                eop = True

            if eop:
                # Checks if curr_path is a better main path than best_path
                if _is_higher(best_path, curr_path[curr_start:], chain):
                    # New best_path
                    best_path = curr_path[curr_start:]
                print("Full path:")
                print_field(
                    chain.field,
                    highlights=[chain.id_pool[id] for id in curr_path[curr_start:]],
                    show_ids=True
                )
            else:
                for con in curr_el:
                    if con.dir != -origin_dir:
                        action_stack.append((len(curr_path), con.to_id, con.dir))
            eop = False


    return best_path


# TODO: Generalize as the main pathfinder?
def follow_subpath(chain: list[Entity], main_path: list[int], group: list[Entity], origin_dir: int, curr_id: int, cont: int = 0):
    """
    Follow subpath recursively (at each junction) until
    all entities outside of main_path are in a group
    
    :param chain: Chain with all entities
    :type chain: list[Entity]
    :param main_path: Main chain path (acts as mask to keep on group)
    :type main_path: list[int]
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
def run_chain(builder: tuple[tuple[str]] | str):
    """
    Go through the entirety of a chain field representation
    """
    chain = Chain(builder)
    
    chain.set_main_path(run_path_iterative(chain))
    
    return chain
