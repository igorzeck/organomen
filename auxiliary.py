from base_structures import Pos, CD_OFFS
# Pathfinding and scout base functions

# -- Functions --
# TODO: Put this function inside scout?
def _is_inside(pos: Pos, n_row:int, n_col:int):
    if (pos.row >= n_row or pos.col >= n_col)\
        or\
        (pos.row < 0 or pos.col < 0):
        return False
    return True


# TODO: Version to accept positions AND id positions
# Scout via id_pool and chain itself perhaps?
# Those params are kinda cluttered!
def scout(field, pos_pool: tuple[Pos], pos_id: int, nxt_id = True, nxt_dir = False, nxt_con_val = False, nxt_str = False):
    """
    Scouts the vicinity of an given position
    field: Chain field
    ids: Chain tuple containing a table of ids to postion in the field
    pos_id: Pos id to investigate around
    direction: If True returns direction
    con_value: If true returns connection value
    returns: list of a tuple (package) with up to the following elements (nxt_id, dir, type)
    """
    nxt_el_l = []
    for dir in CD_OFFS:
        pos = pos_pool[pos_id]
        # Next direction and connectin (to scout)
        nxt_step = CD_OFFS[dir]
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
                selected_id = None
                nxt_pos = None
                for i_pos, _pos in enumerate(pos_pool):
                    if _pos == nxt_el_pos:
                        selected_id = i_pos
                        nxt_pos = _pos
                if selected_id is None:
                    raise ValueError(f"Didn't find position {nxt_el_pos} id's")

                _package = []
                if nxt_id:
                    _package.append(selected_id)
                if nxt_dir:
                    _package.append(dir)
                if nxt_con_val:
                    _package.append(next_con)
                if nxt_str and nxt_pos:
                    _package.append(field[nxt_pos.row][nxt_pos.col])
                nxt_el_l.append(_package)
    return nxt_el_l