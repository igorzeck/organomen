# Test file for the chain resolution algorithm
# -- Definition --
# Stubs: parts of a path
# Path: is a collection of carbon starting
# and ending with a one connection carbon
# or with the same 2 carbons on the path
# Chain: made of all unique combinations
# of paths, which in turn define how
# the chain is named

# TODO: the current way I'm flaggin main chain don't support single "bridges" between cycles!
# TODO: Function that makes ordering of electrons orbital using octate rule so that I can have connection numbers per element!
# TODO: Flag type of amine and type of alcohol
# TODO: Add configuration file with element info, including maximum conenctions!
# TODO: Refactor (or add) docstring for functions!!!
# TODO: Organize in different files and make the integration with Chain class better
# TODO: Shy away from field!!! Using chain representation instead
# TODO: Add config file
# TODO: Support for functional naming (e.g. for alcohols)
# TODO: Make rules for hyphens better!
# TODO: Integate this file extension: https://pubchem.ncbi.nlm.nih.gov/compound/6-Ethyl-3-methylnonane#section=DSSTox-Substance-ID
# -- Imports --
# from base_structures import load_constants
# load_constants(lang = 'pt-br')  # Load constants based on language

from pathfinder import run_chain
from classification import class_chain
from entities import *
from constants import EXTS

import sys

def main():
    # File handling
    recent_files_path: Path = Path('Res') / 'recent_files'
    default_path = ''

    if recent_files_path.is_file():
        with open(recent_files_path, "r") as arq:
            default_path = arq.readline()
    else:
        print("No recent files found!")
    
    curr_f = sys.argv[1]

    chain_path = curr_f if curr_f.rpartition('.')[-1].lower() in EXTS else default_path

    # Pathfinding
    c_main = run_chain(chain_path)
    class_chain(c_main)
    print(c_main)


# -- Start --
if __name__ == "__main__":
    main()