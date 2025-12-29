# Test file for the chain resolution algorithm
# -- Definition --
# Stubs: parts of a path
# Path: is a collection of carbon starting
# and ending with a one connection carbon
# or with the same 2 carbons on the path
# Chain: made of all unique combinations
# of paths, which in turn define how
# the chain is named

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

def main():
    chain_path = "Local/Chains/triclorinefluorinemethane.field"

    # Pathfinding
    c_main = run_chain(chain_path)
    class_chain(c_main)
    print(c_main)

    # c_test = Chain()
    # c_test.load_chain(c_main.chain)
    # print_field(c_test.field)
    # c_test = run_chain(c_test.field)
    # class_chain(c_test)
    # print(c_test)



# -- Start --
if __name__ == "__main__":
    main()