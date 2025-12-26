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
# TODO: make unique hold only unique paths
# TODO: Organize in different files and make the integration with Chain class better
# TODO: Naming arguments in a yaml or simmilar to make it possible translations
# TODO: Shy away from field!!! Using chain representation instead
# TODO: Make it so that main.chain is ordered from carbon 1 to n (in naming order)
# TODO: Add config file
# TODO: Support for functional naming (e.g. for alcohols)
# TODO: Make rules for hyphens better!
# TODO: Ver sobre integrar com formato de arquivos como o do https://pubchem.ncbi.nlm.nih.gov/compound/6-Ethyl-3-methylnonane#section=DSSTox-Substance-ID
# -- Imports --
from pathfinder import run_chain
from classification import class_chain
from entities import *


def main():
    chain_path = "Chains/1-ethyl-3-methylheptane.field"

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