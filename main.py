# Test file for the chain separation algorithm
# -- Definition --
# Stubs are parts of a path

# A path is a collection of carbon starting
# and ending with a one connection carbon
# or with the same 2 carbons on the path

# A chain is made of all unique combinations
# of paths, which in turn define how
# the chain is named

# TODO: Make part division for subnaming
# TODO: make unique hold only unique paths
# TODO: Chain object can be classified as it is read
# TODO: Organize in different files and make the integration with Chain class better
# TODO: Naming arguments in a yaml or simmilar to make it possible translations
# -- Imports --
from pathfinding import per_chain
from classification import name_chain
from entities import *


def main():
    chain_path = "Chains/simple.chain"
    # Pathfinding
    c_main = Chain(chain_path)
    c_main = per_chain(c_main.field)
    c_main.name = name_chain(c_main.field, c_main.id_dict, c_main.paths, c_main.main_path)
    print(c_main)


# -- Start --
if __name__ == "__main__":
    main()