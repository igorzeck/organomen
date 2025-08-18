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
# TODO: Make field tuple?
# TODO: Make it detect edges that are only on Carbon!
# TODO: Shy away from field!!! Using chain representation instead
# TODO: Make it so that main.chain is ordered from carbon 1 to n (in naming order)
# -- Imports --
from pathfinder import run_chain
from classification import class_chain
from entities import *


def main():
    chain_path = "Chains/simple.field"
    # Pathfinding
    # print(Entity(1, 'C', [Connection(1, 3, 1, 1), Connection(1, 2, 2, 1)]))
    c_main = Chain(chain_path)
    print(c_main.chain)
    # for _c in c_main.chain:
    #     print(_c)
    c_main = run_chain(c_main.field)
    class_chain(c_main)
    print(c_main)


# -- Start --
if __name__ == "__main__":
    main()