import entities as ent
import os.path
#
#  File opening
#
main_chain = ent.open_chain_file(os.path.join("Chains", "simple.chain"))
print(main_chain)
