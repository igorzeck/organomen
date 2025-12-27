import yaml

# -- Constants --
# - File related -
COMMENT_WILDCARD = '#'
DEFAULT_LANG = 'pt-br'
EXTS = ['field']

REVERSE_STR = False
CON_STR = ''
PLURAL_STR = 's'
# Enun would be just lovely, but alas...
# Connections directions
NORTHWEST = -4
NORTH = -3
NORTHEAST = -2
WEST = -1
EAST = 1
SOUTHWEST = 2
SOUTH = 3
SOUTHEAST = 4

# Connection type
SIMPLE = 1
DOUBLE = 2
TRIPLE = 3
QUADRUPLE = 4

# Relative to N of Cs in main chain
PREFIXES = [
    '',
    'met',
    'et',
    'prop',
    'but',
    'pent',
    'hex',
    'hept',
    'oct',
    'non',
    'dec',
    'undec',
    'dodec',
    'tridec',
    'tetradec',
    'pentadec',
    'hexadec',
    'heptadec',
    'octadec',
    'nonadec',
    'icos'
]

# Relative to chain functional class
# Maybe use a default dict?
# NOTE: Atoms are always in all caps for ease of comparison
UNRESOLVED = ''

SUFFIXES = {}
INFIXES = []
AFIXES = {}

HALIDES = {}
HETEROATOMS = {}

CLASSIFICATION = {}
# Failsafe for iterative routines
MAX_ITER = 21

def load_constants(lang: str = ""):
    # TODO: get default language from a separete configuration file
    if not lang:
        lang = DEFAULT_LANG
    # TODO: separete into single files
    # In case there isn't the info on that language defaults to the english one
    # Constants to be changed
    global REVERSE_STR
    global CON_STR
    global SUFFIXES
    global INFIXES
    global AFIXES
    global HALIDES
    global HETEROATOMS
    global CLASSIFICATION

    # Read YAML resource file
    res = None
    with open("Res/conf.yaml", 'r') as f:
        data = yaml.load_all(f, Loader=yaml.FullLoader)

        # Somewhat cubersome to access w/ multiple docs, but let it be this way
        for doc in data:
            if lang == doc['lang']:
                res = doc
    if not res:
        raise ValueError("Not a valid language!")
    
    REVERSE_STR = res['reverse_class']
    CON_STR = res['conective_str']

    UNRESOLVED = res['unresolved_str']
    # Merges dictionaries ( Python >= 3.9.0 )
    SUFFIXES = res['suffixes'] | {UNRESOLVED:'NONE'}
    INFIXES = res['infixes']
    AFIXES = res['afixes']

    HALIDES = res['heteroatoms']['halides']
    HETEROATOMS = {}
    for htype in res['heteroatoms']:
        HETEROATOMS |= res['heteroatoms'][htype]

    CLASSIFICATION = res['classifcation']

# By defautl run
load_constants()