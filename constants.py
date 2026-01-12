import yaml
from pathlib import Path # Python 3.4+

# -- Class --
class FallbackDict:
    """
    Class of two "Dictionaries", one act as the primary dictionary
    while the second is the fallback in case the key isn't found on the
    priamry dictionary 
    """
    def __init__(self, primary_dict: dict, fallback_dict: dict):
        self.primary_dict: dict = primary_dict
        self.fallback_dict: dict = fallback_dict
    def __getitem__(self, key):
        if key in self.primary_dict:
            return self.primary_dict[key]
        elif key in self.fallback_dict:
            return self.fallback_dict[key]
        else:
            raise KeyError(f"{key} not in primary nor on default dict!")
    def __repr__(self):
        return self.primary_dict
    def __str__(self):
        return self.__repr__()

# -- Constants --
# - File related -
COMMENT_WILDCARD = '#'
FALLBACK_LANG = 'en'  # Needs to be english, as it's the only fully complete language
DEFAULT_LANG = 'pt-br'
EXTS = ['.field', '.smi']
CONF_PATH = Path('Conf')

# - String building related -
REVERSE_STR = False
CON_STR = ''
PLURAL_STR = 's'

# - Debug related -
PROLIX = False # Shows all functions print output if set to True
AUTOTEST = True

# Connection type
SIMPLE = 1
DOUBLE = 2
TRIPLE = 3
QUADRUPLE = 4

# Entity type (per Carbon connection)
# TODO: padronize so it is translatable
PRIM = 'primary'
SEC = 'secondary'
TERT = 'tertiary'
QUART = 'quaternary'

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
UNDEF = ''

SUFFIXES = {}
INFIXES = []
MULT_PREFS = []
AFFIXES = {}

HALIDES = {}
HETEROATOMS = {}

CLASSIFICATION = {}
SUBCLASSIFICATION = {}

CCLASSIFICATION = {}
# Failsafe for iterative routines
MAX_ITER = 21

def load_constants(lang: str = ''):
    # -- Configurations file --
    # - Conf constants -
    global DEFAULT_LANG
    global PROLIX
    global AUTOTEST
    # First check if there is such file
    conf_file: Path = CONF_PATH / 'conf.yaml'

    if not conf_file.is_file():
        raise FileNotFoundError(f'Configuration file {conf_file} missing!')

    # - Get its content -
    with open(conf_file, 'r') as f:
        doc = yaml.load(f, Loader=yaml.FullLoader)

        DEFAULT_LANG = doc['default_language']
        PROLIX = doc['prolix_debug']
        AUTOTEST = doc['autotest']

    if not lang:
        lang = DEFAULT_LANG

    # - Language resource files -
    file_path: Path = CONF_PATH / (lang + '.yaml')
    fallback_file_path: Path = CONF_PATH / (FALLBACK_LANG + '.yaml')

    # - Checks if language is available -
    # This triggers one error at a time if both the default and current laguage
    # Aren't set correctly
    if not file_path.is_file():
        raise FileNotFoundError(f"Language '{lang}' not supported! (File '{file_path}' not found!)")

    if not fallback_file_path.is_file():
        raise FileNotFoundError(f"Fallback language '{FALLBACK_LANG}' not found! (File '{fallback_file_path}' not found!)")

    # In case there isn't the info on that language defaults to the english one
    # Constants to be changed
    global REVERSE_STR
    global CON_STR
    global UNRESOLVED
    global UNDEF
    global SUFFIXES
    global INFIXES
    global MULT_PREFS
    global AFFIXES
    global HALIDES
    global HETEROATOMS
    global CLASSIFICATION
    global SUBCLASSIFICATION
    global CCLASSIFICATION

    # - Read YAML resource file -
    # 1. Fallback language
    fallback_res: dict = None
    with open(fallback_file_path, 'r') as f:
        doc = yaml.load(f, Loader=yaml.FullLoader)

        if FALLBACK_LANG == doc['lang']:
            fallback_res = doc
    if not fallback_res:
        raise ValueError(f"File is empty or invalid! ({fallback_file_path})")
    
    # 2. Chosen language (only if they are different)
    if FALLBACK_LANG != lang:
        res: FallbackDict = None
        with open(file_path, 'r') as f:
            doc = yaml.load(f, Loader=yaml.FullLoader)

            if lang == doc['lang']:
                res = FallbackDict(primary_dict=doc, fallback_dict=fallback_res)
        if not res:
            raise ValueError(f"File is empty or invalid! ({file_path})")
    else:
        res = fallback_res
    
    REVERSE_STR = res['reverse_class']
    CON_STR = res['conective_str']

    UNRESOLVED = res['unresolved_str']
    UNRESOLVED = res['undefined_str']

    # Merge of dictionaries ( Python >= 3.9.0 )
    SUFFIXES = res['suffixes'] | {UNRESOLVED:'NONE'}
    INFIXES = res['infixes']
    MULT_PREFS = res['multiplier_prefixes']
    AFFIXES = res['afixes']

    HALIDES = res['heteroatoms']['halides']
    HETEROATOMS = {}
    for htype in res['heteroatoms']:
        HETEROATOMS |= res['heteroatoms'][htype]

    CLASSIFICATION = res['classifcation']
    # SUBCLASSIFICATION = res['subclassification']
    CCLASSIFICATION = res['carbon_classification']

# For now, loads it here
load_constants()