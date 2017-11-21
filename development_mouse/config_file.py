import os
import logging


DEFAULT_CONFIG: str = """
# This is an automatically generated config file for cytograph
"""


def create_rc() -> None:
    """Creates default .cytographrc file in the user home folder
    """
    try:
        path = os.path.join(os.path.expanduser("~"), ".cytographrc")
        with open(path, "w") as source:
            pass
    except IOError:
        logging.error(f"Writing to {path} not allowed")


def load_rc() -> str:
    """Load the .cytographrc file, trying for possible loacations

    Adapted from https://stackoverflow.com/questions/7567642/where-to-put-a-configuration-file-in-python
    """
    config_txt: str = None
    for loc in os.curdir, os.path.expanduser("~"), os.environ.get("CYTHOGRAPH_CONF"):
        try:
            with open(os.path.join(loc, ".cytographrc")) as source:
                config_txt = source.read()
            break
        except IOError:
            logging.warn(".cytographrc not found")
    
    return config_txt
