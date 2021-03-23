# TODO: installation via Pip needs to be fixed
__version__ = "0.1.1"

# the following tests if I can import stuff from restring.py
from .restring import files

# the following tests if I can import stuff from restring.py
from restring.restring import working_directory

# the following tests if I can import stuff from gears.py, and where it will end up
from .gears import manzlog

# the following tests if I can import stuff from gears.py, and where it will end up
from gears.gears import get_dirs

__all__ = ("files", "working_directory", "manzlog", "get_dirs")
