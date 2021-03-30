# TODO: installation via Pip needs to be fixed
__version__ = "0.1.12"

# the following tests if and how I can import stuff
#from .restring import files # this works
#from restring.restring import working_directory #this also works
#from .gears import manzlog # this works
#from restring.gears import get_dirs # this also works

from restring.restring import restring_gui

__all__ = ("restring_gui")
