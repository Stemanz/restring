# TODO: installation via Pip needs to be fixed
__version__ = "0.1.8"

# the following tests if I can import stuff from restring.py
#from .restring import files # this works

# the following tests if I can import stuff from restring.py
from restring.restring import working_directory, files #this also works

# the following tests if I can import stuff from gears.py, and where it will end up
from .gears import manzlog # this works

# the following tests Giorgio's advice about importing
from restring.gears import get_dirs

# the following tests if I can import stuff from gears.py, and where it will end up
#from gears.gears import get_dirs # this fails:
#Traceback (most recent call last):
#  File "<stdin>", line 1, in <module>
#  File "/Users/manz/opt/anaconda3/envs/py38/lib/python3.8/site-packages/restring/__init__.py", line 14, in <module>
#    from gears.gears import get_dirs
#ModuleNotFoundError: No module named 'gears'

__all__ = ("files", "working_directory", "manzlog", "get_dirs")
