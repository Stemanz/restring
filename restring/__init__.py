# TODO: installation via Pip needs to be fixed
__version__ = "0.1.6"

# the following tests if I can import stuff from restring.py
from .restring import files

# the following tests if I can import stuff from restring.py
from restring.restring import working_directory #this also works

# the following tests if I can import stuff from gears.py, and where it will end up
from .gears import manzlog, get_dirs

# the following tests if I can import stuff from gears.py, and where it will end up
#from gears.gears import get_dirs # this fails:
#Traceback (most recent call last):
#  File "<stdin>", line 1, in <module>
#  File "/Users/manz/opt/anaconda3/envs/py38/lib/python3.8/site-packages/restring/__init__.py", line 14, in <module>
#    from gears.gears import get_dirs
#ModuleNotFoundError: No module named 'gears'

__all__ = ("files", "working_directory", "manzlog", "get_dirs")
