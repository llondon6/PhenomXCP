#!/usr/bin/env python
'''
init for python development of PhenomXCP
'''

# Import usefuls
from positive import *
alert(yellow('Warm greetings from XCP')+'.',fname='init')

# Import local dependencies
from . import core
from . core import *
from . import io
from . io import *

from os.path import exists
if exists("/Users/book/KOALA/PhenomXCP/xcp/parameter_space_fits.py"):
    # Import the python version of the parameter space fits if it exists
    from .parameter_space_fits import *