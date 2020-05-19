# Licensed under MIT license - see LICENSE.rst

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from .._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from stingray.simulator.base import *
    from stingray.simulator.simulator import *
    from stingray.simulator.transfer import *
    from stingray.simulator.models import *
     
