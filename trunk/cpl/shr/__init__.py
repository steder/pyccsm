"""
pyCPL6 Package

This package contains overloaded MCT types to correspond to the 
CPL6 types (e.g.:  Bundle, Domain, Contract), as well as an F2Py wrapped 
MPH (Multi Program/Component Handshaking) module.

This package contains the following modules:
  * comm
  * fields
"""

def _get_exports_list(module):
    try:
        return list(module.__all__)
    except AttributeError:
        return [n for n in dir(module) if n[0] != '_']

__all__ = []

import const

__all__.append( "const" )

del _get_exports_list
