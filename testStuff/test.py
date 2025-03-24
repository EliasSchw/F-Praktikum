import numpy as np


import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))
from macroswriter import writeLatexMacro # Import macroswriter module


writeLatexMacro('testOhneFehler', 123456789, unit='\\frac{1}{m}')
writeLatexMacro('testMitFehler', 123456789, unit='\\frac{1}{m}', error=20000000)

