from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

# some of these directives increase performance,
# but at the cost of failing in mysterious ways.
directives = {
	"profile" : False,
	"boundscheck" : False,
	"nonecheck" : False,
	"wraparound" : False,
	"embedsignature" : True
	}

ext_modules = [
	Extension("cky",			["cky.pyx",  "cky.pxd"],
								#extra_compile_args=["-g"], extra_link_args=["-g"]
								),
	Extension("kbest",			["kbest.py", "kbest.pxd"],
								#extra_compile_args=["-g"], extra_link_args=["-g"]
								),
	Extension("agenda",			["agenda.pyx", "agenda.pxd"],
								#extra_compile_args=["-g"], extra_link_args=["-g"]
								),
	Extension("containers",		["containers.pyx", "containers.pxd"],
								#extra_compile_args=["-g"], extra_link_args=["-g"]
								),
	Extension("coarsetofine",	["coarsetofine.py", "coarsetofine.pxd"],
								#extra_compile_args=["-g"], extra_link_args=["-g"]
								),
	Extension("disambiguation",	["disambiguation.py", "disambiguation.pxd"]),
	]

for e in ext_modules:
	e.pyrex_directives = directives

setup(
	name = 'seekaywhy',
	cmdclass = {'build_ext': build_ext},
	include_dirs = [np.get_include(), '.'],
	ext_modules = ext_modules
)
