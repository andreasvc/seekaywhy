from distutils.core import setup
from Cython.Build import cythonize
import numpy

# some of these directives increase performance,
# but at the cost of failing in mysterious ways.
directives = {
	"profile" : False,
	"nonecheck" : False,
	"cdivision" : True,
	"wraparound" : True,
	"boundscheck" : False,
	"embedsignature" : True,
	#"extra_compile_args" : ["-O3"],
	#extra_link_args=["-g"]
	}

setup(
	name = 'seekaywhy',
	include_dirs = [numpy.get_include()],
	ext_modules = cythonize(['*.pyx', '*.py'],
		exclude=['setup.py', 'parser.py'],
		pyrex_directives=directives,
		directives=directives,
		**directives)
)
