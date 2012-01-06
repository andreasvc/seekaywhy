cimport numpy as np
from containers cimport ChartItem, Edge, Rule, Terminal, Grammar

cdef extern from "math.h":
	bint isinf(double x)
	bint isnan(double x)
	bint isfinite(double x)
	double log(double x)
	double exp(double x)
