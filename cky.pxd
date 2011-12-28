cimport numpy as np
from containers cimport ChartItem, Edge, Rule, Terminal, Grammar
#from array cimport array

cdef extern from "math.h":
	bint isinf(double x)
	bint isnan(double x)
	bint isfinite(double x)
	double log(double x)
	double exp(double x)

cdef inline ChartItem new_ChartItem(unsigned int label, short left, short right)
cdef inline Edge new_Edge(double inside, Rule rule, short split)
