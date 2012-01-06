cimport cython
cimport numpy as np
from agenda cimport Entry, Agenda, nsmallest
from containers cimport ChartItem, Edge, RankedEdge, new_RankedEdge

@cython.locals(edge=Edge)
cdef inline getcandidates(list chart, ChartItem v, int k)

@cython.locals(ej=RankedEdge, entry=Entry)
cpdef inline lazykthbest(ChartItem v, int k, int k1,
			dict D, dict cand, list chart, set explored)

@cython.locals(ej1=RankedEdge, prob=cython.double)
cdef inline lazynext(RankedEdge ej, int k1,
			dict D, dict cand, list chart, set explored)

@cython.locals(result=cython.double, prob=cython.double, i=cython.int,
			eleft=ChartItem, eright=ChartItem, edge=Edge, entry=Entry)
cdef inline double getprob(list chart, dict D, RankedEdge ej) except -1

@cython.locals(ej=RankedEdge, edge=Edge, children=list, ei=ChartItem,
	i=cython.int, entry=Entry)
cdef inline str getderivation(RankedEdge ej,
			dict D, list chart, dict tolabel, list sent, int n)

@cython.locals(entry=Entry, prob=cython.double)
cpdef list lazykbest(list chart, ChartItem goal, int k, dict tolabel, list sent)
