cimport cython
cimport numpy as np
from kbest cimport lazykthbest, lazykbest
from agenda cimport Entry
from containers cimport ChartItem, Edge, RankedEdge, Grammar, \
			dictcast, getlabel, getvec, itemcast, edgecast
from cky cimport new_ChartItem

@cython.locals(
	kbest=dict,
	d=dict,
	prunetoid=dict,
	prunetolabel=dict,
	toid=dict,
	tolabel=dict,
	prunelist=list,
	Ih=ChartItem)
cpdef whitelistfromchart(dict chart,
		ChartItem goal, Grammar coarse, Grammar fine, int k,
			np.ndarray[np.double_t, ndim=3] whitelist, int maxlen)

@cython.locals(
	entry=Entry,
	e=Edge,
	D=dict,
	outside=dict)
cdef dict kbest_outside(dict chart, ChartItem start, int k)

@cython.locals(
	e=Edge,
	eleft=ChartItem,
	eright=ChartItem,
	eejj=RankedEdge,
	prob=cython.double,
	entry=Entry)
cdef void getitems(RankedEdge ej, double rootprob, dict D,
							dict chart, dict outside)

cpdef filterchart(chart, start)

@cython.locals(
	edge=Edge,
	item=ChartItem)
cdef void filter_subtree(ChartItem start, dict chart, dict chart2)
