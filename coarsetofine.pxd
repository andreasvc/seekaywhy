cimport cython
cimport numpy as np
from kbest cimport lazykthbest, lazykbest
from agenda cimport Entry, Agenda
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
cpdef whitelistfromkbest(dict chart,
		ChartItem goal, Grammar coarse, Grammar fine, int k,
			np.ndarray[np.double_t, ndim=3] whitelist, int maxlen)

@cython.locals(
	label=cython.long,
	lensent=cython.short)
cpdef whitelistfromposteriors(np.ndarray[np.double_t, ndim=3] inside,
	np.ndarray[np.double_t, ndim=3] outside,
	ChartItem goal,
	Grammar coarse,
	Grammar fine,
	np.ndarray[np.double_t, ndim=3] whitelist,
	short maxlen,
	double threshold)

cpdef whitelistfromposteriors1(
	dict chart,
	np.ndarray[np.double_t, ndim=3] viterbi,
	ChartItem goal,
	Grammar coarse,
	Grammar fine,
	np.ndarray[np.double_t, ndim=3] whitelist,
	double threshold)

@cython.locals(
	label=cython.long,
	lensent=cython.short)
cpdef whitelistfromposteriors2(np.ndarray[np.double_t, ndim=3] inside,
	np.ndarray[np.double_t, ndim=3] outside,
	ChartItem goal,
	Grammar coarse,
	Grammar fine,
	dict mapping,
	short maxlen,
	double threshold)


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

@cython.locals(
	item=ChartItem,
	edges=list)
cpdef insidescores(dict chart, ChartItem goal, np.ndarray[np.double_t, ndim=3] inside)

@cython.locals(
	agenda=Agenda,
	visited=set,
	item=ChartItem,
	newitem=ChartItem,
	edge=Edge)
cpdef outsidescores(dict chart, ChartItem goal, np.ndarray[np.double_t, ndim=3] inside, np.ndarray[np.double_t, ndim=3] outside)
