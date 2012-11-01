""" Implementation of Huang & Chiang (2005): Better k-best parsing. """
import logging
from math import exp
from agenda import Agenda
from containers import Edge, RankedEdge
from operator import itemgetter

from agenda cimport Entry, Agenda, nsmallest
from containers cimport Edge, RankedEdge

cdef tuple unarybest = (0, ), binarybest = (0, 0)

cpdef list lazykbest(list chart, unsigned int start, short left, short right,
		int k, dict tolabel):
	""" wrapper function to run lazykthbest and get derivations. """
	cdef Entry entry
	cdef double prob
	cdef list D = [[{} for _ in x] for x in chart]
	cdef list cand = [[{} for _ in x] for x in chart]
	cdef set explored = set()
	lazykthbest(start, left, right, k, k, D, cand, chart, explored)
	logging.debug("%d (sub)derivations considered", len(explored))
	return filter(itemgetter(0), [
		(getderivation(entry.key, D, chart, tolabel, 0, None), entry.value)
		for entry in D[left][right][start]])

cdef inline getcandidates(list chart, unsigned int label,
		short start, short end, int k):
	""" Return a heap with up to k candidate arcs starting from vertex v """
	# NB: the priority queue should either do a stable sort, or should
	# sort on rank vector as well to have ties resolved in FIFO order;
	# otherwise the sequence (0, 0) -> (1, 0) -> (1, 1) -> (0, 1) -> (1, 1)
	# can occur (given that the first two have probability x and the latter
	# three probability y), in which case insertion order should count.
	# Otherwise (1, 1) ends up in D[v] after which (0. 1) generates it
	# as a neighbor and puts it in cand[v] for a second time.
	cdef Edge edge
	if label not in chart[start][end]: return Agenda()
	cell = chart[start][end][label]
	return Agenda(
		[(RankedEdge(label, start, end, edge, 0, 0 if edge.rule is not None
						and edge.rule.rhs2 else -1), edge.inside)
						for edge in nsmallest(k, cell)])

cpdef inline lazykthbest(unsigned int label, short start, short end,
		int k, int k1, list D, list cand, list chart, set explored):
	cdef Entry entry
	cdef RankedEdge ej
	# k1 is the global k
	# first visit of vertex v?
	try:
		if label not in cand[start][end]:
			# initialize the heap
			cand[start][end][label] = getcandidates(chart, label, start, end, k1)
	except IndexError:
		print label, start, end
		print cand
		print len(cand)
		print len(cand[0])
		print len(cand[0][1])
		raise
	while label not in D[start][end] or len(D[start][end][label]) < k:
		if label in D[start][end]:
			# last derivation
			entry = D[start][end][label][-1]
			ej = entry.key
			# update the heap, adding the successors of last derivation
			lazynext(ej, k1, D, cand, chart, explored)
		# get the next best derivation and delete it from the heap
		if cand[start][end][label]:
			D[start][end].setdefault(label, []).append(
				cand[start][end][label].popentry())
		else: break
	return D

cdef inline lazynext(RankedEdge ej, int k1, list D, list cand, list chart,
		set explored):
	cdef RankedEdge ej1
	cdef Edge ec = ej.edge
	cdef double prob
	cdef unsigned int label
	cdef short start, end
	# add the |e| neighbors
	# left child
	label = ec.rule.rhs1; start = ej.left; end = ec.split
	ej1 = RankedEdge(ej.label, ej.left, ej.right, ej.edge,
			ej.leftrank + 1, ej.rightrank)
	# recursively solve a subproblem
	# NB: increment j1[i] again because j is zero-based and k is not
	lazykthbest(label, start, end, ej1.leftrank + 1, k1, D, cand, chart, explored)
	# if it exists and is not in heap yet
	if ((label in D[start][end] and ej1.leftrank < len(D[start][end][label]))
		and ej1 not in explored): #cand[ej1.head]): <= gives duplicates
		prob = getprob(chart, D, ej1)
		# add it to the heap
		cand[ej1.left][ej1.right][ej1.label][ej1] = prob
		explored.add(ej1)
	# right child?
	if ej.rightrank == -1: return
	label = ec.rule.rhs2; start = ec.split; end = ej.right
	ej1 = RankedEdge(ej.label, ej.left, ej.right, ej.edge,
			ej.leftrank, ej.rightrank + 1)
	lazykthbest(label, start, end, ej1.rightrank + 1, k1,
						D, cand, chart, explored)
	# if it exists and is not in heap yet
	if ((label in D[start][end] and ej1.rightrank < len(D[start][end][label]))
		and ej1 not in explored): #cand[ej1.head]): <= gives duplicates
		prob = getprob(chart, D, ej1)
		# add it to the heap
		cand[ej1.left][ej1.right][ej1.label][ej1] = prob
		explored.add(ej1)

cdef inline double getprob(list chart, list D, RankedEdge ej) except -1.0:
	cdef Edge ec, edge
	cdef Entry entry
	cdef double result, prob
	ec = ej.edge
	label = ec.rule.rhs1; start = ej.left; end = ec.split
	if label in D[start][end]:
		entry = D[start][end][label][ej.leftrank]
		prob = entry.value
	elif ej.leftrank == 0: edge = chart[start][end][label][0]; prob = edge.inside
	else: raise ValueError(
		"non-zero rank vector not part of explored derivations")
	# NB: edge.inside if preterminal, 0.0 for terminal
	result = ec.rule.prob + prob
	if ej.rightrank >= 0: #if e.right.label:
		label = ec.rule.rhs2
		start = ec.split; end = ej.right
		if label in D[start][end]:
			entry = D[start][end][label][ej.rightrank]
			prob = entry.value
		elif ej.rightrank == 0:
			edge = chart[start][end][label][0]
			prob = edge.inside
		else: raise ValueError(
			"non-zero rank vector not part of explored derivations")
		result += prob
	return result

cdef inline getderivation(RankedEdge ej, list D, list chart,
		dict tolabel, int n, str debin):
	""" Translate the (e, j) notation to an actual tree string in
	bracket notation.  e is an edge, j is a vector prescribing the rank of the
	corresponding tail node. For example, given the edge <S, [NP, VP], 1.0> and
	vector [2, 1], this points to the derivation headed by S and having the 2nd
	best NP and the 1st best VP as children.
	If `debin' is specified, will perform on-the-fly debinarization of nodes
	with labels containing `debin' an a substring. """
	cdef RankedEdge rankededge
	cdef Entry entry
	cdef Edge edge
	cdef str child, children = ""
	cdef int i = ej.leftrank
	cdef unsigned int label
	cdef short start, end
	if n > 100: return ""	#hardcoded limit to prevent cycles
	label = ej.edge.rule.rhs1
	start = ej.left; end = ej.edge.split
	while i != -1:
		if label not in chart[start][end]:
			# this must be a terminal
			children = " %s" % ej.edge.rule.word.encode('utf-8')
			break
		if label in D[start][end]:
			rankededge = (<Entry>D[start][end][label][i]).key
		else:
			assert i == 0, "non-best edge missing in derivations"
			entry = getcandidates(chart, label, start, end, 1).popentry()
			D[start][end][label] = [entry]
			rankededge = entry.key
		child = getderivation(rankededge, D, chart, tolabel, n + 1, debin)
		if child == "": return ""
		children += " %s" % child
		if end == ej.right: break
		label = ej.edge.rule.rhs2
		start, end = end, ej.right
		i = ej.rightrank
	if debin is not None and debin in tolabel[ej.label]:
		return children
	return "(%s%s)" % (tolabel[ej.label], children)

def getderiv(RankedEdge ej, dict D, dict chart, dict tolabel, int n, str debin):
	return getderivation(ej, D, chart, tolabel, n, debin)
