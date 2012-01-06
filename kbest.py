""" Implementation of Huang & Chiang (2005): Better k-best parsing
"""
from math import exp, fsum
import numpy as np
from agenda import Agenda, Entry
from containers import ChartItem, Edge, RankedEdge
from operator import itemgetter
try: assert nsmallest(1, [1]) == [1]
except NameError: from heapq import *

unarybest = (0, )
binarybest = (0, 0)

#def getkbestparses(k, sent, goal, viterbi, grammar):
#	start = 0
#	end = len(sent)
#	kbesttrees = []
#	explored = set()
#	D = {}
#	cand = {}
#	for i in range(k + 1):
#		tree = gettree(goal, i, k, D, grammar)
#		if tree is None: break
#		kbesttrees.append((tree, D[goal][-1]))
#	return kbesttrees
#
#def gettree(v, k, k1, D, grammar):
#	lazykthbest(v, k, k1, D, cand, viterbi, explored)
#	if v.rhs1 == 0: #is tag
#		# ...
#		return tagnode
#
#	if k - 1 >= len(D[v]): return None
#
#	d = D[v][0]
#	child = ChartItem(d.edge.rule.rhs1, v.left, d.edge.split)
#	t = gettree(child, d.left, k1, D, grammar)
#	children = [t]
#	if d.right != -1:
#		child = ChartItem(d.edge.rule.rhs2, d.edge.split, v.right)
#		t = gettree(child, d.right, k1, D, grammar)
#		children.append(t)
#	return "(%s %s)" % (grammar.tolabel[v.label], "".join(children))

def getcandidates(chart, v, k):
	""" Return a heap with up to k candidate arcs starting from vertex v """
	# NB: the priority queue should either do a stable sort, or should
	# sort on rank vector as well to have ties resolved in FIFO order;
	# otherwise the sequence (0, 0) -> (1, 0) -> (1, 1) -> (0, 1) -> (1, 1)
	# can occur (given that the first two have probability x and the latter
	# three probability y), in which case insertion order should count.
	# Otherwise (1, 1) ends up in D[v] after which (0. 1) generates it
	# as a neighbor and puts it in cand[v] for a second time.
	if not chart[v.left][v.right].get(v.label, False):
		#shouldn't raise error because terminals should end here
		#raise ValueError("%r not in chart" % v)
		return Agenda()
	return Agenda(
		[(RankedEdge(v, edge, 0, 0 if edge.rule.rhs2 else -1), edge.inside)
					for edge in nsmallest(k, chart[v.left][v.right][v.label])])

def lazykthbest(v, k, k1, D, cand, chart, explored):
	# k1 is the global k
	## first visit of vertex v?
	# initialize the heap
	if v not in cand:
		cand[v] = getcandidates(chart, v, k1)
	if v not in D: D[v] = []
	while len(D[v]) < k:
		if D[v]:
			# last derivation
			ej = D[v][-1].key
			# update the heap, adding the successors of last derivation
			lazynext(ej, k1, D, cand, chart, explored)
		# get the next best derivation and delete it from the heap
		if cand[v]:
			D[v].append(cand[v].popentry())
		else: break
	return D

def lazynext(ej, k1, D, cand, chart, explored):
	# add the |e| neighbors
	for i in range(2):
		if i == 0:
			ei = ChartItem(ej.edge.rule.rhs1, ej.head.left, ej.edge.split)
			ej1 = RankedEdge(ej.head, ej.edge, ej.left + 1, ej.right)
		elif i == 1 and ej.right >= 0: #edge.right.label:
			ei = ChartItem(ej.edge.rule.rhs2, ej.edge.split, ej.head.right)
			ej1 = RankedEdge(ej.head, ej.edge, ej.left, ej.right + 1)
		else: break
		# recursively solve a subproblem
		# NB: increment j1[i] again because j is zero-based and k is not
		lazykthbest(ei, (ej1.right if i else ej1.left) + 1, k1,
							D, cand, chart, explored)
		# if it exists and is not in heap yet
		if ((ei in D and (ej1.right if i else ej1.left) < len(D[ei]))
			and ej1 not in explored): #cand[ej1.head]): <= gives duplicates
			prob = getprob(chart, D, ej1)
			# add it to the heap
			cand[ej1.head][ej1] = prob
			explored.add(ej1)

def getprob(chart, D, ej):
	e = ej.edge
	eleft = ChartItem(e.rule.rhs1, ej.head.left, e.split)
	if eleft in D: # and ej.left < len(D[eleft]):
		entry = D[eleft][ej.left]; prob = entry.value
	elif ej.left == 0:
		if chart[eleft.left][eleft.right].get(eleft.label, False):
			edge = chart[ej.head.left][e.split][e.rule.rhs1][0]
			prob = edge.inside
		else:
			prob = np.inf
			print "not found left:", ej
	else: raise ValueError("non-zero rank vector not part of explored derivations")
	result = e.rule.prob + prob
	if ej.right >= 0: #if e.right.label:
		eright = ChartItem(e.rule.rhs2, e.split, ej.head.right)
		if eright in D: # and ej.right < len(D[eright]):
			entry = D[eright][ej.right]; prob = entry.value
		elif ej.right == 0:
			if chart[eright.left][eright.right].get(eright.label, False):
				edge = chart[e.split][ej.head.right][e.rule.rhs2][0]
				prob = edge.inside
			else:
				prob = np.inf
				print "not found right:", ej

		else: raise ValueError("non-zero rank vector not part of explored derivations")
		result += prob
	return result

def getderivation(ej, D, chart, tolabel, sent, n):
	""" Translate the (e, j) notation to an actual tree string in
	bracket notation.  e is an edge, j is a vector prescribing the rank of the
	corresponding tail node. For example, given the edge <S, [NP, VP], 1.0> and
	vector [2, 1], this points to the derivation headed by S and having the 2nd
	best NP and the 1st best VP as children.
	"""
	if n > 100: return ""	#hardcoded limit to prevent cycles
	e = ej.edge
	children = []
	eleft = ChartItem(e.rule.rhs1, ej.head.left, e.split)
	eright = ChartItem(e.rule.rhs2, e.split, ej.head.right)
	for ei, i in ((eleft, ej.left), (eright, ej.right)):
		if i == -1: break
		if chart[ei.left][ei.right].get(ei.label, False):
			if ei in D:
				entry = D[ei][i]
				children.append(
					getderivation(entry.key, D, chart, tolabel, sent, n + 1))
			else:
				if i == 0:
					edge = chart[ei.left][ei.right][ei.label][0]
					children.append(getderivation(
						RankedEdge(ei, edge, 0, 0 if edge.rule.rhs2 else -1),
						D, chart, tolabel, sent, n + 1))
				else: raise ValueError("non-best edge missing in derivations")
		else:
			# this must be a terminal
			children.append(sent[ei.left])

	if "" in children: return ""
	return "(%s %s)" % (tolabel[ej.head.label], "".join(children))

def lazykbest(chart, goal, k, tolabel, sent):
	""" wrapper function to run lazykthbest and get the actual derivations.
	chart is a monotone hypergraph; should be acyclic unless probabilities
	resolve the cycles (maybe nonzero weights for unary productions are
	sufficient?).
	maps ChartItems to lists of tuples with ChartItems and a weight. The
	items in each list are to be ordered as they were added by the viterbi
	parse, with the best item first.
	goal is a ChartItem that is to be the root node of the derivations.
	k is the number of derivations desired.
	tolabel is a dictionary mapping numeric IDs to the original nonterminal
	labels.  """
	import logging
	D = {}
	cand = {}
	explored = set()
	assert k # sanity check; don't ask for zero derivations
	lazykthbest(goal, k, k, D, cand, chart, explored)
	#for a in D:
	#	print a.left, a.right
	#	print tolabel[a.label] if a.label in tolabel else a.label
	#	for entry in D[a]:
	#		print entry.key, entry.value
	logging.debug("%d (sub)derivations considered" % len(explored))
	return filter(itemgetter(0), [
			(getderivation(entry.key, D, chart, tolabel, sent, 0), entry.value)
			for entry in D[goal]])

