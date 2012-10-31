"""Assorted functions to project items from a coarse chart to corresponding
items for a fine grammar.
"""
import re, logging
from collections import defaultdict
from math import exp, log
import numpy as np
from nltk import Tree
from agenda import Agenda
from containers import ChartItem, Edge, RankedEdge
from kbest import lazykthbest, lazykbest

cimport cython
cimport numpy as np
from agenda cimport Entry, Agenda
from containers cimport ChartItem, Edge, RankedEdge, Grammar, new_ChartItem

infinity = float('infinity')
removeids = re.compile("@[0-9_]+")
removeparentannot = re.compile("\^<.*>")
reducemarkov = [re.compile("\|<[^>]*>"),
				re.compile("\|<([^->]*)-?[^>]*>"),
				re.compile("\|<([^->]*-[^->]*)-?[^>]*>"),
				re.compile("\|<([^->]*-[^->]*-[^->]*)-?[^>]*>")]

def whitelistfromposteriors2(np.ndarray[np.double_t, ndim=3] inside,
	np.ndarray[np.double_t, ndim=3] outside, ChartItem goal, Grammar coarse,
	Grammar fine, dict mapping, short maxlen, double threshold):
	""" compute posterior probabilities and prune away cells below some
	threshold. this version is for use with parse_sparse(). """
	cdef long label
	cdef short lensent = goal.right
	sentprob = inside[0, lensent, goal.label]
	print "sentprob=%g" % sentprob,
	posterior = (inside[:lensent, :lensent+1]
		* outside[:lensent, :lensent+1]) / sentprob

	print " ", (posterior > threshold).sum(),
	print "of", (posterior != 0.0).sum(),
	print "nonzero coarse items left"
	finechart = [[{} for _ in range(maxlen)] for _ in range(maxlen)]
	leftidx, rightidx, labels = (posterior[:lensent, :lensent+1]
		> threshold).nonzero()
	#labels, leftidx, rightidx = (posterior[:,:lensent,:lensent+1] > threshold).nonzero()
	for label, left, right in zip(labels, leftidx, rightidx):
		finechart[left][right].update((lhs, []) for lhs in mapping[label])
	print "copied chart."
	return finechart

def whitelistfromposteriors(np.ndarray[np.double_t, ndim=3] inside,
	np.ndarray[np.double_t, ndim=3] outside, ChartItem goal, Grammar coarse,
	Grammar fine, np.ndarray[np.double_t, ndim=3] finechart, short maxlen,
	double threshold):
	""" compute posterior probabilities and prune away cells below some
	threshold. this version produces a matrix with pruned spans having NaN as
	value. """
	cdef long label
	cdef short lensent = goal.right
	sentprob = inside[0, lensent, goal.label]
	print "sentprob=%g" % sentprob
	posterior = (inside[:lensent, :lensent+1, :]
			* outside[:lensent, :lensent+1, :]) / sentprob
	inside[:lensent, :lensent + 1, :] = np.NAN
	inside[posterior > threshold] = np.inf
	print " ", (posterior > threshold).sum(),
	print "of", (posterior != 0.0).sum(),
	print "nonzero coarse items left",
	#labels, leftidx, rightidx = (posterior[:lensent, :lensent+1, :]
	#	> threshold).nonzero()
	#for left, right, label in zip(leftidx, rightidx, labels):
	#	for x in mapping[label]:
	#		finechart[left, right, x] = inside[left, right, label]
	for label in range(len(fine.toid)):
		finechart[:lensent, :lensent+1, label] = inside[:lensent,:lensent+1,
			coarse.toid[removeids.sub("", fine.tolabel[label])]]
	print "copied chart."

cpdef whitelistfromposteriors1(
	dict chart,
	np.ndarray[np.double_t, ndim=3] viterbi,
	ChartItem goal,
	Grammar coarse,
	Grammar fine,
	np.ndarray[np.double_t, ndim=3] whitelist,
	double threshold):
	""" compute posterior probabilities and prune away cells below some
	threshold. """
	cdef short lensent = goal.right
	inside = np.zeros((lensent, lensent+1, len(coarse.toid)), dtype='d')
	outside = np.zeros_like(inside)
	outside[0, lensent, goal.label] = 1.0
	insidescores(chart, goal, inside)
	outsidescores(chart, goal, inside, outside)
	sentprob = inside[0, lensent, goal.label]
	posterior = inside * outside / sentprob
	viterbi.fill(np.inf)
	viterbi[posterior < threshold] = np.NAN
	print (posterior >= threshold).sum(), "coarse items left",
	for label, id in fine.toid.iteritems():
		whitelist[id] = viterbi[:, :, coarse.toid[removeids.sub("", label)]]

def whitelistfromkbest(dict chart, ChartItem goal,
	Grammar coarse, Grammar fine, int k,
	np.ndarray[np.double_t, ndim=3] whitelist, int maxlen):
	""" Produce a white list of chart items occurring in the k-best derivations
	of chart, where labels X in the coarse grammar are projected to the labels
	X and X@n in 'toid', for possible values of n.  When k==0, the chart is
	merely filtered to contain only items that contribute to a complete
	derivation."""
	cdef ChartItem Ih
	#l = [{} for a in coarse.toid]
	#for a, label in fine.toid.iteritems():
	#	for left, right in l[coarse.toid.get(removeids.sub("", a), 1)]:
	#		whitelist[left, right, label] = np.inf
	cdef np.ndarray[np.double_t, ndim=3] l = np.empty(
		(len(coarse.toid), maxlen, (maxlen+1)), dtype='d')
	l.fill(np.NAN)
	# construct a table mapping each nonterminal A or A@x
	# to the outside score for A in the chart to prune with
	kbest = kbest_outside(chart, goal, k)
	# uses ids of labels in coarse chart
	for Ih in kbest:
		l[(<ChartItem>Ih).label, (<ChartItem>Ih).left,
			(<ChartItem>Ih).right] = np.inf #kbest[Ih]
	print (l == np.inf).sum(), #"of", len(chart), 
	print "coarse items left"
	for a, label in fine.toid.iteritems():
		whitelist[label] = l[coarse.toid[removeids.sub("", a)]]
	#logging.debug('pruning with %d nonterminals, %d items' % (
	#	len(filter(None, whitelist)), len(kbest)))

cdef dict kbest_outside(dict chart, ChartItem start, int k):
	""" produce a dictionary of ChartItems with the best outside score
	according to the k-best derivations in a chart. """
	cdef Entry entry
	cdef Edge edge
	cdef dict D = {}, outside = { start : 0.0 }
	if k == 0:
		outside = filterchart(chart, start)
		for a in outside:
			# use if probabilities matter
			#e = min(outside[a])
			#outside[a] = e.inside
			outside[a] = 0.0
	else:
		lazykthbest(start, k, k, D, {}, chart, set())
		for entry in D[start]:
			getitems(entry.key, entry.value, D, chart, outside)
	return outside

cdef void getitems(RankedEdge ej, double rootprob, dict D,
		dict chart, dict outside):
	""" Traverse a derivation e,j, noting outside costs relative to its root
	edge """
	cdef Edge e
	cdef Entry entry
	cdef ChartItem, eleft, eright
	cdef RankedEdge eejj
	cdef double prob
	e = ej.edge
	eleft = new_ChartItem(e.rule.rhs1, ej.head.left, e.split)
	if chart[ej.head.left][e.split].get(e.rule.rhs1):
		if eleft in D:
			entry = D[eleft][ej.left]
			eejj = entry.key; prob = entry.value
		elif ej.left == 0:
			eejj = RankedEdge(eleft,
				chart[ej.head.left][e.split][e.rule.rhs1][0], 0, 0)
			prob = eejj.edge.inside
		else: raise ValueError
		if eleft not in outside:
			outside[eleft] = rootprob - prob
		getitems(eejj, rootprob, D, chart, outside)
	if e.rule.rhs2:
		eright = new_ChartItem(e.rule.rhs2, e.split, ej.head.right)
		if eright in D:
			entry = D[eright][ej.right]
			eejj = entry.key; prob = entry.value
		elif ej.right == 0:
			eejj = RankedEdge(eright,
				chart[e.split][ej.head.right][e.rule.rhs2][0], 0, 0)
			prob = eejj.edge.inside
		else: raise ValueError
		if eright not in outside:
			outside[eright] = rootprob - prob
		getitems(eejj, rootprob, D, chart, outside)

def filterchart(chart, start):
	""" remove all entries that do not contribute to a complete derivation
	headed by "start" """
	chart2 = {}
	filter_subtree(start, <dict>chart, chart2)
	return chart2

cdef void filter_subtree(ChartItem start, dict chart, dict chart2):
	cdef ChartItem item
	cdef Edge edge
	chart2[start] = chart[start]
	for edge in chart[start]:
		item = new_ChartItem(edge.rule.rhs1, start.left, edge.split)
		if item.label and item not in chart2:
			filter_subtree(item, chart, chart2)
		item = new_ChartItem(edge.rule.rhs2, edge.split, start.right)
		if item.label and item not in chart2:
			filter_subtree(item, chart, chart2)

def insidescores(dict chart, ChartItem goal, np.ndarray[np.double_t, ndim=3] inside):
	# this is probably incorrect ...
	# I guess edges which share the same children should share probability 
	# mass in the total inside probability
	cdef ChartItem item
	cdef list edges
	for item, edges in chart.iteritems():
		inside[item.left, item.right, item.label] = logsumexp(
				[edge.inside for edge in edges])

def outsidescores(dict chart, ChartItem goal, np.ndarray[np.double_t, ndim=3] inside, np.ndarray[np.double_t, ndim=3] outside):
	cdef Agenda agenda = Agenda()
	cdef ChartItem item, newitem
	cdef Edge edge
	cdef set visited
	seq = 0
	agenda[goal] = seq
	visited = set()
	while agenda:
		item = agenda.popentry().key
		seq += 1
		visited.add(item)
		for edge in chart[item]:
			if edge.rule.rhs2 == 0:
				outside[item.left, edge.split, edge.rule.rhs1] = (
					outside[item.left, edge.split, edge.rule.rhs1]
					+ outside[item.left, item.right, item.label]
					* exp(-edge.rule.prob))
			else:
				outside[item.left, edge.split, edge.rule.rhs1] = (
					outside[item.left, edge.split, edge.rule.rhs1] +
					outside[item.left, item.right, item.label]
					* inside[edge.split, item.right, edge.rule.rhs2]
					* exp(-edge.rule.prob))
				outside[edge.split, item.right, edge.rule.rhs2] = (
					outside[edge.split, item.right, edge.rule.rhs2] +
					outside[item.left, item.right, item.label]
					* inside[item.left, edge.split, edge.rule.rhs1]
					* exp(-edge.rule.prob))
				newitem = ChartItem(edge.rule.rhs2, edge.split, item.right)
				if newitem not in visited: agenda[newitem] = seq
			if edge.rule.rhs1:
				newitem = ChartItem(edge.rule.rhs1, item.left, edge.split)
				if newitem not in visited: agenda[newitem] = seq
			# do we need to order the agenda more specifically?

def logsumexp(logprobs):
	#NB: this expects a list of negative log probabilities and
	#returns a normal probability in the interval [0-1].
	maxprob = min(logprobs)
	return sum([exp(maxprob - prob) for prob in logprobs]) / exp(maxprob)
