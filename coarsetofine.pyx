"""Assorted functions to project items from a coarse chart to corresponding
items for a fine grammar.
"""
import re, logging
from sys import stderr
from collections import defaultdict
from math import exp, log
import numpy as np
from nltk import Tree
from agenda import Agenda
from containers import ChartItem, Edge, RankedEdge
from kbest import lazykthbest

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
	np.ndarray[np.double_t, ndim=3] outside, start, lensent, Grammar coarse,
	Grammar fine, dict mapping, double threshold):
	""" compute posterior probabilities and prune away cells below some
	threshold. this version is for use with parse_sparse(). """
	cdef long label
	sentprob = inside[0, lensent, start]
	print >>stderr, "sentprob=%g" % sentprob,
	posterior = (inside[:lensent, :lensent+1]
		* outside[:lensent, :lensent+1]) / sentprob

	print >>stderr, " ", (posterior > threshold).sum(),
	print >>stderr, "of", (posterior != 0.0).sum(),
	print >>stderr, "nonzero coarse items left"
	finechart = [[{} for _ in range(lensent+1)] for _ in range(lensent)]
	leftidx, rightidx, labels = (posterior[:lensent, :lensent+1]
		> threshold).nonzero()
	#labels, leftidx, rightidx = (posterior[:,:lensent,:lensent+1] > threshold).nonzero()
	for label, left, right in zip(labels, leftidx, rightidx):
		finechart[left][right].update((lhs, []) for lhs in mapping[label])
	print >>stderr, "copied chart."
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
	print >>stderr, "sentprob=%g" % sentprob
	posterior = (inside[:lensent, :lensent+1, :]
			* outside[:lensent, :lensent+1, :]) / sentprob
	inside[:lensent, :lensent + 1, :] = np.NAN
	inside[posterior > threshold] = np.inf
	print >>stderr, " ", (posterior > threshold).sum(),
	print >>stderr, "of", (posterior != 0.0).sum(),
	print >>stderr, "nonzero coarse items left",
	#labels, leftidx, rightidx = (posterior[:lensent, :lensent+1, :]
	#	> threshold).nonzero()
	#for left, right, label in zip(leftidx, rightidx, labels):
	#	for x in mapping[label]:
	#		finechart[left, right, x] = inside[left, right, label]
	for label in range(len(fine.toid)):
		finechart[:lensent, :lensent+1, label] = inside[:lensent,:lensent+1,
			coarse.toid[removeids.sub("", fine.tolabel[label])]]
	print >>stderr, "copied chart."

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
	print >>stderr, (posterior >= threshold).sum(), "coarse items left",
	for label, id in fine.toid.iteritems():
		whitelist[id] = viterbi[:, :, coarse.toid[removeids.sub("", label)]]

def whitelistfromkbest(list chart, start, lensent,
	Grammar coarse, Grammar fine, int k, dict mapping):
	""" Produce a white list of chart items occurring in the k-best derivations
	of chart, where labels X in the coarse grammar are projected to the labels
	X and X@n in 'toid', for possible values of n.  When k==0, the chart is
	merely filtered to contain only items that contribute to a complete
	derivation."""
	cdef Entry entry
	cdef Edge edge
	cdef list D, kbestspans = [[{} for _ in range(lensent+1)]
			for _ in range(lensent)]
	kbestspans[0][lensent][start] = 0.0
	if k == 0:
		raise NotImplemented
		#kbestspans = filterchart(chart, goal)
		#for a in kbestspans:
		#	# use if probabilities matter
		#	#e = min(kbestspans[a])
		#	#kbestspans[a] = e.inside
		#	kbestspans[a] = 0.0
	else:
		cand = [[{} for _ in x] for x in chart]
		D = [[{} for _ in x] for x in chart]
		lazykthbest(start, 0, lensent, k, k, D, cand, chart, set())
		for entry in D[0][lensent][start]:
			getitems(entry.key, entry.value, D, chart, kbestspans)

	print >>stderr, ''
	#print >>stderr, (kbestspans == np.inf).sum(), #"of", len(chart), 
	#print >>stderr, "coarse items left"
	finechart = [[{} for _ in range(lensent+1)] for _ in range(lensent)]
	for left in range(lensent):
		for right in range(left + 1, lensent+1):
			for label in kbestspans[left][right]:
				finechart[left][right].update((finelabel, [])
					for finelabel in mapping[label])
	#logging.debug('pruning with %d nonterminals, %d items' % (
	#	len(filter(None, whitelist)), len(kbest)))
	return finechart

cdef void getitems(RankedEdge ej, double rootprob, list D,
		list chart, list items):
	""" Traverse a derivation e,j, noting outside costs relative to its root
	edge """
	cdef Edge e
	cdef Entry entry
	cdef RankedEdge eejj
	cdef double prob
	e = ej.edge
	if chart[ej.left][e.split].get(e.rule.rhs1):
		if e.rule.rhs1 in D[ej.left][e.split]:
			entry = D[ej.left][e.split][e.rule.rhs1][ej.leftrank]
			eejj = entry.key; prob = entry.value
		elif ej.left == 0:
			eejj = RankedEdge(e.rule.rhs1, ej.left, e.split,
				chart[ej.left][e.split][e.rule.rhs1][0], 0, 0)
			prob = eejj.edge.inside
		else: raise ValueError
		if e.rule.rhs1 not in items[ej.left][e.split]:
			items[ej.left][e.split][e.rule.rhs1] = rootprob - prob
		getitems(eejj, rootprob, D, chart, items)
	if e.rule.rhs2:
		if e.rule.rhs2 in D[e.split][ej.right]:
			entry = D[e.split][ej.right][e.rule.rhs2][ej.rightrank]
			eejj = entry.key; prob = entry.value
		elif ej.rightrank == 0:
			eejj = RankedEdge(e.rule.rhs2, e.split, ej.right,
				chart[e.split][ej.right][e.rule.rhs2][0], 0, 0)
			prob = eejj.edge.inside
		else: raise ValueError
		if e.rule.rhs2 not in items[e.split][ej.right]:
			items[e.split][ej.right][e.rule.rhs2] = rootprob - prob
		getitems(eejj, rootprob, D, chart, items)

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
