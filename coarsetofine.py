"""Assorted functions to project items from a coarse chart to corresponding
items for a fine grammar.
"""
import re, logging
from collections import defaultdict
import numpy as np
from nltk import Tree
try: dictcast({})
except NameError:
	from containers import * #ChartItem, Edge, dictcast, itemcast, edgecast
	from kbest import * #lazykthbest 

infinity = float('infinity')
removeids = re.compile("@[0-9]+")
removeparentannot = re.compile("\^<.*>")
reducemarkov = [re.compile("\|<[^>]*>"),
				re.compile("\|<([^->]*)-?[^>]*>"),
				re.compile("\|<([^->]*-[^->]*)-?[^>]*>"),
				re.compile("\|<([^->]*-[^->]*-[^->]*)-?[^>]*>")]

def whitelistfromchart(chart, goal, coarse, fine, k, whitelist, maxlen):
	""" Produce a white list of chart items occurring in the k-best derivations
	of chart, where labels X in the coarse grammar are projected to the labels
	X and X@n in 'toid', for possible values of n.  When k==0, the chart is
	merely filtered to contain only items that contribute to a complete
	derivation."""

	#l = [{} for a in coarse.toid]
	#for a, label in fine.toid.iteritems():
	#	for left, right in l[coarse.toid.get(removeids.sub("", a), 1)]:
	#		whitelist[left, right, label] = np.inf
	l = np.array([np.NAN], dtype='d').repeat(
		len(coarse.toid) * maxlen * (maxlen+1)).reshape(
		(len(coarse.toid), maxlen, (maxlen+1)))
	# construct a table mapping each nonterminal A or A@x
	# to the outside score for A in the chart to prune with
	kbest = kbest_outside(chart, goal, k)
	# uses ids of labels in coarse chart
	for Ih in kbest:
		l[itemcast(Ih).label, itemcast(Ih).left,
			itemcast(Ih).right] = np.inf #kbest[Ih]
	for a, label in fine.toid.iteritems():
		whitelist[label] = l[coarse.toid[removeids.sub("", a)]]
	#logging.debug('pruning with %d nonterminals, %d items' % (
	#	len(filter(None, whitelist)), len(kbest)))

def kbest_outside(chart, start, k):
	""" produce a dictionary of ChartItems with the best outside score
	according to the k-best derivations in a chart. """
	D = {}
	outside = { start : 0.0 }
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

def getitems(ej, rootprob, D, chart, outside):
	""" Traverse a derivation e,j, noting outside costs relative to its root
	edge """
	e = ej.edge
	eleft = new_ChartItem(e.rule.rhs1, ej.head.left, e.split)
	if eleft in chart:
		if eleft in D:
			entry = D[eleft][ej.left]
			eejj = entry.key; prob = entry.value
		elif ej.left == 0:
			eejj = RankedEdge(eleft, min(chart[eleft]), 0, 0)
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
			eejj = RankedEdge(eright, min(chart[eright]), 0, 0)
			prob = eejj.edge.inside
		else: raise ValueError
		if eright not in outside:
			outside[eright] = rootprob - prob
		getitems(eejj, rootprob, D, chart, outside)

def filterchart(chart, start):
	""" remove all entries that do not contribute to a complete derivation
	headed by "start" """
	chart2 = {}
	filter_subtree(start, dictcast(chart), chart2)
	return chart2

def filter_subtree(start, chart, chart2):
	chart2[start] = chart[start]
	for edge in chart[start]:
		item = new_ChartItem(edge.rule.rhs1, start.i, edge.split)
		if item.label and item not in chart2:
			filter_subtree(item, chart, chart2)
		item = new_ChartItem(edge.rule.rhs2, edge.split, start.right)
		if item.label and item not in chart2:
			filter_subtree(item, chart, chart2)


"""
def doctf(coarse, fine, sent, tree, k, doph, headrules, pa, split,
			verbose=False):
	try: from plcfrs import parse, pprint_chart
	except ImportError: from oldplcfrs import parse, pprint_chart
	#from coarsetofine import kbest_outside, merged_kbest, prunelist_fromchart
	from disambiguation import marginalize
	from treetransforms import mergediscnodes, un_collinize
	from containers import getlabel, getvec
	from grammar import canonicalize, rem_marks
	from math import exp
	sent, tags = zip(*sent)
	print " C O A R S E ",
	p, start = parse(sent, coarse, start=coarse.toid['ROOT'], tags=tags)
	if start:
		mpp = marginalize(p, start, coarse.tolabel)
		for t in mpp:
			print exp(-mpp[t]),
			t = Tree.parse(t, parse_leaf=int)
			#print t.pprint(margin=999)
			un_collinize(t)
			if split: mergediscnodes(t)
			t = canonicalize(rem_marks(t))
			print "exact match" if t == canonicalize(tree) else "no match"
	else:
		print "no parse"
		return
		pprint_chart(p, sent, coarse.tolabel)
	l = prunelist_fromchart(p, start, coarse, fine, k, #0 if split else k,
				removeparentannotation=pa, mergesplitnodes=False,
				reduceh=doph)
	if verbose:
		print "\nitems in 50-best of coarse chart"
		if split:
			d = merged_kbest(p, start, k, coarse)
			for label in d:
				print label, map(bin, d[label].keys())
		else:
			kbest = kbest_outside(p, start, k)
			for a,b in kbest.items():
				print coarse.tolabel[getlabel(a)], bin(getvec(a)), b
		print "\nprunelist:"
		for n,x in enumerate(l):
			print fine.tolabel[n], [(bin(v), s) for v,s in x.items()]
	print " F I N E ",
	pp, start = parse(sent, fine, start=fine.toid['ROOT'], tags=tags, prunelist=None, neverblockmarkovized=pa, neverblockdiscontinuous=False, splitprune=split)
	if start:
		mpp = marginalize(pp, start, fine.tolabel)
		for t in mpp:
			print exp(-mpp[t]),
			t = Tree.parse(t, parse_leaf=int)
			un_collinize(t)
			t = canonicalize(rem_marks(t))
			#print t.pprint(margin=999)
			print "exact match" if t == canonicalize(tree) else "no match"
	else:
		print "no parse. problem."; return
		#if verbose:
		xp = set((coarse.tolabel[a.label], a.vec) for a in p.keys() if p[a])
		xpp = set((fine.tolabel[a.label], a.vec) for a in pp.keys() if pp[a])
		print "difference:"
		for a in xp - xpp:
			if "*" not in a[0]: print a[0], bin(a[1])
		print "\nfine chart:"
		for a in xpp:
			if "@" not in a[0]: print a[0], bin(a[1])

		#pprint_chart(pp, sent, fine.tolabel)

def main():
	from treetransforms import splitdiscnodes, collinize
	from negra import NegraCorpusReader, readheadrules
	from grammar import splitgrammar, induce_srcg, dop_srcg_rules,\
			printrule, subsetgrammar
	headrules = readheadrules()
	k = 50
	#corpus = NegraCorpusReader(".", "sample2\.export",
	#	encoding="iso-8859-1", headorder=True,
	#	headfinal=True, headreverse=False)
	corpus = NegraCorpusReader("../rparse", "negraproc\.export",
		encoding="utf-8", headorder=True, headfinal=True, headreverse=False)
	train = 400
	test = 40
	testmaxlen = 15
	trees = list(corpus.parsed_sents()[:train])
	sents = corpus.sents()[:train]
	
	dtrees = [t.copy(True) for t in trees]
	parenttrees = [t.copy(True) for t in trees]
	for t in trees: collinize(t, vertMarkov=0)
	cftrees = [splitdiscnodes(t.copy(True), markorigin=True) for t in trees]
	for t in cftrees:
			#t.chomsky_normal_form(childChar=":")
			collinize(t, horzMarkov=999, tailMarker='', leftMostUnary=True, childChar=":")
	for t in parenttrees: collinize(t, vertMarkov=2)
	for t in dtrees: collinize(t, vertMarkov=0, horzMarkov=1)
	# mark heads, canonicalize, binarize head outward
	normalsrcg = induce_srcg(trees, sents)
	normal = splitgrammar(normalsrcg)
	parent = splitgrammar(induce_srcg(parenttrees, sents))
	split = splitgrammar(induce_srcg(cftrees, sents))
	for t,s in zip(cftrees, sents):
		for (r,yf),w in induce_srcg([t], [s]): assert len(yf) == 1
	fine999srcg = dop_srcg_rules(trees, sents)
	fine999 = splitgrammar(fine999srcg)
	fine1 = splitgrammar(dop_srcg_rules(dtrees, sents))
	trees = list(corpus.parsed_sents()[train:train+test])
	sents = corpus.tagged_sents()[train:train+test]
	if subsetgrammar(normalsrcg, fine999srcg):
		print "DOP grammar is superset"
	else: print "DOP grammar is NOT a superset!"
	for msg, coarse, fine, settings in zip(
		"normal parentannot cf-split".split(),
		(normal, parent, split),
		(fine999, fine1, fine1),
		((False, False), (True, False), (False, True))):
		print "coarse grammar:", msg
		for n, (sent, tree) in enumerate(zip(sents, trees)):
			if len(sent) > testmaxlen: continue
			print n,
			doctf(coarse, fine, sent, tree, k,
					1 if msg == "parentannot" else 999, headrules, *settings)

if __name__ == '__main__': main()
"""
