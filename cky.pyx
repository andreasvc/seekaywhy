""" Probabilistic CKY parser for PCFGs """
import codecs
from sys import stdout, stderr
from collections import defaultdict
from nltk import Tree
import numpy as np

from agenda cimport EdgeAgenda, Entry
cimport numpy as np
from containers cimport Edge, Rule, Terminal, Grammar, new_Edge

cdef extern from "math.h":
	bint isinf(double)
	bint isnan(double)
	bint isfinite(double)
	double log(double)
	double exp(double)

encoding = "utf-8"

def parse(list sent, Grammar grammar, whitelist):
	""" A CKY parser modeled after Bodenstab's `fast grammar loop.'
		and the Stanford parser. """
	cdef short left, right, mid, span, lensent = len(sent)
	cdef short narrowl, narrowr, minmid, maxmid
	cdef long numsymbols = len(grammar.toid), lhs, rhs1
	cdef double prob
	cdef bint foundbetter, newspan
	cdef Rule rule
	cdef Terminal terminal
	cdef EdgeAgenda unaryagenda = EdgeAgenda()
	cdef list chart = [[{} for _ in range(lensent+1)] for _ in range(lensent)]
	cdef dict cell
	cdef unicode word
	# the viterbi chart is initially filled with infinite log probabilities,
	# spans which are to be blocked contain NaN.
	cdef np.ndarray[np.double_t, ndim=3] viterbi
	# matrices for the filter which gives minima and maxima for splits
	cdef np.ndarray[np.int16_t, ndim=2] minleft, maxleft, minright, maxright
	minleft = np.empty((numsymbols, lensent+1), dtype='int16'); minleft.fill(-1)
	maxleft = np.empty_like(minleft); maxleft.fill(lensent+1)
	minright = np.empty_like(minleft); minright.fill(lensent+1)
	maxright = np.empty_like(minleft); maxright.fill(-1)
	assert grammar.logprob, "Expecting grammar with log probabilities."

	if whitelist is None:
		viterbi = np.empty((lensent, lensent + 1, numsymbols), dtype='d')
		viterbi.fill(np.inf)
	else: viterbi = whitelist
	# assign POS tags
	print >>stderr, 1, # == span
	for left in range(lensent):
		right = left + 1
		cell = chart[left][right]
		word = unicode(sent[left]) if sent[left] in grammar.lexicon else u""
		for terminal in <list>grammar.lexical:
			if (not isnan(viterbi[left, right, terminal.lhs]) and
				terminal.word == word):
				lhs = terminal.lhs
				viterbi[left, right, lhs] = terminal.prob
				cell[lhs] = [new_Edge(terminal.prob, terminal, right)]
				# update filter
				if left > minleft[lhs, right]: minleft[lhs, right] = left
				if left < maxleft[lhs, right]: maxleft[lhs, right] = left
				if right < minright[lhs, left]: minright[lhs, left] = right
				if right > maxright[lhs, left]: maxright[lhs, left] = right

		# unary rules on the span of this POS tag
		# NB: for this agenda, only the probabilities of the edges matter
		unaryagenda.update([(
			rhs1, new_Edge(viterbi[left, right, rhs1], None, 0))
			for rhs1 in range(numsymbols)
			if isfinite(viterbi[left, right, rhs1])])
		while unaryagenda.length:
			rhs1 = unaryagenda.popentry().key
			for rule in <list>grammar.unarybyrhs[rhs1]:
				if isnan(viterbi[left, right, rule.lhs]): continue
				lhs = rule.lhs
				prob = rule.prob + viterbi[left, right, rhs1]
				if not isfinite(viterbi[left, right, lhs]): cell[lhs] = []
				edge = new_Edge(prob, rule, right)
				cell[lhs].append(edge)
				if prob >= viterbi[left, right, lhs]: continue
				unaryagenda.setitem(lhs, edge)
				viterbi[left, right, lhs] = prob
				# update filter
				if left > minleft[lhs, right]: minleft[lhs, right] = left
				if left < maxleft[lhs, right]: maxleft[lhs, right] = left
				if right < minright[lhs, left]: minright[lhs, left] = right
				if right > maxright[lhs, left]: maxright[lhs, left] = right

	for span in range(2, lensent + 1):
		print >>stderr, span,
		stdout.flush()

		# loop over all non-pruned indices. this appears to be slow.
		#labels, leftidx, rightidx = np.isinf(
		#		viterbi[:lensent,:lensent+1,:]).nonzero()
		#indices = (rightidx - leftidx).argsort()
		#prevspan = 0
		#for idx in indices:
		#	lhs = labels[idx]
		#	left = leftidx[idx]
		#	right = rightidx[idx]
		#	span = right - left
		#	if span > prevspan:
		#		prevspan = span
		#		print >>stderr, span,
		#		stdout.flush()
		#	for rule in <list>(grammar.binary[lhs]):

		# constituents from left to right
		for left in range(0, lensent - span + 1):
			right = left + span
			cell = chart[left][right]
			# binary rules 
			for lhs, rules in enumerate(<list>grammar.binary):
				if isnan(viterbi[left, right, lhs]): continue
				for rule in <list>rules:
					#if not (np.isfinite(viterbi[left, left+1:right, rule.rhs1]
					#	).any() and np.isfinite(viterbi[
					#	left:right-1, right, rule.rhs2]).any()):
					#	continue
					narrowr = minright[rule.rhs1, left]
					narrowl = minleft[rule.rhs2, right]
					if narrowr >= right or narrowl < narrowr: continue
					minmid = smax(narrowr, maxleft[rule.rhs2, right])
					maxmid = smin(maxright[rule.rhs1, left], narrowl) + 1
					newspan = isinf(viterbi[left, right, lhs])
					foundbetter = False
					for mid in range(minmid, maxmid):
						if (isfinite(viterbi[left, mid, rule.rhs1])
							and isfinite(viterbi[mid, right, rule.rhs2])):
							prob = (rule.prob + viterbi[left, mid, rule.rhs1]
									+ viterbi[mid, right, rule.rhs2])
							if not isfinite(viterbi[left, right, lhs]):
								cell[lhs] = []
							edge = new_Edge(prob, rule, mid)
							cell[lhs].append(edge)
							if prob < viterbi[left, right, lhs]:
								foundbetter = True
								viterbi[left, right, lhs] = prob
					# update filter
					if not foundbetter or not newspan: continue
					if left > minleft[lhs, right]: minleft[lhs, right] = left
					if left < maxleft[lhs, right]: maxleft[lhs, right] = left
					if right < minright[lhs, left]: minright[lhs, left] = right
					if right > maxright[lhs, left]: maxright[lhs, left] = right

			# unary rules on this span
			unaryagenda.update([(rhs1, new_Edge(
					viterbi[left, right, rhs1], None, 0))
				for rhs1 in range(numsymbols)
				if isfinite(viterbi[left, right, rhs1])])
			while unaryagenda.length:
				rhs1 = unaryagenda.popentry().key
				for rule in <list>grammar.unarybyrhs[rhs1]:
					if isnan(viterbi[left, right, rule.lhs]): continue
					lhs = rule.lhs
					prob = rule.prob + viterbi[left, right, rhs1]
					if not isfinite(viterbi[left, right, lhs]): cell[lhs] = []
					edge = new_Edge(prob, rule, right)
					cell[lhs].append(edge)
					if prob >= viterbi[left, right, lhs]: continue
					viterbi[left, right, lhs] = prob
					unaryagenda.setitem(lhs, edge)
					# update filter
					if left > minleft[lhs, right]: minleft[lhs, right] = left
					if left < maxleft[lhs, right]: maxleft[lhs, right] = left
					if right < minright[lhs, left]: minright[lhs, left] = right
					if right > maxright[lhs, left]: maxright[lhs, left] = right
	print >>stderr, ''
	return chart, viterbi

def parse_sparse(list sent, Grammar grammar, chart):
	""" A CKY parser modeled after Bodenstab's `fast grammar loop.'
		and the Stanford parser.
	This version keeps the viterbi probabilities and the rest of chart
	in a single hash table, useful for large grammars. """
	cdef short left, right, mid, span, lensent = len(sent)
	cdef short narrowl, narrowr, minmid, maxmid
	cdef long numsymbols = len(grammar.toid), lhs
	cdef double oldscore, prob, infinity = float('infinity')
	cdef unicode word
	cdef dict cell
	cdef bint foundbetter, newspan
	cdef Rule rule
	cdef Terminal terminal
	cdef EdgeAgenda unaryagenda = EdgeAgenda()
	cdef Entry entry
	# matrices for the filter which gives minima and maxima for splits
	cdef np.ndarray[np.int16_t, ndim=2] minleft, maxleft, minright, maxright
	minleft = np.empty((numsymbols, lensent+1), dtype='int16'); minleft.fill(-1)
	maxleft = np.empty_like(minleft); maxleft.fill(lensent+1)
	minright = np.empty_like(minleft); minright.fill(lensent+1)
	maxright = np.empty_like(minleft); maxright.fill(-1)
	assert grammar.logprob, "Expecting grammar with log probabilities."
	if chart is None:
		chart = [[{} for _ in range(lensent)] for _ in range(lensent)]
	# assign POS tags
	print >>stderr, 1, # == span
	for left in range(lensent):
		right = left + 1
		cell = chart[left][right]
		word = unicode(sent[left]) if sent[left] in grammar.lexicon else u""
		for terminal in <list>grammar.lexical:
			if terminal.lhs in cell and terminal.word == word:
				lhs = terminal.lhs
				cell[lhs] = [new_Edge(terminal.prob, terminal, right)]
				# update filter
				if left > minleft[lhs, right]: minleft[lhs, right] = left
				if left < maxleft[lhs, right]: maxleft[lhs, right] = left
				if right < minright[lhs, left]: minright[lhs, left] = right
				if right > maxright[lhs, left]: maxright[lhs, left] = right
		# unary rules on the span of this POS tag
		# NB: for this agenda, only the probabilities of the edges matter
		unaryagenda.update([(rhs1, edges[0]) for rhs1, edges in cell.iteritems() if edges])
		while unaryagenda.length:
			entry = unaryagenda.popentry()
			for rule in <list>grammar.unarybyrhs[entry.key]:
				if rule.lhs not in cell: continue
				lhs = rule.lhs
				prob = rule.prob + (<Edge>entry.value).inside
				edge = new_Edge(prob, rule, right)
				if cell[lhs]:
					cell[lhs].append(edge)
					if prob >= (<Edge>cell[lhs][0]).inside: continue
					unaryagenda.setitem(lhs, edge)
					# switch previous best & new best
					#cell[lhs][0], cell[lhs][-1] = cell[lhs][-1], <Edge>entry.value
					cell[lhs][0], cell[lhs][-1] = cell[lhs][-1], cell[lhs][0]
					continue
				cell[lhs] = [edge]
				# update filter
				if left > minleft[lhs, right]: minleft[lhs, right] = left
				if left < maxleft[lhs, right]: maxleft[lhs, right] = left
				if right < minright[lhs, left]: minright[lhs, left] = right
				if right > maxright[lhs, left]: maxright[lhs, left] = right
	for span in range(2, lensent + 1):
		print >>stderr, span,
		stdout.flush()

		# constituents from left to right
		for left in range(0, lensent - span + 1):
			right = left + span
			cell = chart[left][right]
			# binary rules 
			for lhs, rules in enumerate(<list>grammar.binary):
				if lhs not in cell: continue
				for rule in <list>rules:
					narrowr = minright[rule.rhs1, left]
					narrowl = minleft[rule.rhs2, right]
					if narrowr >= right or narrowl < narrowr: continue
					minmid = smax(narrowr, maxleft[rule.rhs2, right])
					maxmid = smin(maxright[rule.rhs1, left], narrowl) + 1
					newspan = not cell[lhs]
					foundbetter = False
					for mid in range(minmid, maxmid):
						if (chart[left][mid].get(rule.rhs1)
							and chart[mid][right].get(rule.rhs2)):
							prob = (rule.prob
								+ (<Edge>chart[left][mid][rule.rhs1][0]).inside
								+ (<Edge>chart[mid][right][rule.rhs2][0]).inside)
							if not cell[lhs] or prob < (<Edge>cell[lhs][0]).inside:
								foundbetter = True
								cell[lhs].append(new_Edge(prob, rule, mid))
								# switch previous best & new best
								cell[lhs][0], cell[lhs][-1] = (
										cell[lhs][-1], cell[lhs][0])
							else: cell[lhs] = [new_Edge(prob, rule, mid)]
					# update filter
					if not foundbetter or not newspan: continue
					if left > minleft[lhs, right]: minleft[lhs, right] = left
					if left < maxleft[lhs, right]: maxleft[lhs, right] = left
					if right < minright[lhs, left]: minright[lhs, left] = right
					if right > maxright[lhs, left]: maxright[lhs, left] = right

			# unary rules
			unaryagenda.update([(rhs1, edges[0]) for rhs1, edges in cell.iteritems() if edges])
			while unaryagenda.length:
				entry = unaryagenda.popentry()
				for rule in <list>grammar.unarybyrhs[entry.key]:
					if rule.lhs not in cell: continue
					lhs = rule.lhs
					prob = rule.prob + (<Edge>entry.value).inside
					edge = new_Edge(prob, rule, right)
					if cell[lhs]:
						cell[lhs].append(edge)
						if prob >= (<Edge>cell[lhs][0]).inside: continue
						unaryagenda.setitem(lhs, edge)
						# switch previous best & new best
						#cell[lhs][0], cell[lhs][-1] = cell[lhs][-1], <Edge>entry.value
						cell[lhs][0], cell[lhs][-1] = cell[lhs][-1], cell[lhs][0]
						continue
					cell[lhs] = [edge]
					# update filter
					if left > minleft[lhs, right]: minleft[lhs, right] = left
					if left < maxleft[lhs, right]: maxleft[lhs, right] = left
					if right < minright[lhs, left]: minright[lhs, left] = right
					if right > maxright[lhs, left]: maxright[lhs, left] = right
	print >>stderr, ''
	return chart

def doinsideoutside(list sent, Grammar grammar, inside, outside):
	start = grammar.toid["TOP"]
	lensent = len(sent)
	numsymbols = len(grammar.toid)
	if inside == None:
		inside = np.zeros((lensent, lensent+1, numsymbols), dtype='d')
	else: inside[:len(sent), :len(sent)+1, :] = 0.0
	if outside == None:
		outside = np.zeros((lensent, lensent+1, numsymbols), dtype='d')
	else: outside[:len(sent), :len(sent)+1, :] = 0.0
	minmaxlr = insidescores(sent, grammar, inside)
	outside = outsidescores(grammar, start, lensent, inside, outside, *minmaxlr)
	return inside, outside

def insidescores(list sent, Grammar grammar,
	np.ndarray[np.double_t, ndim=3] inside):
	""" Compute (viterbi?) inside scores. """
	cdef short left, right, mid, span, lensent = len(sent)
	cdef short narrowl, narrowr, minmid, maxmid
	cdef long numsymbols = len(grammar.toid), lhs
	cdef double oldscore, prob, ls, rs, ins
	cdef bint foundbetter = False
	cdef Rule rule
	cdef Terminal terminal
	cdef unicode word
	# matrices for the filter which give minima and maxima for splits
	cdef np.ndarray[np.int16_t, ndim=2] minleft, maxleft, minright, maxright
	minleft = np.empty((numsymbols, lensent+1), dtype='int16'); minleft.fill(-1)
	maxleft = np.empty_like(minleft); maxleft.fill(lensent+1)
	minright = np.empty_like(minleft); minright.fill(lensent+1)
	maxright = np.empty_like(minleft); maxright.fill(-1)
	inside[:lensent, :lensent+1, :] = 0.0
	assert not grammar.logprob, "Grammar must not have log probabilities."
	print >>stderr, "inside ",
	# assign POS tags
	for left in range(lensent):
		right = left + 1
		word = unicode(sent[left]) if sent[left] in grammar.lexicon else u""
		for terminal in <list>grammar.lexical:
			if terminal.word == word:
				lhs = terminal.lhs
				inside[left, right, lhs] = terminal.prob
				# update filter
				if left > minleft[lhs, right]: minleft[lhs, right] = left
				if left < maxleft[lhs, right]: maxleft[lhs, right] = left
				if right < minright[lhs, left]: minright[lhs, left] = right
				if right > maxright[lhs, left]: maxright[lhs, left] = right
		# unary rules on POS tags 
		for rule in <list>grammar.unary:
			if inside[left, right, rule.rhs1 ] != 0.0:
				lhs = rule.lhs
				prob = rule.prob * inside[left, right, rule.rhs1]
				inside[left, right, lhs] += prob
				# update filter
				if left > minleft[lhs, right]: minleft[lhs, right] = left
				if left < maxleft[lhs, right]: maxleft[lhs, right] = left
				if right < minright[lhs, left]: minright[lhs, left] = right
				if right > maxright[lhs, left]: maxright[lhs, left] = right
	for span in range(1, lensent + 1):
		print >>stderr, span,
		stdout.flush()
		# constituents from left to right
		for left in range(0, lensent - span + 1):
			right = left + span
			# binary
			for lhs, rules in enumerate(grammar.binary):
				for rule in rules:
					narrowr = minright[rule.rhs1, left]
					narrowl = minleft[rule.rhs2, right]
					if narrowr >= right or narrowl < narrowr: continue
					minmid = smax(narrowr, maxleft[rule.rhs2, right])
					maxmid = smin(maxright[rule.rhs1, left], narrowl) + 1
					#oldscore = inside[left, right, lhs]
					foundbetter = False
					for split in range(minmid, maxmid):
						ls = inside[left, split, rule.rhs1]
						if ls == 0.0: continue
						rs = inside[split, right, rule.rhs2]
						if rs == 0.0: continue
						foundbetter = True
						inside[left, right, rule.lhs] += rule.prob * ls * rs
					if foundbetter: #and oldscore == 0.0:
						if left > minleft[lhs,right]: minleft[lhs,right] = left
						if left < maxleft[lhs,right]: maxleft[lhs,right] = left
						if right < minright[lhs,left]: minright[lhs,left] = right
						if right > maxright[lhs,left]: maxright[lhs,left] = right
			# unary
			for rule in grammar.unary:
				lhs = rule.lhs
				ins = inside[left, right, rule.rhs1]
				if ins == 0.0: continue
				inside[left, right, lhs] += ins * rule.prob
				if left > minleft[lhs, right]: minleft[lhs, right] = left
				if left < maxleft[lhs, right]: maxleft[lhs, right] = left
				if right < minright[lhs, left]: minright[lhs, left] = right
				if right > maxright[lhs, left]: maxright[lhs, left] = right
	print >>stderr, ''
	return minleft, maxleft, minright, maxright

def outsidescores(Grammar grammar, long start, short lensent,
	np.ndarray[np.double_t, ndim=3] inside,
	np.ndarray[np.double_t, ndim=3] outside,
	np.ndarray[np.int16_t, ndim=2] minleft,
	np.ndarray[np.int16_t, ndim=2] maxleft,
	np.ndarray[np.int16_t, ndim=2] minright,
	np.ndarray[np.int16_t, ndim=2] maxright):
	cdef short left, right, mid, span
	cdef short narrowl, narrowr, minmid, maxmid
	cdef long numsymbols = len(grammar.toid), lhs, rhs1, rhs2
	cdef double ls, rs, os
	cdef bint foundbetter = False
	cdef Rule rule
	cdef Terminal terminal
	cdef unicode word
	assert not grammar.logprob, "Grammar must not have log probabilities."
	outside[0, lensent, start] = 1.0
	print >>stderr, "outside",
	for span in range(lensent, 0, -1):
		print >>stderr, span,
		stdout.flush()
		for left in range(1 + lensent - span):
			if left == lensent: continue
			right = left + span
			# unary
			for rule in grammar.unary:
				os = outside[left, right, rule.lhs]
				if os == 0.0: continue
				outside[left, right, rule.rhs1] += os * rule.prob
			# binary
			for lhs, rules in enumerate(grammar.binary):
				for rule in rules:
					os = outside[left, right, lhs]
					if os == 0.0: continue
					narrowr = minright[rule.rhs1, left]
					narrowl = minleft[rule.rhs2, right]
					if narrowr >= right or narrowl < narrowr: continue
					minmid = smax(narrowr, maxleft[rule.rhs2, right])
					maxmid = smin(maxright[rule.rhs1, left], narrowl) + 1
					for split in range(minmid, maxmid):
						ls = inside[left, split, rule.rhs1]
						if ls == 0.0: continue
						rs = inside[split, right, rule.rhs2]
						if rs == 0.0: continue
						outside[left, split, rule.rhs1] += rule.prob * rs * os
						outside[split, right, rule.rhs2] += rule.prob * ls * os
	print >>stderr, ''
	return outside

def dopparseprob(tree, Grammar grammar, dict mapping, lexchart):
	""" Given an NLTK tree, compute the exact DOP parse probability given
	a DOP reduction.

	This follows up on a suggestion made by Goodman (2003, p. 20)
	of calculating DOP probabilities of given parse trees, although I'm not
	sure it has complexity O(nP) as he suggests (with n as number of nodes in input,
	and P as max number of rules consistent with a node in the input).
	Furthermore, the idea of sampling trees "long enough" until we have the MPP
	is no faster than sampling without applying this procedure, because there
	is no way to determine that some probability p is the maximal probability,
	except in the unlikely case that p > 0.5. Hence, this method is mostly
	useful in a reranking framework where it is known in advance that a small
	set of trees is of interest.

	expects a mapping which gives a list of consistent rules from the reduction
	(as Rule objects) given a rule as key (as a tuple of strings); e.g. ('NP',
	'DT', 'NN') -> [Rule(...), Rule(...), ...]

	NB: this algorithm could also be used to determine the probability of
	derivations, but then the input would have to distinguish whether nodes are
	internal nodes of fragments, or whether they join two fragments. """
	neginf = float('-inf')
	cdef dict chart = {}	#chart[left, right][label]
	cdef tuple a, b, c
	cdef Rule rule
	assert grammar.logprob, "Grammar should have log probabilities."
	# log probabilities are not ideal here because we do lots of additions,
	# but the probabilities are very small. a possible alternative is to scale
	# them somehow.

	# add all possible POS tags
	chart.update(lexchart)
	for n, word in enumerate(tree.leaves()):
		# replace leaves with indices so as to easily find spans
		tree[tree.leaf_treeposition(n)] = n

	# do post-order traversal (bottom-up)
	for node in list(tree.subtrees())[::-1]:
		if not isinstance(node[0], Tree): continue
		prod = (node.node,) + tuple(a.node for a in node)
		left = min(node.leaves())
		right = max(node.leaves()) + 1
		if len(node) == 1: #unary node
			for rule in mapping[prod]:
				b = (rule.rhs1, left, right)
				if b in chart:
					a = (rule.lhs, left, right)
					if a in chart:
						chart[a] = logprobadd(chart[a], -rule.prob + chart[b])
					else:
						chart[a] = (-rule.prob + chart[b])
		elif len(node) == 2: #binary node
			split = min(node[1].leaves())
			for rule in mapping[prod]:
				b = (rule.rhs1, left, split)
				c = (rule.rhs2, split, right)
				if b in chart and c in chart:
					a = (rule.lhs, left, right)
					if a in chart:
						chart[a] = logprobadd(chart[a],
							(-rule.prob + chart[b] + chart[c]))
					else:
						chart[a] = -rule.prob + chart[b] + chart[c]
		else: raise ValueError("expected binary tree.")
	return chart.get((grammar.toid[tree.node], 0, len(tree.leaves())), neginf)

cdef inline short smax(short a, short b): return a if a >= b else b
cdef inline short smin(short a, short b): return a if a <= b else b

def readbitpargrammar(rules, lexiconfile, unknownwords, logprob=True, freqs=True):
	""" Reads grammars in bitpar's format; this one is actually less fussy
	about the difference between tabs and spaces. Does require a binarized
	grammar."""
	nonterminals = set(); unary = []; binary = []; lexical = []
	lexicon = set()
	print >>stderr, "reading the grammar...",
	stdout.flush()
	for a in open(rules):
		rule = a.rstrip().split()
		if len(rule) > 4: raise ValueError("rule is not binarized: %s" % rule)
		if len(rule) < 3:
			raise ValueError("malformed rule: %s while reading %s" % (
				rule, rules))
		nonterminals.update(rule[1:])
	print >>stderr, len(nonterminals), "non-terminals"
	symbols = ["Epsilon"] + sorted(nonterminals)
	toid = dict((a, n) for n,a in enumerate(symbols))
	tolabel = dict(enumerate(symbols))
	unarybyrhs = [[] for a in symbols]
	binary = [[] for a in symbols]
	fd = defaultdict(float)
	for a in open(rules):
		rule = a.rstrip().split()
		freq = float(rule[0])
		lhs = toid[rule[1]]
		rhs1 = toid[rule[2]]
		if len(rule) == 3:
			rule = Rule(lhs, rhs1, 0, freq)
			unary.append(rule)
			unarybyrhs[rhs1].append(rule)
		elif len(rule) == 4:
			rhs2 = toid[rule[3]]
			binary[lhs].append(Rule(lhs, rhs1, rhs2, freq))
		fd[lhs] += freq
	print >>stderr, "reading the lexicon...",
	stdout.flush()
	for a in codecs.open(lexiconfile, encoding=encoding):
		rule = a.split()
		word = rule[0]
		lexicon.add(word)
		for tag, freq in zip(rule[1::2], rule[2::2]):
			if tag not in nonterminals:
				raise ValueError("unknown symbol %s in lexicon" % tag)
			lhs = toid[tag]
			freq = float(freq)
			lexical.append(Terminal(lhs, word, freq))
			fd[lhs] += freq
	if unknownwords:
		for a in open(unknownwords):
			tag, freq = a.split()
			lhs = toid[tag]
			freq = float(freq)
			lexical.append(Terminal(lhs, u"", float(freq)))
			fd[lhs] += freq
	print >>stderr, "%d words\nparameter estimation..." % (len(lexicon)),
	stdout.flush()
	if freqs:
		# turn rule frequencies into relative frequencies
		if logprob:
			for prods in [unary, lexical] + binary:
				for a in prods: a.prob = -log(a.prob / fd[a.lhs])
		else:
			for prods in [unary, lexical] + binary:
				for a in prods: a.prob /= fd[a.lhs]
	for a in unary:
		#prevent cycles
		if a.prob == 0.0: a.prob = 0.0001
	print >>stderr, "finished"
	return Grammar(lexical, unary, unarybyrhs, binary, tolabel, toid,
			lexicon, logprob)

def reestimate(Grammar coarse, Grammar fine):
	""" Modify probabilities of coarse grammar such that they are the sum
	of probabilities of rules in the fine grammar that map to the same
	nonterminals.  """
	cdef Rule rule
	cdef dict lhs = <dict>defaultdict(list)
	cdef dict rhs = <dict>defaultdict(list)
	cdef dict mapping = dict((b, coarse.toid[a.rsplit("@", 1)[0]])
							for a, b in fine.toid.iteritems())
	assert fine.logprob
	for rule in fine.unary:
		lhs[mapping[rule.lhs]].append(rule.prob)
		rhs[mapping[rule.lhs], mapping[rule.rhs1]].append(rule.prob)
	for rules in fine.binary:
		for rule in rules:
			lhs[mapping[rule.lhs]].append(rule.prob)
			rhs[mapping[rule.lhs], mapping[rule.rhs1],
				mapping[rule.rhs2]].append(rule.prob)
	for rule in coarse.unary:
		rule.prob = logsum(rhs[rule.lhs, rule.rhs1]) - logsum(lhs[rule.lhs])
		if not coarse.logprob: rule.prob = exp(-rule.prob)
	for rules in coarse.binary:
		for rule in rules:
			rule.prob = (logsum(rhs[rule.lhs, rule.rhs1, rule.rhs2])
						- logsum(lhs[rule.lhs]))
			if not coarse.logprob: rule.prob = exp(-rule.prob)

cdef inline double logsum(list logprobs) except 1.0:
	# Adding probabilities in log space
	# http://blog.smola.org/post/987977550/log-probabilities-semirings-and-floating-point-numbers
	# https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations
	#NB: this version deals with and returns negative logprobs.
	maxprob = min(logprobs)
	return maxprob - log(sum([exp(maxprob - prob) for prob in logprobs]))

def getgrammarmapping(coarse, fine):
	""" producing a mapping of coarse rules to sets of fine rules;
	e.g. mapping["S", "NP", "VP"] == set(Rule(...), Rule(...), ...) """
	mapping = {}
	for rule in coarse.unary:
		mapping[coarse.tolabel[rule.lhs], coarse.tolabel[rule.rhs1]] = []
	for rules in coarse.binary:
		for rule in rules:
			mapping[coarse.tolabel[rule.lhs],
				coarse.tolabel[rule.rhs1], coarse.tolabel[rule.rhs2]] = []
	for rule in fine.unary:
		mapping[fine.tolabel[rule.lhs].rsplit("@", 1)[0],
			fine.tolabel[rule.rhs1].rsplit("@", 1)[0]].append(rule)
	for rules in fine.binary:
		for rule in rules:
			mapping[fine.tolabel[rule.lhs].rsplit("@", 1)[0],
				fine.tolabel[rule.rhs1].rsplit("@", 1)[0],
				fine.tolabel[rule.rhs2].rsplit("@", 1)[0]].append(rule)
	return mapping

def cachingdopparseprob(tree, Grammar grammar, dict mapping, dict cache):
	""" Given an NLTK tree, compute the DOP parse probability given a
	DOP reduction.
	Experimental version which maintains a memoization table. """
	#from bigfloat import BigFloat, setcontext, quadruple_precision
	#setcontext(quadruple_precision)
	#chart = defaultdict(lambda: BigFloat(0))
	#chart = defaultdict(float)
	neginf = float('-inf')
	cdef dict chart = <dict>defaultdict(lambda: neginf)
	cdef Rule rule
	cdef Terminal terminal

	for n, word in enumerate(tree.leaves()):
		for terminal in grammar.lexical:
			if terminal.word == word:
				chart[terminal.lhs, n, n+1] = logprobadd(
					chart[terminal.lhs, n, n+1], -terminal.prob)
				#chart[terminal.lhs, n, n+1] += terminal.prob
		# replace leaves with indices so as to easily find spans
		tree[tree.leaf_treeposition(n)] = n
	tree = tree.freeze()

	for node in list(tree.subtrees())[::-1]:
		if not isinstance(node[0], Tree): continue
		if node in cache:
			chart.update(cache[node])
			continue
		prod = (node.node,) + tuple(a.node for a in node)
		left = min(node.leaves())
		right = max(node.leaves()) + 1
		if len(node) == 2:
			split = min(node[1].leaves())
			for rule in mapping[prod]:
				chart[rule.lhs, left, right] = logprobadd(
					chart[rule.lhs, left, right],
					(-rule.prob
					+ chart[rule.rhs1, left, split]
					+ chart[rule.rhs2, split, right]))
				#chart[rule.lhs, left, right] += (rule.prob
				#	* chart[rule.rhs1, left, split]
				#	* chart[rule.rhs2, split, right])
		elif len(node) == 1:
			for rule in mapping[prod]:
				chart[rule.lhs, left, right] = logprobadd(
					chart[rule.lhs, left, right],
					-rule.prob + chart[rule.rhs1, left, right])
				#chart[rule.lhs, left, right] += (rule.prob
				#		* chart[rule.rhs1, left, right])
		else: raise ValueError("expected binary tree.")
		cache[node] = chart.items()
		#print prod[0], left, right, chart[grammar.toid[prod[0]], left, right]
	return chart[grammar.toid[tree.node], 0, len(tree.leaves())]

def doplexprobs(tree, Grammar grammar):
	neginf = float('-inf')
	cdef dict chart = <dict>defaultdict(lambda: neginf)
	cdef Terminal terminal

	for n, word in enumerate(tree.leaves()):
		for terminal in grammar.lexical:
			if terminal.word == word:
				chart[terminal.lhs, n, n+1] = logprobadd(
					chart[terminal.lhs, n, n+1], -terminal.prob)
	return chart

cdef double log1e200 = log(1e200)
cdef inline logprobadd(double x, double y):
	""" add two log probabilities in log space.
	i.e., logprobadd(log(a), log(b)) == log(a + b) """
	if isinf(x): return y
	elif isinf(y): return x
	# If one value is much smaller than the other, keep the larger value.
	elif x < (y - log1e200): return y
	elif y < (x - log1e200): return x
	diff = y - x
	assert not isinf(diff)
	if isinf(exp(diff)):	# difference is too large
		return x if x > y else y
	# otherwise return the sum.
	return x + log(1.0 + exp(diff))

def pprint_matrix(matrix, sent, tolabel):
	for lhs in tolabel:
		for left in range(len(sent)):
			for right in range(left + 1, len(sent) + 1):
				if matrix[left, right, lhs]:
					print "%s[%d:%d] = %f" % (
						tolabel[lhs], left, right, matrix[left, right, lhs])


def pprint_chart(chart, sent, tolabel):
	""" `pretty print' a chart. """
	cdef Edge edge
	print "chart:"
	for left in range(len(sent)):
		for right in range(left + 1, len(sent) + 1):
			for label, edges in chart[left][right].items():
				print "%s[%d:%d] =>" % (tolabel[label], left, right),
				if not edges: print "[]"
				else: print
				for edge in edges:
					print "%g\t%g" % (exp(-edge.inside), exp(-edge.rule.prob)),
					if edge.rule.rhs1:
						print "\t%s[%d:%d]" % (tolabel[edge.rule.rhs1],
													left, edge.split),
					else:
						print "\t", repr(sent[left]),
					if edge.rule.rhs2:
						print "\t%s[%d:%d]" % (tolabel[edge.rule.rhs2],
												edge.split, right),
					print
