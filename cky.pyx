""" Probabilistic CKY parser for PCFGs """
import sys, codecs
from collections import defaultdict
from math import log, exp
import numpy as np
from agenda cimport Agenda
encoding = "utf-8"

def parse(list sent, Grammar grammar, whitelist):
	""" A CKY parser modeled after Bodenstab's `fast grammar loop.'
		and the Stanford parser. """
	cdef short left, right, mid, span, lensent = len(sent)
	cdef short narrowr, narrowl, widel, wider, minmid, maxmid
	cdef long numsymbols = len(grammar.toid), lhs
	cdef double oldscore, prob
	cdef bint foundbetter = False
	cdef Rule rule
	cdef Terminal terminal
	cdef dict chart = {}
	# the viterbi chart is initially filled with infinite log probabilities,
	# cells which are to be blocked contain NaN.
	cdef np.ndarray[np.double_t, ndim=3] viterbi = np.array([np.inf],
		dtype='d').repeat(lensent * (lensent+1) * numsymbols).reshape(
		(numsymbols, lensent, (lensent+1))) if whitelist is None else whitelist
	cdef np.ndarray[np.int16_t, ndim=2] minsplitleft = np.array([-1],
		dtype='int16').repeat(numsymbols * (lensent + 1)
		).reshape(numsymbols, lensent + 1)
	cdef np.ndarray[np.int16_t, ndim=2] maxsplitleft = np.array([lensent+1],
		dtype='int16').repeat(numsymbols * (lensent + 1)).reshape(
		numsymbols, lensent + 1)
	cdef np.ndarray[np.int16_t, ndim=2] minsplitright = np.array([lensent + 1],
		dtype='int16').repeat(numsymbols * (lensent + 1)
		).reshape(numsymbols, lensent + 1)
	cdef np.ndarray[np.int16_t, ndim=2] maxsplitright = np.array([-1],
		dtype='int16').repeat(numsymbols * (lensent + 1)).reshape(
		numsymbols, lensent + 1)
	# assign POS tags
	for left in range(lensent):
		right = left + 1
		for terminal in <list>grammar.lexical:
			if terminal.word == (sent[left] if sent[left] in grammar.lexicon else u""):
				viterbi[terminal.lhs, left, right] = terminal.prob
				chart[new_ChartItem(terminal.lhs, left, right)] = [
					new_Edge(terminal.prob, terminal, right)]
				# update filter
				if left > minsplitleft[terminal.lhs, right]:
					minsplitleft[terminal.lhs, right] = left
				if left < maxsplitleft[terminal.lhs, right]:
					maxsplitleft[terminal.lhs, right] = left
				if right < minsplitright[terminal.lhs, left]:
					minsplitright[terminal.lhs, left] = right
				if right > maxsplitright[terminal.lhs, left]:
					maxsplitright[terminal.lhs, left] = right
		# unary rules on POS tags 
		for rule in <list>grammar.unary:
			if (not isnan(viterbi[rule.lhs, left, right])
				and isfinite(viterbi[rule.rhs1, left, right])):
				prob = rule.prob + viterbi[rule.rhs1, left, right]
				if isfinite(viterbi[rule.lhs, left, right]):
					(<list>chart[new_ChartItem(rule.lhs, left,
						right)]).append(new_Edge(prob, rule, right))
				else:
					chart[new_ChartItem(rule.lhs, left, right)] = [
						new_Edge(prob, rule, right)]
				if (prob < viterbi[rule.lhs, left, right]):
					viterbi[rule.lhs, left, right] = prob
					# update filter
					if left > minsplitleft[rule.lhs, right]:
						minsplitleft[rule.lhs, right] = left
					if left < maxsplitleft[rule.lhs, right]:
						maxsplitleft[rule.lhs, right] = left
					if right < minsplitright[rule.lhs, left]:
						minsplitright[rule.lhs, left] = right
					if right > maxsplitright[rule.lhs, left]:
						maxsplitright[rule.lhs, left] = right
	for span in range(1, lensent + 1):
		print span,
		sys.stdout.flush()
		# constituents from left to right
		for left in range(0, lensent - span + 1):
			right = left + span
			# binary rules 
			for lhs, rules in enumerate(<list>grammar.binary):
				if isnan(viterbi[lhs, left, right]): continue
				for rule in <list>rules:
					narrowr = minsplitright[rule.rhs1, left]
					if narrowr >= right: continue
					narrowl = minsplitleft[rule.rhs2, right]
					if narrowl < narrowr: continue
					widel = maxsplitleft[rule.rhs2, right]
					minmid = smax(narrowr, widel)
					wider = maxsplitright[rule.rhs1, left]
					maxmid = smin(wider, narrowl) + 1
					oldscore = viterbi[lhs, left, right]
					foundbetter = False
					for mid in range(minmid, maxmid):
						if (isfinite(viterbi[rule.rhs1, left, mid])
							and isfinite(viterbi[rule.rhs2, mid, right])):
							prob = (rule.prob + viterbi[rule.rhs1, left, mid]
									+ viterbi[rule.rhs2, mid, right])
							if isfinite(viterbi[lhs, left, right]):
								(<list>chart[new_ChartItem(lhs, left, right)
									]).append(new_Edge(prob, rule, mid))
							else:
								chart[new_ChartItem(lhs, left, right)] = [
									new_Edge(prob, rule, mid)]
							if prob < viterbi[lhs, left, right]:
								foundbetter = True
								viterbi[lhs, left, right] = prob
					# update filter
					if foundbetter and isinf(oldscore):
						if left > minsplitleft[lhs, right]:
							minsplitleft[lhs, right] = left
						if left < maxsplitleft[lhs, right]:
							maxsplitleft[lhs, right] = left
						if right < minsplitright[lhs, left]:
							minsplitright[lhs, left] = right
						if right > maxsplitright[lhs, left]:
							maxsplitright[lhs, left] = right

			# unary rules
			if span == 1: continue
			for rule in <list>grammar.unary:
				if (not isnan(viterbi[rule.lhs, left, right])
					and isfinite(viterbi[rule.rhs1, left, right])):
					prob = rule.prob + viterbi[rule.rhs1, left, right]
					if isfinite(viterbi[rule.lhs, left, right]):
						(<list>chart[new_ChartItem(rule.lhs, left,
							right)]).append(new_Edge(prob, rule, right))
					else:
						chart[new_ChartItem(rule.lhs, left, right)] = [
							new_Edge(prob, rule, right)]
					if prob < viterbi[rule.lhs, left, right]:
						viterbi[rule.lhs, left, right] = prob
						# update filter
						if left > minsplitleft[rule.lhs, right]:
							minsplitleft[rule.lhs, right] = left
						if left < maxsplitleft[rule.lhs, right]:
							maxsplitleft[rule.lhs, right] = left
						if right < minsplitright[rule.lhs, left]:
							minsplitright[rule.lhs, left] = right
						if right > maxsplitright[rule.lhs, left]:
							maxsplitright[rule.lhs, left] = right
	print
	return chart, viterbi

def insidescores1(chart, start, grammar, minsplitleft, maxsplitleft, minsplitright, maxsplitright):
	lensent = start.right
	inside = np.array([np.inf], dtype='d').repeat(
		len(grammar.toid) * lensent * (lensent+1)).reshape(
		(len(grammar.toid), lensent, (lensent+1)))
	for span in range(lensent, -1, -1):
		for left in range(lensent - span):
			right = left + span
			# unary
			for rhs in range(len(grammar.toid)):
				ins = inside[rhs, left, right]
				if isinf(ins): continue
				for rule in grammar.unary:
					if rule.rhs1 != rhs: continue
					prob = ins + rule.prob
					inside[rule.lhs, left, right] += exp(-prob)
			# binary
			for lhs, rules in enumerate(grammar.binary):
				min1 = minsplitright[lhs, left]
				if right < min1: continue
				for rule in rules:
					max1 = minsplitleft[rule.rhs2][right]
					if max1 < min1: continue
					if max1 - min1 > 2:
						min2 = maxsplitleft[rule.rhs2, right]
						if min2 > min1: min1 = min2
						if max1 < min1: continue
						max2 = maxsplitright[rule.rhs1, left]
						if max2 < max1: max1 = max2
						if max1 < min1: continue
					for split in range(min1, max1 + 1):
						ls = chart[rule.rhs1, left, split]
						if isinf(ls): continue
						rs = chart[rule.rhs2, split, right]
						if isinf(rs): continue
						tot = rule.prob + ls + rs
						inside[rule.lhs1, left, right] += exp(-tot)

def outsidescores1(inside, start, grammar, minsplitleft, maxsplitleft, minsplitright, maxsplitright):
	# expects grammar.binaryleft, grammar.binaryright
	lensent = start.right
	outside = np.array([np.inf], dtype='d').repeat(
		len(grammar.toid) * lensent * (lensent+1)).reshape(
		(len(grammar.toid), lensent, (lensent+1)))
	outside[start.label, 0, lensent] = 0.0
	for span in range(lensent, -1, -1):
		for left in range(lensent - span):
			right = left + span
			# unary
			for lhs in range(len(grammar.toid)):
				os = outside[lhs, left, right]
				if isinf(os): continue
				for rule in grammar.unary:
					if rule.lhs != lhs: continue
					prob = os + rule.prob
					if (prob < outside[rule.rhs1, left, right]
						and isfinite(inside[rule.rhs1, left, right])):
						outside[rule.rhs1, left, right] = prob
			# binary-left
			for lhs, rules in enumerate(grammar.binaryleft):
				min1 = minsplitright[lhs, left]
				if right < min1: continue
				for rule in rules:
					os = outside[rule.lhs, left, right]
					if isinf(os): continue
					max1 = minsplitleft[rule.rhs2][right]
					if max1 < min1: continue
					if max1 - min1 > 2:
						min2 = maxsplitleft[rule.rhs2, right]
						if min2 > min1: min1 = min2
						if max1 < min1: continue
						max2 = maxsplitright[rule.rhs1, left]
						if max2 < max1: max1 = max2
						if max1 < min1: continue
					for split in range(min1, max1 + 1):
						ls = inside[rule.rhs1, left, split]
						if isinf(ls): continue
						rs = inside[rule.rhs2, split, right]
						if isinf(rs): continue
						outside[rule.rhs1, left, split] += rule.prob + rs + os
						outside[rule.rhs2, split, right] += rule.prob + ls + os
			# binary-right
			for lhs, rules in enumerate(grammar.binaryright):
				max1 = minsplitleft[lhs, right]
				if max1 < left: continue
				for rule in rules:
					os = outside[lhs, left, right]
					if isinf(os): continue
					min1 = minsplitright[rule.rhs1][left]
					if max1 < min1: continue
					if max1 - min1 > 2:
						min2 = maxsplitleft[rule.rhs2, right]
						if min2 > min1: min1 = min2
						if max1 < min1: continue
						max2 = maxsplitright[rule.rhs1, left]
						if max2 < max1: max1 = max2
						if max1 < min1: continue
					for split in range(min1, max1 + 1):
						ls = inside[rule.rhs1, left, split]
						if isinf(ls): continue
						rs = inside[rule.rhs2, split, right]
						if isinf(rs): continue
						outside[rule.rhs1, left, split] += rule.prob + rs + os
						outside[rule.rhs2, split, right] += rule.prob + ls + os
	return outside

# to avoid overhead of __init__ and __cinit__ constructors
# belongs in containers but putting it here gives
# a better chance of successful inlining
cdef inline ChartItem new_ChartItem(unsigned int label, short left, short right):
	cdef ChartItem item = ChartItem.__new__(ChartItem)
	item.label = label; item.left = left; item.right = right
	return item

cdef inline Edge new_Edge(double inside, Rule rule, short split):
	cdef Edge edge = Edge.__new__(Edge)
	edge.inside = inside; edge.rule = rule; edge.split = split
	return edge

cdef inline short smax(short a, short b): return a if a >= b else b
cdef inline short smin(short a, short b): return a if a <= b else b

def readbitpargrammar(rules, lexiconfile, unknownwords, normalize=True):
	""" Reads grammars in bitpar's format; this one is actually less fussy
	about the difference between tabs and spaces. Does require a binarized
	grammar."""
	nonterminals = set(); unary = []; binary = []; lexical = []
	lexicon = set()
	print "reading the grammar...",
	sys.stdout.flush()
	for a in open(rules):
		rule = a.rstrip().split()
		if len(rule) > 4: raise ValueError("rule is not binarized: %s" % rule)
		if len(rule) < 3:
			raise ValueError("malformed rule: %s while reading %s" % (
				rule, rules))
		nonterminals.update(rule[1:])
	print len(nonterminals), "non-terminals"
	symbols = ["Epsilon"] + sorted(nonterminals)
	toid = dict((a, n) for n,a in enumerate(symbols))
	tolabel = dict(enumerate(symbols))
	binary = [[] for a in symbols]
	fd = defaultdict(float)
	for a in open(rules):
		rule = a.rstrip().split()
		freq = float(rule[0])
		lhs = toid[rule[1]]
		rhs1 = toid[rule[2]]
		if len(rule) == 3:
			unary.append(Rule(lhs, rhs1, 0, freq))
		elif len(rule) == 4:
			rhs2 = toid[rule[3]]
			binary[lhs].append(Rule(lhs, rhs1, rhs2, freq))
		fd[lhs] += freq
	print "reading the lexicon...",
	sys.stdout.flush()
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
	print "%d words\nparameter estimation..." % (len(lexicon)),
	sys.stdout.flush()
	if normalize:
		# turn rule frequencies into relative frequencies (as log probs)
		for prods in [unary, lexical] + binary:
			for a in prods: a.prob = -log(a.prob / fd[a.lhs])
	for a in unary:
		#prevent cycles
		if a.prob == 0.0: a.prob = 0.0001
	print "finished"
	return Grammar(lexical, unary, binary, tolabel, toid, lexicon)

def reestimate(Grammar coarse, Grammar fine):
	""" Modify probabilities of coarse grammar such that the are the sum
	of probabilities of rules in the fine grammar that map to the same
	nonterminals.  """
	cdef dict lhs = <dict>defaultdict(list)
	cdef dict rhs = <dict>defaultdict(list)
	cdef dict mapping = dict((b, coarse.toid[a.split("@")[0]])
							for a,b in fine.toid.iteritems())
	for rule in fine.unary:
		lhs[mapping[rule.lhs]].append(rule.prob)
		rhs[mapping[rule.lhs], mapping[rule.rhs1]].append(rule.prob)
	for rules in fine.binary:
		for rule in rules:
			lhs[mapping[rule.lhs]].append(rule.prob)
			rhs[mapping[rule.lhs], mapping[rule.rhs1],
				mapping[rule.rhs2]].append(rule.prob)
	for rule in coarse.unary:
		rule.prob = (logsum(rhs[mapping[rule.lhs], mapping[rule.rhs1]])
					- logsum(lhs[mapping[rule.lhs]]))
	for rules in coarse.binary:
		for rule in rules:
			rule.prob = (logsum(rhs[mapping[rule.lhs],
						mapping[rule.rhs1], mapping[rule.rhs2]])
						- logsum(lhs[mapping[rule.lhs]]))

def logsum(list logprobs):
	# Adding probabilities in log space
	# http://blog.smola.org/post/987977550/log-probabilities-semirings-and-floating-point-numbers
	# https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations
	#NB: this version deals with and returns negative logprobs.
	maxprob = min(logprobs)
	return maxprob - log(sum([exp(maxprob - prob) for prob in logprobs]))

def pprint_chart(chart, sent, tolabel):
	""" `pretty print' a chart. """
	cdef ChartItem a
	cdef Edge edge
	print "chart:"
	for a in sorted(chart, key=lambda a: a.left):
		if not chart[a]: continue
		print "%s[%d:%d] =>" % (tolabel[a.label], a.left, a.right)
		for edge in chart[a]:
			print "%g\t%g" % (exp(-edge.inside), exp(-edge.rule.prob)),
			if edge.rule.rhs1:
				print "\t%s[%d:%d]" % (tolabel[edge.rule.rhs1],
											a.left, edge.split),
			else:
				print "\t", repr(sent[a.left]),
			if edge.rule.rhs2:
				print "\t%s[%d:%d]" % (tolabel[edge.rule.rhs2], edge.split, a.right),
			print
