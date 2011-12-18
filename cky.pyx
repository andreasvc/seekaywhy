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
	cdef unicode word
	# the viterbi chart is initially filled with infinite log probabilities,
	# cells which are to be blocked contain NaN.
	cdef np.ndarray[np.double_t, ndim=3] viterbi
	# matrices for the filter which gives minima and maxima for splits
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
	if whitelist is None:
		viterbi = np.array([np.inf],
		dtype='d').repeat(lensent * (lensent+1) * numsymbols).reshape(
		(numsymbols, lensent, (lensent+1)))
	else: viterbi = whitelist
	# assign POS tags
	for left in range(lensent):
		right = left + 1
		for terminal in <list>grammar.lexical:
			word = unicode(sent[left]) if sent[left] in grammar.lexicon else u""
			if terminal.word == word:
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

def doinsideoutside(list sent, Grammar grammar, inside, outside):
	lensent = len(sent); numsymbols = len(grammar.toid)
	start = ChartItem(grammar.toid["TOP"], 0, lensent)
	if inside == None:
		inside = np.array([0.0], dtype='d'
			).repeat(lensent * (lensent+1) * numsymbols
			).reshape((numsymbols, lensent, (lensent+1)))
	else:
		inside[:,:len(sent),:len(sent)+1] = 0.0
	if outside == None:
		outside = np.array([0.0], dtype='d'
			).repeat(lensent * (lensent+1) * numsymbols
			).reshape((numsymbols, lensent, (lensent+1)))
	else:
		outside[:,:len(sent),:len(sent)+1] = 0.0
	isl, asl, isr, asr = insidescores(sent, grammar, inside)
	outside = outsidescores(inside, start, grammar, outside, isl, asl, isr, asr)
	return inside, outside

cdef insidescores(list sent, Grammar grammar,
	np.ndarray[np.double_t, ndim=3] inside):
	""" Grammar must not have log probabilities. """
	cdef short left, right, mid, span, lensent = len(sent)
	cdef short narrowr, narrowl, widel, wider, minmid, maxmid
	cdef long numsymbols = len(grammar.toid), lhs
	cdef double oldscore, prob, ls, rs, ins
	cdef bint foundbetter = False
	cdef Rule rule
	cdef Terminal terminal
	cdef unicode word
	# matrices for the filter which give minima and maxima for splits
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
	inside[:,:lensent,:lensent+1] = 0.0
	print "inside",
	# assign POS tags
	for left in range(lensent):
		right = left + 1
		word = unicode(sent[left]) if sent[left] in grammar.lexicon else u""
		for terminal in <list>grammar.lexical:
			if terminal.word == word:
				inside[terminal.lhs, left, right] = terminal.prob
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
			if inside[rule.rhs1, left, right] != 0.0:
				prob = rule.prob * inside[rule.rhs1, left, right]
				inside[rule.lhs, left, right] += prob
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
			# binary
			for lhs, rules in enumerate(grammar.binary):
				for rule in rules:
					narrowr = minsplitright[rule.rhs1, left]
					if narrowr >= right: continue
					narrowl = minsplitleft[rule.rhs2, right]
					if narrowl < narrowr: continue
					widel = maxsplitleft[rule.rhs2, right]
					minmid = smax(narrowr, widel)
					wider = maxsplitright[rule.rhs1, left]
					maxmid = smin(wider, narrowl) + 1
					#oldscore = inside[lhs, left, right]
					foundbetter = False
					for split in range(minmid, maxmid):
						ls = inside[rule.rhs1, left, split]
						if ls == 0.0: continue
						rs = inside[rule.rhs2, split, right]
						if rs == 0.0: continue
						foundbetter = True
						inside[rule.lhs, left, right] += rule.prob * ls * rs
					if foundbetter: #and oldscore == 0.0:
						if left > minsplitleft[rule.lhs, right]:
							minsplitleft[rule.lhs, right] = left
						if left < maxsplitleft[rule.lhs, right]:
							maxsplitleft[rule.lhs, right] = left
						if right < minsplitright[rule.lhs, left]:
							minsplitright[rule.lhs, left] = right
						if right > maxsplitright[rule.lhs, left]:
							maxsplitright[rule.lhs, left] = right
			# unary
			for rule in grammar.unary:
				ins = inside[rule.rhs1, left, right]
				if ins == 0.0: continue
				inside[rule.lhs, left, right] += ins * rule.prob
				if left > minsplitleft[rule.lhs, right]:
					minsplitleft[rule.lhs, right] = left
				if left < maxsplitleft[rule.lhs, right]:
					maxsplitleft[rule.lhs, right] = left
				if right < minsplitright[rule.lhs, left]:
					minsplitright[rule.lhs, left] = right
				if right > maxsplitright[rule.lhs, left]:
					maxsplitright[rule.lhs, left] = right
	print
	return maxsplitleft, minsplitleft, maxsplitright, minsplitright

cdef outsidescores(np.ndarray[np.double_t, ndim=3] inside,
	ChartItem start,
	Grammar grammar,
	np.ndarray[np.double_t, ndim=3] outside,
	np.ndarray[np.int16_t, ndim=2] minsplitleft,
	np.ndarray[np.int16_t, ndim=2] maxsplitleft,
	np.ndarray[np.int16_t, ndim=2] minsplitright,
	np.ndarray[np.int16_t, ndim=2] maxsplitright):
	cdef short left, right, mid, span, lensent = start.right
	cdef short narrowr, narrowl, widel, wider, minmid, maxmid
	cdef long numsymbols = len(grammar.toid), lhs, rhs1, rhs2
	cdef double oldscore, ls, rs, os
	cdef bint foundbetter = False
	cdef Rule rule
	cdef Terminal terminal
	cdef unicode word
	outside[start.label, 0, lensent] = 1.0
	print "outside",
	for span in range(lensent, 0, -1):
		print span,
		sys.stdout.flush()
		for left in range(1 + lensent - span):
			if left == lensent: continue
			right = left + span
			# unary
			for rule in grammar.unary:
				os = outside[rule.lhs, left, right]
				if os == 0.0: continue
				outside[rule.rhs1, left, right] += os * rule.prob
			# binary
			for lhs, rules in enumerate(grammar.binary):
				for rule in rules:
					os = outside[lhs, left, right]
					if os == 0.0: continue
					narrowr = minsplitright[rule.rhs1, left]
					#if narrowr >= right: continue
					narrowl = minsplitleft[rule.rhs2, right]
					#if narrowl < narrowr: continue
					widel = maxsplitleft[rule.rhs2, right]
					minmid = smax(narrowr, widel)
					wider = maxsplitright[rule.rhs1, left]
					maxmid = smin(wider, narrowl) + 1
					#for split in range(minmid, maxmid):
					for split in range(left + 1, right):
						ls = inside[rule.rhs1, left, split]
						if ls == 0.0: continue
						rs = inside[rule.rhs2, split, right]
						if rs == 0.0: continue
						outside[rule.rhs1, left, split] += rule.prob * rs * os
						outside[rule.rhs2, split, right] += rule.prob * ls * os
	print
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

def readbitpargrammar(rules, lexiconfile, unknownwords, logprob=True, normalize=True):
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

def pprint_matrix(matrix, sent, tolabel):
	for lhs in tolabel:
		for left in range(len(sent)):
			for right in range(left + 1, len(sent) + 1):
				if matrix[lhs, left, right]:
					print "%s[%d:%d] = %f" % (
						tolabel[lhs], left, right, matrix[lhs, left, right])

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
