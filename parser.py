import sys, time, re
from getopt import gnu_getopt
from operator import itemgetter
from math import exp, log, isinf
from heapq import nlargest
import numpy as np
from nltk import Tree
from cky import parse, parse_sparse, readbitpargrammar, doinsideoutside, \
	pprint_chart, pprint_matrix, dopparseprob, getgrammarmapping, \
	cachingdopparseprob, doplexprobs
from kbest import lazykbest
from containers import ChartItem
from coarsetofine import whitelistfromkbest, whitelistfromposteriors, \
	whitelistfromposteriors2
from disambiguation import marginalize
removeids = re.compile("@[0-9_]+")

usage = """usage: %s [options] rules lexicon [input [output]]
or: %s [--rerank] [options] coarserules coarselexicon finerules finelexicon \
[input [output]]

Grammars need to be binarized, and are in bitpar format.
When no file is given, output is written to standard output; when additionally
no input is given, it is raed from standard input.

Options:
-u file         handle unknown words using given open class tags list
-b k            return the k-best parses instead of just 1
--prob          print probabilities as well as parse trees
--threshold p   in coarse-to-fine mode, the posterior threshold (as log prob)
--rerank k      enable DOP reranking mode: find DOP parse probabilities
                for k coarse derivations.
""" % (sys.argv[0], sys.argv[0])
#TODO:
#--quiet         produce no output except parses
#--kbestctf k    use k-best coarse-to-fine; instead of posterior threshold,
#                use k-best derivations as threshold

def main():
	print "SeeKayWhy PCFG parser - Andreas van Cranenburgh"
	options = ("quiet", "prob", "rerank=", "threshold=", "kbestctf=")
	opts, args = gnu_getopt(sys.argv[1:], "u:b:", options)
	opts = dict(opts)
	unknownwords = opts.get("-u")
	k = opts.get("-b", 1)
	prob = "--prob" in opts
	if 2 <= len(args) <= 4:
		simple(args, unknownwords, k, prob)
	elif 4 <= len(args) <= 6 and "--rerank" in opts:
			rerank(args, unknownwords, k, prob, int(opts['--rerank']))
	elif 4 <= len(args) <= 6:
		if "--kbestctf" in opts:
			threshold = int(opts.get("--kbestctf"))
		else: threshold = float(opts.get("--threshold", -6.2)) #0.01
		ctf(args, unknownwords, k, prob, threshold, posterior=not opts.get("--kbestctf"))
	else: print usage

def simple(args, unknownwords, k, printprob):
	grammar = readbitpargrammar(args[0], args[1], unknownwords)
	input = sys.stdin; out = sys.stdout
	if len(args) >= 3: input = open(args[2])
	if len(args) == 4: out = open(args[3], "w")
	times = [time.clock()]
	for n, a in enumerate(input.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		for word in sent:
			assert word in grammar.lexicon or unknownwords, (
				"unknown word %r and no open class tags supplied" % word)
		print "parsing:", n, " ".join(sent),
		sys.stdout.flush()
		chart, viterbi = parse(sent, grammar, None)
		start = ChartItem(grammar.toid["TOP"], 0, len(sent))
		#pprint_chart(chart, sent, grammar.tolabel)
		if chart[0][len(sent)].get(grammar.toid["TOP"]):
			parsetrees = lazykbest(chart, start, k, grammar.tolabel, sent)
			assert len(parsetrees) == len(set(parsetrees))
			assert len(parsetrees) == len(set(tree for tree, prob in parsetrees))
			if printprob:
				out.writelines("vitprob=%.16g\n%s\n" % (exp(-prob), tree)
					for tree, prob in parsetrees)
			else: out.writelines("%s\n" % tree for tree, _ in parsetrees)
		else:
			out.write("(NP %s)\n" % "".join("(%s %s)" % (a,a) for a in sent))
			#out.write("No parse for \"%s\"\n" % " ".join(sent))
		#out.write("\n")
		out.flush()
		times.append(time.clock())
		print times[-1] - times[-2], "s"
	print "raw cpu time", time.clock() - times[0]
	times = [a - b for a, b in zip(times[1::2], times[::2])]
	print "average time per sentence", sum(times) / len(times)
	print "finished"
	out.close()

def ctf(args, unknownwords, k, printprob, threshold, posterior=True, mpd=False):
	m = 10000
	derivthreshold = 1000	# kbest derivations to prune with
	if posterior: threshold = exp(threshold)
	maxlen = 65
	unparsed = 0
	coarse = readbitpargrammar(args[0], args[1], unknownwords, logprob=not posterior)
	fine = readbitpargrammar(args[2], args[3], unknownwords, freqs=False)
	for a in fine.toid:
		assert a.rsplit("@", 1)[0] in coarse.toid, "%s not in coarse grammar" % a
	input = sys.stdin; out = sys.stdout
	if len(args) >= 5: input = open(args[4])
	if len(args) == 6: out = open(args[5], "w")
	times = [time.clock()]
	if posterior:
		inside = np.zeros((maxlen, maxlen + 1, len(coarse.toid)), dtype='d')
		outside = np.zeros_like(inside)
	else:
		coarsechart = np.empty((maxlen, maxlen + 1, len(coarse.toid)), dtype='d')
		coarsechart.fill(np.inf)
		finechart = np.empty_like(coarsechart)
		finechart.fill(np.NAN)
	mapping = nonterminalmapping(coarse, fine)

	for n, a in enumerate(input.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		if len(sent) > maxlen: continue
		for word in sent:
			assert unknownwords or (
				word in coarse.lexicon and word in fine.lexicon), (
				"unknown word and no open class tags supplied")
		print "parsing:", n, " ".join(sent)
		start = ChartItem(coarse.toid["TOP"], 0, len(sent))
		if posterior:
			inside, outside = doinsideoutside(sent, coarse, inside, outside)
			#print "inside"; pprint_matrix(inside, sent, coarse.tolabel)
			#print "outside"; pprint_matrix(outside, sent, coarse.tolabel)
		else:
			coarsechart[:len(sent), :len(sent)+1, :] = np.inf
			chart, coarsechart = parse(sent, coarse, coarsechart)
		if posterior: goalitem = inside[0, len(sent), coarse.toid["TOP"]]
		else: goalitem = chart[0][len(sent)].get(coarse.toid["TOP"])
		if goalitem:
			print "pruning ...",
			sys.stdout.flush()
			if posterior:
				finechart = whitelistfromposteriors2(inside, outside, start,
					coarse, fine, mapping, maxlen, threshold)
			else:
				whitelistfromkbest(chart, start, coarse, fine, threshold,
					finechart, maxlen)
			#chart, finechart = parse(sent, fine, finechart)
			chart = parse_sparse(sent, fine, finechart)
			start = ChartItem(fine.toid["TOP"], 0, len(sent))

			try:
				assert chart[0][len(sent)][fine.toid["TOP"]], (
				"sentence covered by coarse grammar could not be parsed "\
				"by fine grammar")
			except (AssertionError, KeyError):
				pprint_chart(chart, sent, fine.tolabel)
				raise
			parsetrees = marginalize(chart, start, fine.tolabel, sent, n=m, mpd=mpd)
			results = nlargest(k, parsetrees, key=parsetrees.get)
			# print k-best parsetrees
			if printprob:
				label = "derivprob" if mpd else "parseprob"
				out.writelines("%s=%.16g\n%s\n" % (label, parsetrees[tree], tree)
					for tree in results)
			else:
				out.writelines("%s\n" % tree for tree in results)
		else:
			unparsed += 1
			print "No parse"
			out.write("No parse for \"%s\"\n" % " ".join(sent))
		out.write("\n")
		times.append(time.clock())
		print times[-1] - times[-2], "s"
		out.flush()
	print "raw cpu time", time.clock() - times[0]
	times = [a - b for a, b in zip(times[1::2], times[::2])]
	print "average time per sentence", sum(times) / len(times)
	print "unparsed sentences:", unparsed
	print "finished"
	out.close()

def rerank(args, unknownwords, printk, printprob, k):
	maxlen = 999 #??
	unparsed = 0
	coarse = readbitpargrammar(args[0], args[1], unknownwords, logprob=True)
	fine = readbitpargrammar(args[2], args[3], unknownwords, freqs=False, logprob=True)
	for a in fine.toid:
		assert a.rsplit("@", 1)[0] in coarse.toid, "%s not in coarse grammar" % a
	if len(args) >= 5: input = open(args[4])
	else: input = sys.stdin
	if len(args) == 6: out = open(args[5], "w")
	else: out = sys.stdout
	times = [time.clock()]
	mapping = getgrammarmapping(coarse, fine)
	for n, a in enumerate(input.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		if len(sent) > maxlen: continue
		for word in sent:
			assert unknownwords or (
				word in coarse.lexicon and word in fine.lexicon), (
				"unknown word and no open class tags supplied")
		print "parsing:", n, " ".join(sent)
		start = ChartItem(coarse.toid["TOP"], 0, len(sent))
		chart, _ = parse(sent, coarse, None); print
		if chart[0][len(sent)].get(coarse.toid["TOP"]):
			trees = []
			candidates = lazykbest(chart, start, k, coarse.tolabel, sent)
			lexchart = doplexprobs(Tree(candidates[0][0]), fine)
			for m, (tree, prob) in enumerate(candidates):
				trees.append((dopparseprob(Tree(tree),
						fine, mapping, lexchart), tree))
				print m, exp(-prob), exp(trees[-1][0])
				sys.stdout.flush()
			results = nlargest(printk, trees)
			# print k-best parsetrees
			if printprob:
				out.writelines("parseprob=%.16g\n%s\n" % (exp(prob), tree)
						for prob, tree in results)
			else:
				out.writelines("%s\n" % tree for _, tree in results)
		else:
			unparsed += 1
			print "No parse"
			out.write("No parse for \"%s\"\n" % " ".join(sent))
		out.write("\n")
		times.append(time.clock())
		print times[-1] - times[-2], "s"
		out.flush()
	print "raw cpu time", time.clock() - times[0]
	times = [a - b for a, b in zip(times[1::2], times[::2])]
	print "average time per sentence", sum(times) / len(times)
	print "unparsed sentences:", unparsed
	print "finished"
	out.close()

def nonterminalmapping(coarse, fine):
	mapping = {}
	for a in coarse.tolabel:
		mapping[a] = set()
	for a, b in fine.toid.items():
		mapping[coarse.toid[removeids.sub("", a)]].add(b)
	return mapping

if __name__ == '__main__': main()
