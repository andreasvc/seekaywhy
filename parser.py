import re, time
from sys import argv, stdin, stdout, stderr
from math import exp
from getopt import gnu_getopt
from heapq import nlargest
import numpy as np
from nltk import Tree
from cky import parse, parse_sparse, readbitpargrammar, doinsideoutside, \
	getgrammarmapping, doplexprobs, dopparseprob, reestimate
#	pprint_matrix, cachingdopparseprob
from kbest import lazykbest
from coarsetofine import whitelistfromkbest, whitelistfromposteriors2
from disambiguation import marginalize
removeids = re.compile("@[0-9_]+")

usage = """usage: %s [options] rules lexicon [input [output]]
or: %s [options] coarserules coarselexicon finerules finelexicon \
[input [output]]

Grammars need to be binarized, and are in bitpar format.
When no file is given, output is written to standard output; when additionally
no input is given, it is raed from standard input.

Options:
-u file       Handle unknown words using given open class tags list.
-b k          Return the k-best parses instead of just 1.
-s x          Use "x" as start symbol instead of default "TOP".
--threshold p In coarse-to-fine mode, set the posterior threshold (as logprob).
--kbestctf k  Use k-best coarse-to-fine; instead of posterior threshold,
              prune items not in k-best derivations.
--rerank k    Enable DOP reranking mode: find DOP parse probabilities
              for k parse trees from the coarse grammar.
--prob        Print probabilities as well as parse trees.
--reestimate  Make probabilities of coarse grammar reflect sum of the relevant
              fine grammar probabilities.
--mpd         In coarse-to-fine mode, produce the most probable derivation (MPD)
              instead of the most probable parse (MPP).
""" % (argv[0], argv[0])

def main():
	print >>stderr, "SeeKayWhy PCFG parser - Andreas van Cranenburgh"
	options = "rerank= threshold= kbestctf= prob reestimate mpd".split()
	opts, args = gnu_getopt(argv[1:], "u:b:s:", options)
	opts = dict(opts)
	k = int(opts.get("-b", 1))
	start = opts.get("-s", "TOP")
	unknownwords = opts.get("-u")
	prob = "--prob" in opts
	if 2 <= len(args) <= 4:
		simple(args, unknownwords, k, prob, start)
	elif 4 <= len(args) <= 6 and "--rerank" in opts:
		rerank(args, unknownwords, k, prob, start, int(opts['--rerank']),
				"--reestimate" in opts)
	elif 4 <= len(args) <= 6:
		posterior = "--kbestctf" not in opts
		if posterior: threshold = float(opts.get("--threshold", -6.2)) #0.01
		else: threshold = int(opts["--kbestctf"])
		ctf(args, unknownwords, k, prob, start, threshold, posterior,
				"--reestimate" in opts, "--mpd" in opts)
	else: print usage

def simple(args, unknownwords, k, printprob, start):
	grammar = readbitpargrammar(args[0], args[1], unknownwords)
	assert start in grammar.toid, "Start symbol %r not in grammar." % start
	infile = open(args[2]) if len(args) >= 3 else stdin
	out = open(args[3], "w") if len(args) == 4 else stdout
	times = [time.clock()]
	for n, a in enumerate(infile.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		for word in sent:
			assert word in grammar.lexicon or unknownwords, (
				"unknown word %r and no open class tags supplied" % word)
		print >>stderr, "parsing:", n, " ".join(sent),
		stdout.flush()
		chart, _ = parse(sent, grammar, None)
		if chart[0][len(sent)].get(grammar.toid[start]):
			parsetrees = lazykbest(chart, grammar.toid[start], 0, len(sent),
					k, grammar.tolabel)
			assert len(parsetrees) == len(set(parsetrees))
			assert len(parsetrees) == len(set(tree for tree, _ in parsetrees))
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
		print >>stderr, times[-1] - times[-2], "s"
	print >>stderr, "raw cpu time", time.clock() - times[0]
	times = [a - b for a, b in zip(times[1::2], times[::2])]
	print >>stderr, "average time per sentence", sum(times) / len(times)
	print >>stderr, "finished"
	out.close()

def ctf(args, unknownwords, k, printprob, start, threshold, posterior,
		doreestimate, mpd):
	m = 10000 # number of derivations from fine grammar to marginalize
	if posterior: threshold = exp(threshold)
	maxlen = 999 #65
	unparsed = 0
	coarse = readbitpargrammar(args[0], args[1], unknownwords, logprob=not posterior)
	fine = readbitpargrammar(args[2], args[3], unknownwords, freqs=False)
	for a in fine.toid:
		assert a.rsplit("@", 1)[0] in coarse.toid, "%s not in coarse grammar" % a
	assert start in fine.toid, "Start symbol %r not in grammar." % start
	if doreestimate: reestimate(coarse, fine)
	infile = open(args[4]) if len(args) >= 5 else stdin
	out = open(args[5], "w") if len(args) == 6 else stdout
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

	for n, a in enumerate(infile.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		if len(sent) > maxlen: continue
		for word in sent:
			assert unknownwords or (
				word in coarse.lexicon and word in fine.lexicon), (
				"unknown word and no open class tags supplied")
		print >>stderr, "parsing:", n, " ".join(sent)
		if posterior:
			inside, outside = doinsideoutside(sent, coarse, inside, outside)
			#print "inside"; pprint_matrix(inside, sent, coarse.tolabel)
			#print "outside"; pprint_matrix(outside, sent, coarse.tolabel)
		else:
			coarsechart[:len(sent), :len(sent)+1, :] = np.inf
			chart, coarsechart = parse(sent, coarse, coarsechart)
		if posterior: goalitem = inside[0, len(sent), coarse.toid[start]]
		else: goalitem = chart[0][len(sent)].get(coarse.toid[start])
		if goalitem:
			print >>stderr, "pruning ...",
			stdout.flush()
			if posterior:
				finechart = whitelistfromposteriors2(inside, outside,
					coarse.toid[start], len(sent), coarse, fine, mapping, threshold)
			else:
				finechart = whitelistfromkbest(chart, coarse.toid[start], len(sent),
					coarse, fine, threshold, mapping)
			#chart, finechart = parse(sent, fine, finechart)
			chart = parse_sparse(sent, fine, finechart)

			assert chart[0][len(sent)][fine.toid[start]], (
				"sentence covered by coarse grammar could not be parsed "\
				"by fine grammar")
			parsetrees = marginalize(chart, fine.toid[start], fine.tolabel,
					sent, n=m, mpd=mpd)
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
			print >>stderr, "No parse"
			out.write("No parse for \"%s\"\n" % " ".join(sent))
		out.write("\n")
		times.append(time.clock())
		print >>stderr, times[-1] - times[-2], "s"
		out.flush()
	print >>stderr, "raw cpu time", time.clock() - times[0]
	times = [a - b for a, b in zip(times[1::2], times[::2])]
	print >>stderr, "average time per sentence", sum(times) / len(times)
	print >>stderr, "unparsed sentences:", unparsed
	print >>stderr, "finished"
	out.close()

def rerank(args, unknownwords, printk, printprob, start, k, doreestimate):
	maxlen = 999 #??
	unparsed = 0
	coarse = readbitpargrammar(args[0], args[1], unknownwords, logprob=True)
	fine = readbitpargrammar(args[2], args[3], unknownwords, freqs=False, logprob=True)
	for a in fine.toid:
		assert a.rsplit("@", 1)[0] in coarse.toid, "%s not in coarse grammar" % a
	assert start in fine.toid, "Start symbol %r not in grammar." % start
	if doreestimate: reestimate(coarse, fine)
	infile = open(args[4]) if len(args) >= 5 else stdin
	out = open(args[5], "w") if len(args) == 6 else stdout
	times = [time.clock()]
	mapping = getgrammarmapping(coarse, fine)
	for n, a in enumerate(infile.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		if len(sent) > maxlen: continue
		for word in sent:
			assert unknownwords or (
				word in coarse.lexicon and word in fine.lexicon), (
				"unknown word and no open class tags supplied")
		print >>stderr, "parsing:", n, " ".join(sent)
		chart, _ = parse(sent, coarse, None); print >>stderr, ''
		if chart[0][len(sent)].get(coarse.toid[start]):
			trees = []
			candidates = lazykbest(chart, coarse.toid[start], 0, len(sent),
					k, coarse.tolabel)
			lexchart = doplexprobs(Tree(candidates[0][0]), fine)
			for m, (tree, prob) in enumerate(candidates):
				trees.append((dopparseprob(Tree(tree),
						fine, mapping, lexchart), tree))
				print >>stderr, m, exp(-prob), exp(trees[-1][0])
				stdout.flush()
			results = nlargest(printk, trees)
			# print k-best parsetrees
			if printprob:
				out.writelines("parseprob=%.16g\n%s\n" % (exp(prob), tree)
						for prob, tree in results)
			else:
				out.writelines("%s\n" % tree for _, tree in results)
		else:
			unparsed += 1
			print >>stderr, "No parse"
			out.write("No parse for \"%s\"\n" % " ".join(sent))
		out.write("\n")
		times.append(time.clock())
		print >>stderr, times[-1] - times[-2], "s"
		out.flush()
	print >>stderr, "raw cpu time", time.clock() - times[0]
	times = [a - b for a, b in zip(times[1::2], times[::2])]
	print >>stderr, "average time per sentence", sum(times) / len(times)
	print >>stderr, "unparsed sentences:", unparsed
	print >>stderr, "finished"
	out.close()

def nonterminalmapping(coarse, fine):
	mapping = {}
	for a in coarse.tolabel:
		mapping[a] = set()
	for a, b in fine.toid.items():
		mapping[coarse.toid[removeids.sub("", a)]].add(b)
	return mapping

if __name__ == '__main__': main()
