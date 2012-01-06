import sys, time, re
from operator import itemgetter
from math import exp, log, isinf
from heapq import nlargest
import numpy as np
from nltk import Tree
from cky import parse, parse_nomatrix, readbitpargrammar, doinsideoutside, pprint_chart, pprint_matrix, dopparseprob, getgrammarmapping, cachingdopparseprob, doplexprobs
from kbest import lazykbest
from containers import ChartItem
from coarsetofine import whitelistfromkbest, whitelistfromposteriors, whitelistfromposteriors2
from disambiguation import marginalize
removeids = re.compile("@[0-9_]+")

def mainsimple(unknownwords):
	k = 1
	grammar = readbitpargrammar(sys.argv[1], sys.argv[2], unknownwords)
	if len(sys.argv) >= 4: input = open(sys.argv[3])
	else: input = sys.stdin
	if len(sys.argv) == 5: out = open(sys.argv[4], "w")
	else: out = sys.stdout
	times = [time.clock()]
	for n, a in enumerate(input.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		for word in sent:
			assert word in grammar.lexicon or unknownwords, "unknown word \"%s\" and no open class tags supplied" % word
		print "parsing:", n, " ".join(sent),
		sys.stdout.flush()
		chart, viterbi = parse(sent, grammar, None)
		start = ChartItem(grammar.toid["TOP"], 0, len(sent))
		#pprint_chart(chart, sent, grammar.tolabel)
		if chart[0][len(sent)].get(grammar.toid["TOP"], False):
			parsetrees = lazykbest(chart, start, k, grammar.tolabel, sent)
			assert len(parsetrees) == len(set(parsetrees))
			assert len(parsetrees) == len(set(tree for tree, prob in parsetrees))
			#out.writelines("vitprob=%.16g\n%s\n" % (exp(-prob), tree)
			out.writelines("%s\n" % tree
				for tree, prob in parsetrees)
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

def mainctf(unknownwords):
	k = 1000
	m = 10000
	threshold = exp(-6.2) #0.01
	maxlen = 65
	unparsed = 0
	coarse = readbitpargrammar(sys.argv[1], sys.argv[2], unknownwords, logprob=False)
	fine = readbitpargrammar(sys.argv[3], sys.argv[4], unknownwords, freqs=False)
	for a in fine.toid:
		assert a.split("@")[0] in coarse.toid, "%s not in coarse grammar" % a
	if len(sys.argv) >= 6: input = open(sys.argv[5])
	else: input = sys.stdin
	if len(sys.argv) == 7: out = open(sys.argv[6], "w")
	else: out = sys.stdout
	times = [time.clock()]
	coarsechart = np.array([np.inf], dtype='d').repeat(
					maxlen * (maxlen + 1) * len(coarse.toid)
					).reshape((len(coarse.toid), maxlen, maxlen + 1))
	inside = np.array([0.0], dtype='d').repeat(
					maxlen * (maxlen + 1) * len(coarse.toid)
					).reshape((len(coarse.toid), maxlen, maxlen + 1))
	outside = np.array([0.0], dtype='d').repeat(
					maxlen * (maxlen + 1) * len(coarse.toid)
					).reshape((len(coarse.toid), maxlen, maxlen + 1))
	#finechart = np.array([np.NAN], dtype='d').repeat(
	#				maxlen * (maxlen + 1) * len(fine.toid)
	#				).reshape((len(fine.toid), maxlen, maxlen + 1))
	mapping = nonterminalmapping(coarse, fine)

	for n, a in enumerate(input.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		if len(sent) > maxlen: continue
		for word in sent:
			assert unknownwords or (word in coarse.lexicon and word in fine.lexicon), "unknown word and no open class tags supplied"
		print "parsing:", n, " ".join(sent)
		start = ChartItem(coarse.toid["TOP"], 0, len(sent))
		#coarsechart[:,:len(sent),:len(sent)+1] = np.inf
		#chart, coarsechart = parse(sent, coarse, coarsechart)
		inside, outside = doinsideoutside(sent, coarse, inside, outside)
		#print "inside"; pprint_matrix(inside, sent, coarse.tolabel)
		#print "outside"; pprint_matrix(outside, sent, coarse.tolabel)
		#if chart[0][len(sent)].get(coarse.toid["TOP"], False):
		if inside[coarse.toid["TOP"], 0, len(sent)] != 0.0:
			print "pruning ...",
			sys.stdout.flush()
			#whitelistfromkbest(chart, start, coarse, fine, k, finechart, maxlen)
			finechart = whitelistfromposteriors2(inside, outside, start, coarse, fine, mapping, maxlen, threshold)
			#chart, finechart = parse(sent, fine, finechart)
			chart = parse_nomatrix(sent, fine, finechart)
			start = ChartItem(fine.toid["TOP"], 0, len(sent))
			#pprint_chart(chart, sent, fine.tolabel)
			#for l, _ in enumerate(chart):
			#	for r, _ in enumerate(chart[l]):
			#		for label in chart[l][r]:
			#			print fine.tolabel[label], label, l, r,
			#			if chart[l][r][label]: print chart[l][r][label][0]
			#			else: print []

			#assert start in chart, "sentence covered by coarse grammar could not be parsed by fine grammar"
			assert chart[0][len(sent)][fine.toid["TOP"]], "sentence covered by coarse grammar could not be parsed by fine grammar"
			# MPP
			parsetrees = marginalize(chart, start, fine.tolabel, sent, n=m).items()
			# print all parsetrees
			#out.writelines("parseprob=%.16g\n%s\n" % (prob, tree) for tree, prob in parsetrees)
			# most probable parse
			out.write("prob=%.16g\n%s\n" % max(parsetrees, key=itemgetter(1))[::-1])
			# MPD
			#parsetrees = marginalize(chart, start, fine.tolabel, sent, n=m, mpd=True).items()
			# print 10-best derivations
			#out.writelines("derivprob=%.16g\n%s\n" % (prob, tree) for tree, prob in nlargest(10, parsetrees, key=itemgetter(1)))

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

def mainrerank(unknownwords):
	k = 50
	maxlen = 65
	unparsed = 0
	coarse = readbitpargrammar(sys.argv[1], sys.argv[2], unknownwords, logprob=True)
	fine = readbitpargrammar(sys.argv[3], sys.argv[4], unknownwords, freqs=False, logprob=True)
	for a in fine.toid:
		assert a.split("@")[0] in coarse.toid, "%s not in coarse grammar" % a
	if len(sys.argv) >= 6: input = open(sys.argv[5])
	else: input = sys.stdin
	if len(sys.argv) == 7: out = open(sys.argv[6], "w")
	else: out = sys.stdout
	times = [time.clock()]
	mapping = getgrammarmapping(coarse, fine)
	for n, a in enumerate(input.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		if len(sent) > maxlen: continue
		for word in sent:
			assert unknownwords or (word in coarse.lexicon and word in fine.lexicon), "unknown word and no open class tags supplied"
		print "parsing:", n, " ".join(sent)
		start = ChartItem(coarse.toid["TOP"], 0, len(sent))
		chart, _ = parse(sent, coarse, None); print
		if chart[0][len(sent)].get(coarse.toid["TOP"], False):
			trees = []
			for m, (tree, prob) in enumerate(lazykbest(chart, start, k, coarse.tolabel, sent)):
				if m == 0: lexchart = doplexprobs(Tree(tree), fine)
				trees.append((dopparseprob(Tree(tree), fine, mapping, lexchart), tree))
				print m, exp(-prob), trees[-1][0] #, exp(trees[-1][0])
				sys.stdout.flush()
			prob, tree = max(trees)
			prob = exp(prob)
			print "\nprob=%.16g\n%s\n" % (prob, tree)
			out.write("prob=%.16g\n%s\n" % (prob, tree))
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

def main():
	if "-u" in sys.argv:
		i = sys.argv.index("-u")
		unknownwords = sys.argv[i+1]
		sys.argv[i:i+2] = []
	else: unknownwords = None
	if 2 <= len(sys.argv) <= 5: mainsimple(unknownwords)
	elif 4 <= len(sys.argv) <= 8:
		if sys.argv[1] == "rerank":
			del sys.argv[1]
			mainrerank(unknownwords)
		elif len(sys.argv) <= 7:
			mainctf(unknownwords)
	else: #if len(sys.argv) < 3:
		print "usage: %s rules lexicon [input [output]]" % sys.argv[0]
		print "or: %s [rerank] coarserules coarselexicon finerules finelexicon [input [output]]" % sys.argv[0]
		return

if __name__ == '__main__': main()
