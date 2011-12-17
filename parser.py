import sys, time
from operator import itemgetter
from math import exp
from heapq import nlargest
import numpy as np
from cky import parse, readbitpargrammar, pprint_chart
from kbest import lazykbest
from containers import ChartItem
from coarsetofine import whitelistfromkbest, whitelistfromposteriors
from disambiguation import marginalize

def mainsimple(unknownwords):
	grammar = readbitpargrammar(sys.argv[1], sys.argv[2], unknownwords)
	if len(sys.argv) >= 4: input = open(sys.argv[3])
	else: input = sys.stdin
	if len(sys.argv) == 5: out = open(sys.argv[4], "w")
	else: out = sys.stdout
	times = [time.clock()]
	for n, a in enumerate(input.read().split("\n\n")):
		if n == 1: break
		if not a.strip(): continue
		sent = a.splitlines()
		for word in sent:
			assert word in grammar.lexicon or unknownwords, "unknown word and no open class tags supplied"
		print "parsing:", n, " ".join(sent),
		sys.stdout.flush()
		chart, viterbi = parse(sent, grammar, None)
		start = ChartItem(grammar.toid["TOP"], 0, len(sent))
		#pprint_chart(chart, sent, grammar.tolabel)
		parsetrees = lazykbest(chart, start, 50, grammar.tolabel, sent)
		assert len(parsetrees) == len(set(parsetrees))
		assert len(parsetrees) == len(set(tree for tree, prob in parsetrees))
		if start in chart:
			out.writelines("vitprob=%.16g\n%s\n" % (exp(-prob), tree)
				for tree, prob in parsetrees)
		else:
			out.write("No parse for \"%s\"\n" % " ".join(sent))
		out.write("\n")
		out.flush()
		times.append(time.clock())
		print times[-1] - times[-2], "s"
	print "raw cpu time", time.clock() - times[0]
	times = [a - b for a, b in zip(times[1::2], times[::2])]
	print "average time per sentence", sum(times) / len(times)
	print "finished"
	out.close()

def mainctf(unknownwords):
	k = 50
	m = 10000
	threshold = 1e-5
	maxlen = 65
	unparsed = 0
	coarse = readbitpargrammar(sys.argv[1], sys.argv[2], unknownwords)
	fine = readbitpargrammar(sys.argv[3], sys.argv[4], unknownwords, normalize=False)
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
	finechart = np.array([np.NAN], dtype='d').repeat(
					maxlen * (maxlen + 1) * len(fine.toid)
					).reshape((len(fine.toid), maxlen, maxlen + 1))
	l=[]
	for n, a in enumerate(input.read().split("\n\n")):
		if not a.strip(): continue
		sent = a.splitlines()
		if len(sent) > maxlen: continue
		for word in sent:
			assert unknownwords or (word in coarse.lexicon and word in fine.lexicon), "unknown word and no open class tags supplied"
		coarsechart.fill(np.inf)
		print "parsing:", n, " ".join(sent)
		chart, coarsechart = parse(sent, coarse, coarsechart)
		start = ChartItem(coarse.toid["TOP"], 0, len(sent))
		if start in chart:
			print "getting whitelist",
			sys.stdout.flush()
			#finechart.fill(np.NAN)
			finechart[:,:len(sent),:len(sent)+1] = np.NAN
			#whitelistfromkbest(chart, start, coarse, fine, k, finechart, maxlen)
			whitelistfromposteriors(chart, coarsechart, start, coarse, fine, finechart, threshold)
			print "done"
			chart, finechart = parse(sent, fine, finechart)
			start = ChartItem(fine.toid["TOP"], 0, len(sent))
			assert start in chart, "sentence covered by coarse grammar could not be parsed by fine grammar"
			#pprint_chart(chart, sent, fine.tolabel)
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

def main():
	if "-u" in sys.argv:
		i = sys.argv.index("-u")
		unknownwords = sys.argv[i+1]
		sys.argv[i:i+2] = []
	else: unknownwords = None
	if 2 <= len(sys.argv) <= 5: mainsimple(unknownwords)
	elif 4 <= len(sys.argv) <= 7: mainctf(unknownwords)
	else: #if len(sys.argv) < 3:
		print "usage: %s rules lexicon [input [output]]" % sys.argv[0]
		print "or: %s coarserules coarselexicon finerules finelexicon [input [output]]" % sys.argv[0]
		return

if __name__ == '__main__': main()
