SeeKayWhy
=========

A probabilistic CKY parser for PCFGs.

- accepts arbitrary binary grammars, in the same format as used by Bitpar.
- Input and output is also similar to Bitpar: input is one word per line,
  sentences separated by a blank line.
- produces exact k-best lists for arbitrary k, on the basis of
  an exhaustive chart
- coarse-to-fine parsing: parse first with a coarse grammar with symbols of
  the form A, then parse with fine grammar with symbols of the form A@x (WIP).
- reranking with DOP reductions: compute the exact DOP parse probability of
  the derivations from a coarse grammar.

NB: all of these features are now incorporated in disco-dop;
cf. https://github.com/andreasvc/disco-dop

Requirements:
-------------
- Python 2.6+   http://www.python.org (need headers, e.g. python-dev package)
- Cython 0.15+  http://www.cython.org
- GCC           http://gcc.gnu.org/
- NLTK          http://www.nltk.org
- Numpy         http://numpy.scipy.org/

For example, to install these dependencies and compile the code on Ubuntu
(tested on 12.04), run the following sequence of commands:

    sudo apt-get install cython python-dev python-nltk python-numpy build-essential
    git clone --depth 1 git://github.com/andreasvc/seekaywhy.git
    cd seekaywhy
    make

Alternatively Cython, NLTK, and Numpy can all be installed with
`pip install cython nltk numpy`,
which does not require root rights and may be more up-to-date.

Example invocation:
-------------------

	$ python parser.py t.rules t.lexicon t.test
	[...]
	vitprob=0.002057407417695131
	(TOP (S (NP (PN John))(VP (V sees)(NP (NP (DT the)(NN boy))(PP (IN with)(NP (DT the)(NN telescope)))))))

An example of coarse-to-fine with a PCFG reduction of DOP as the fine grammar:

	$ python parser.py t.rules t.lexicon tdop.rules tdop.lexicon t.test
	[...]
	prob=5.715931040919541e-45
	(TOP (S (NP (PN John))(VP (VP (V sees)(NP (DT the)(NN boy)))(PP (IN with)(NP (DT the)(NN telescope))))))

Reranking PCFG parses using a DOP reduction as the fine grammar:

	python parser.py --rerank t.rules t.lexicon tdop.rules tdop.lexicon t.test
	[...]
	prob=5.055473074786365e-16
	(TOP (S (NP (PN John))(VP (VP (V sees)(NP (DT the)(NN boy)))(PP (IN with)(NP (DT the)(NN telescope))))))

