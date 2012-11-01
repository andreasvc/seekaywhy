
all:
	python setup.py build_ext --inplace

clean:
	rm -f *.c *.so

test: all
	python parser.py -b 5 --prob t.rules t.lexicon t.test
	python parser.py -b 5 --prob --threshold -6.2 t.rules t.lexicon tdop.rules tdop.lexicon t.test
	python parser.py -b 5 --prob --kbestctf 50 t.rules t.lexicon tdop.rules tdop.lexicon t.test
	python parser.py -b 5 --prob --rerank 50 t.rules t.lexicon tdop.rules tdop.lexicon t.test
	python parser.py -b 5 --prob --reestimate t.rules t.lexicon tdop.rules tdop.lexicon t.test
	python parser.py -b 5 --prob --rerank 50 --reestimate t.rules t.lexicon tdop.rules tdop.lexicon t.test
	python parser.py -b 5 --prob --mpd t.rules t.lexicon tdop.rules tdop.lexicon t.test

