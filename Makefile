
all:
	python setup.py build_ext --inplace

clean:
	rm -f *.c *.so

test: all
	python parser.py t.rules t.lexicon t.test
	python parser.py t.rules t.lexicon tdop.rules tdop.lexicon t.test
	python parser.py t.rules t.lexicon tdop.rules tdop.lexicon t.test --rerank 50
	#python parser.py t.rules t.lexicon tdop.rules tdop.lexicon t.test --kbestctf 50

