
cdef class ChartItem:
	def __init__(self, label, left, right):
		self.label = label
		self.left = left
		self.right = right
	def __hash__(ChartItem self):
		cdef long h
		# juxtapose bits of label and left, right
		h = self.label ^ (<long>self.left << 32UL) ^ (<long>self.right << 40UL)
		return -2 if h == -1 else h
	def __richcmp__(ChartItem self, ChartItem other, int op):
		if op == 2: return self.label == other.label and self.left == other.left and self.right == other.right
		elif op == 3: return self.label != other.label or self.left != other.left or self.right != other.right
		elif op == 5: return self.label >= other.label or self.left >= other.left or self.right >= other.right
		elif op == 1: return self.label <= other.label or self.left <= other.left or self.right <= other.right
		elif op == 0: return self.label < other.label or self.left < other.left or self.right < other.right
		elif op == 4: return self.label > other.label or self.left > other.left or self.right > other.right
	def __nonzero__(ChartItem self):
		return self.label != 0 and self.right != 0
	def __repr__(ChartItem self):
		return "ChartItem(%d, %d, %d)" % (self.label, self.left, self.right)

cdef class Edge:
	def __init__(self, inside, rule, split):
		self.inside = inside; self.rule = rule; self.split = split
	def __hash__(self):
		cdef long h
		#self._hash = hash((inside, prob, left, right))
		# this is the hash function used for tuples, apparently
		h = (1000003UL * 0x345678UL) ^ <long>self.inside
		h = (1000003UL * h) ^ <long>(<Rule>self.rule).prob
		h = (1000003UL * h) ^ (<Rule>self.rule).rhs1
		h = (1000003UL * h) ^ (<Rule>self.rule).rhs2
		h = (1000003UL * h) ^ <long>self.split
		return -2 if h == -1 else h
	def __richcmp__(Edge self, other, int op):
		# the ordering only depends on the estimate / inside score
		if op == 0: return self.inside < (<Edge>other).inside
		elif op == 1: return self.inside <= (<Edge>other).inside
		# (in)equality compares all elements
		# boolean trick: equality and inequality in one expression i.e., the
		# equality between the two boolean expressions acts as biconditional
		elif op == 2 or op == 3:
			return (op == 2) == (
				self.inside == (<Edge>other).inside
				and self.rule == (<Edge>other).rule
				and self.split == (<Edge>other).split)
		elif op == 4: return self.inside > other.inside
		elif op == 5: return self.inside >= other.inside
	def __repr__(self):
		return "Edge(%g, %r, %r)" % (
					self.inside, self.rule, self.split)

cdef class RankedEdge:
	def __init__(self, head, edge, left, right):
		self.head = head; self.edge = edge
		self.left = left; self.right = right
	#def __cinit__(self, ChartItem head, Edge edge, int j1, int j2):
	#	self.head = head; self.edge = edge
	#	self.left = j1; self.right = j2
	def __hash__(self):
		cdef long h
		#h = hash((head, edge, j1, j2))
		h = (1000003UL * 0x345678UL) ^ hash(self.head)
		h = (1000003UL * h) ^ hash(self.edge)
		h = (1000003UL * h) ^ self.left
		h = (1000003UL * h) ^ self.right
		if h == -1: h = -2
		return h
	def __richcmp__(self, RankedEdge other, int op):
		if op == 2 or op == 3:
			return (op == 2) == (
				self.left == other.left
				and self.right == other.right
				and self.head == other.head
				and self.edge == other.edge)
		else:
			raise NotImplemented
	def __repr__(self):
		return "RankedEdge(%r, %r, %d, %d)" % (
					self.head, self.edge, self.left, self.right)

cdef inline RankedEdge new_RankedEdge(ChartItem head, Edge edge, short j1, short j2):
	cdef RankedEdge rankededge = RankedEdge.__new__(RankedEdge)
	rankededge.head = head; rankededge.edge = edge;
	rankededge.left = j1; rankededge.right = j2
	return rankededge

cdef inline ChartItem new_ChartItem(unsigned int label, short left, short right):
	cdef ChartItem item = ChartItem.__new__(ChartItem)
	item.label = label; item.left = left; item.right = right
	return item

cdef class Grammar:
	def __init__(self, lexical, unary, unarybyrhs, binary, tolabel, toid,
			lexicon, logprob):
		self.lexical = lexical; self.unary = unary; self.binary = binary
		self.tolabel = tolabel; self.toid = toid; self.lexicon = lexicon
		self.unarybyrhs = unarybyrhs; self.logprob = logprob
	def __repr__(self):
		return ("Grammar with %d nonterminals, %d unary productions, "
			"%d binary productions, %d lexical productions, and "
			"% words in lexicon. Logprobs: %r" % (len(self.toid), len(self.unary),
			sum(map(len, self.binary)), len(self.lexical), len(self.lexicon),
			self.logprob))

cdef class Rule:
	def __init__(self, lhs, rhs1, rhs2, prob):
		self.lhs = lhs; self.rhs1 = rhs1; self.rhs2 = rhs2
		self.prob = prob
	def __repr__(self):
		return "Rule(%r, %r, %r, %r)" % (
			self.lhs, self.rhs1, self.rhs2, self.prob)

cdef class Terminal(Rule):
	def __init__(self, lhs, word, prob):
		self.lhs = lhs; self.word = word; self.prob = prob
		self.rhs1 = 0; self.rhs2 = 0
	def __repr__(self):
		return "Terminal(%r, %r, %r, %r)" % (
			self.lhs, self.rhs1, self.rhs2, self.prob)
