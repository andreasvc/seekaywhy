
cdef class ChartItem:
	cdef public unsigned int label
	cdef public short left, right

cdef class Edge:
	cdef public double inside
	cdef public Rule rule
	cdef public short split

cdef class RankedEdge:
	cdef public ChartItem head
	cdef public Edge edge
	cdef public short left, right

cdef class Grammar:
	cdef public list lexical, unary, unarybyrhs, binary
	cdef public dict toid, tolabel
	cdef public set lexicon
	cdef public bint logprob

cdef class Rule:
	cdef public unsigned int lhs, rhs1, rhs2
	cdef public double prob

cdef class Terminal(Rule):
	cdef public unicode word

cdef inline RankedEdge new_RankedEdge(ChartItem head, Edge edge, short j1, short j2)
cdef inline ChartItem new_ChartItem(unsigned int label, short left, short right)
