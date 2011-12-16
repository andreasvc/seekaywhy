
cdef class ChartItem:
	cdef public unsigned int label
	cdef public short left, right

cdef class Edge:
	cdef public double score
	cdef public double inside
	cdef public Rule rule
	cdef public short split

cdef class RankedEdge:
	cdef public ChartItem head
	cdef public Edge edge
	cdef public short left, right

cdef class Grammar:
	cdef public list lexical, unary, binary
	cdef public dict toid, tolabel
	cdef public set lexicon

cdef class Rule:
	cdef public unsigned int lhs, rhs1, rhs2
	cdef public double prob

cdef class Terminal(Rule):
	cdef public unicode word

cpdef inline unsigned int getlabel(ChartItem a)
cpdef inline unsigned long long getvec(ChartItem a)
cpdef inline double getscore(Edge a)
cpdef inline dict dictcast(d)
cpdef inline ChartItem itemcast(i)
cpdef inline Edge edgecast(e)
cdef inline RankedEdge new_RankedEdge(ChartItem head, Edge edge, short j1, short j2)
