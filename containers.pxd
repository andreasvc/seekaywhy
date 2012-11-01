
cdef class ChartItem:
	cdef public unsigned int label
	cdef public short left, right

cdef class Edge:
	cdef public double inside
	cdef public Rule rule
	cdef public short split

cdef class RankedEdge:
	cdef public unsigned int label
	cdef public short left, right
	cdef public short leftrank, rightrank
	cdef public Edge edge

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

cdef inline ChartItem new_ChartItem(unsigned int label,
		short left, short right):
	cdef ChartItem item = ChartItem.__new__(ChartItem)
	item.label = label; item.left = left; item.right = right
	return item

cdef inline Edge new_Edge(double inside, Rule rule, short split):
	cdef Edge edge = Edge.__new__(Edge)
	edge.inside = inside; edge.rule = rule; edge.split = split
	return edge

cdef inline RankedEdge new_RankedEdge(unsigned int label,
		short left, short right, Edge edge, short j1, short j2):
	cdef RankedEdge rankededge = RankedEdge.__new__(RankedEdge)
	rankededge.label = label; rankededge.left = left; rankededge.right = right
	rankededge.leftrank = j1; rankededge.rightrank = j2
	rankededge.edge = edge;
	return rankededge

