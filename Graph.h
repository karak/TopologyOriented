#pragma once

#include "Vector2.h"
#include "Geometry.h"

typedef std::size_t index_t;
static const index_t NOTHING = -1;

class RegionItemBase
{
public:
	RegionItemBase() : he_(NOTHING)
	{
	}

	RegionItemBase(index_t he) : he_(he)
	{
	}


	index_t boundingHalfedge() const
	{
		return he_;
	}

private:
	index_t he_;	//one of bounding halfedge CCW
};


class RegionItem : public RegionItemBase
{
public:
	RegionItem() : RegionItemBase(), site_()
	{
	}

	//of Euclid plane
	RegionItem(index_t he, const Vector2& site) : RegionItemBase(he), site_(site)
	{
	}

	const Vector2C site() const
	{
		return Vector2C(site_);
	}

private:
	Vector2 site_;
};

class RegionItemC : public RegionItemBase
{
public:
	RegionItemC() : RegionItemBase(), site_()
	{
	}

	//of infinite plane
	RegionItemC(index_t he) : RegionItemBase(he), site_(Vector2C::infinity())
	{
	}

	const Vector2C& site() const
	{
		return site_;
	}

private:
	Vector2C site_;
};

class EdgeItem
{
public:
	/***
	 *  @param vs start vertex id
	 *  @param ve end vertex id
	 *  @param rl left region id
	 *  @param rr right region id
	 *  @param b  bisector cache
	 */
	EdgeItem(index_t vs, index_t ve, index_t rl, index_t rr, const Line2& b) : bisector_(b)
	{
		v[0] = vs;
		v[1] = ve;
		r[0] = rl;
		r[1] = rr;
	}

	EdgeItem() : bisector_()
	{
		v[0] = v[1] = NOTHING;
		r[0] = r[0] = NOTHING;
	}

	EdgeItem(const EdgeItem& rhs) : bisector_(rhs.bisector_)
	{
		v[0] = rhs.v[0];
		v[1] = rhs.v[1];
		r[0] = rhs.r[0];
		r[1] = rhs.r[1];
	}

	index_t startVertex() const { return v[0]; }
	index_t endVertex() const { return v[1]; }
	const Line2& bisector() const { return bisector_; }
	
	static index_t forwardHalfedge(index_t e)
	{
		return e << 1;
	}

	static index_t backwardHalfedge(index_t e)
	{
		return e << 1 | 1;
	}

	static index_t forwardBit(index_t he)
	{
		return he & 1;
	}

	static index_t backwardBit(index_t he)
	{
		return forwardBit(he) ^ 1;
	}

	static index_t parentEdge(index_t he)
	{
		return he >> 1;
	}

private:
	friend class RootGraph;
	index_t v[2];
	index_t r[2];
	Line2 bisector_;
};

class VertexItem
{
public:
	static const std::size_t DEGREE = 3;

	VertexItem()
	{
		std::fill(halfedges, halfedges + DEGREE, NOTHING);
	}

	VertexItem(const VertexItem& rhs) 
	{
		std::copy(rhs.halfedges, rhs.halfedges + DEGREE, halfedges);
	}
	//out-halfedges in CCW-order	
	VertexItem(index_t he0, index_t he1, index_t he2, const Circumcircle& circumcircle) : circumcircle(circumcircle)
	{
		halfedges[0] = he0;
		halfedges[1] = he1;
		halfedges[2] = he2;
	}

	index_t halfedgeCW(index_t he) const
	{
		return halfedges[(indexOf(he) + 2) % DEGREE];
	}

	index_t halfedgeCCW(index_t he) const
	{
		return halfedges[(indexOf(he) + 1) % DEGREE];
	}

	index_t halfedgeAround(std::size_t i) const
	{
		return halfedges[i];
	}

	template<typename FuncT>
	FuncT forEachHalfedge(FuncT f) const
	{
		return std::for_each(halfedges, halfedges + DEGREE, f);
	}

	//out-halfedges in CCW order
	index_t halfedges[DEGREE];
	Circumcircle circumcircle;

private:
	int indexOf(index_t he) const;
};

#include <array>

class RootGraph
{
public:
	RootGraph(const Vector2& p0, const Vector2& p1, const Vector2& p2);
	RootGraph(const RootGraph& rhs);

public:
	/* edge/halfedge */
	//structure
	index_t mateHalfedge(index_t he) const
	{
		return (he ^ 0x1);
	}

	index_t startVertexOfEdge(index_t e) const
	{
		return E[e].startVertex();
	}

	index_t endVertexOfEdge(index_t e) const
	{
		return E[e].endVertex();
	}

	index_t startVertexOfHalfedge(index_t he) const
	{
		index_t edge = he >> 1;
		index_t outBit = he & 0x1;
		return E[edge].v[outBit];
	}

	index_t endVertexOfHalfedge(index_t he) const
	{
		return startVertexOfHalfedge(mateHalfedge(he));
	}

	index_t nextHalfedge(index_t he) const
	{
		return mateHalfedge(
			halfedgeCWAroundVertex(startVertexOfHalfedge(he), mateHalfedge(he))
		);
	}

	index_t parentEdge(index_t he) const
	{
		return EdgeItem::parentEdge(he);
	}

	index_t forwardHalfedge(index_t e) const
	{
		return EdgeItem::forwardHalfedge(e);
	}

	index_t backwardHalfedge(index_t e) const
	{
		return EdgeItem::backwardHalfedge(e);
	}

	index_t leftRegionOfHalfedge(index_t he) const
	{
		index_t e = parentEdge(he);
		return E[e].r[EdgeItem::forwardBit(he)];
	}

	index_t rightRegionOfHalfedge(index_t he) const
	{
		index_t e = parentEdge(he);
		return E[e].r[EdgeItem::backwardBit(he)];
	}

	template<typename FuncT>
	FuncT forEachEdge(FuncT f) const
	{
		for (index_t e = 0; e < Ne; ++e)	f(e);
		return f;
	}

	/* property */
	const Line2& bisector(index_t e) const
	{
		return E[e].bisector();
	}

	index_t maxIdOfEdge() const
	{
		return Ne;
	}

public:
	/* vertex */
	index_t halfedgeCWAroundVertex(index_t v, index_t he) const
	{
		return V[v].halfedgeCW(he);
	}

	index_t halfedgeCCWAroundVertex(index_t v, index_t he) const
	{
		return V[v].halfedgeCCW(he);
	}

	index_t halfEdgeAroundVertex(index_t v, std::size_t i) const
	{
		return V[v].halfedgeAround(i);
	}

	index_t regionAroundVertex(index_t v, std::size_t i) const
	{
		return leftRegionOfHalfedge(halfEdgeAroundVertex(v, i));
	}

	template<typename FuncT>
	FuncT forEachHalfedgeAroundVertex(index_t v, FuncT f) const
	{
		return V[v].forEachHalfedge(f);
	}

	/* property */
	const Circumcircle& circumcircle(index_t v) const
	{
		return V[v].circumcircle;
	}

	index_t maxIdOfVertex() const
	{
		return Nv;
	}
public:
	/* region */
	template<typename FuncT>
	FuncT forEachHalfedgeAroundRegion(index_t r, FuncT f) const
	{
		const auto h0 = boundingHalfedge(r);
		auto h = h0;
		do
		{
			f(h);
			h = halfedgeCCWAroundRegion(r, h);
		} while (h != h0);
		return f;
	}

	template<typename FuncT>
	FuncT forEachVertexAroundRegion(index_t r, FuncT f) const
	{
		forEachHalfedgeAroundRegion(r, [&](index_t he)
		{
			f(startVertexOfHalfedge(he));
		});
		return f;
	}
	
	index_t halfedgeCCWAroundRegion(index_t r, index_t he) const
	{
		return halfedgeCWAroundVertex(
			endVertexOfHalfedge(he),
			mateHalfedge(he)
		);
	}

	template<typename FuncT>
	FuncT forEachRegion(FuncT f) const
	{
		for (index_t r = 0; r < Nr; ++r)	f(r);
		return f;
	}

	const Vector2C site(index_t r) const
	{
		return r == 0? RInf_.site() : R_[r-1].site();
	}

	index_t maxIdOfRegion() const
	{
		return Nr;
	}

private:
	const index_t boundingHalfedge(index_t r) const
	{
		return regionItem(r).boundingHalfedge();
	}

	const RegionItemBase& regionItem(index_t r) const
	{
		return r == 0? static_cast<const RegionItemBase&>(RInf_) : R_[r-1];
	}

public:
	static const std::size_t Nr = 4;
	static const std::size_t NrEuclid = Nr - 1;
	static const std::size_t Ne = 6;
	static const std::size_t Nv = 4;
	RegionItemC RInf_;
	std::array<RegionItem, NrEuclid> R_;
	std::array<EdgeItem, Ne> E;
	std::array<VertexItem, Nv> V;
};


#include <hash_set>
#include <cassert>


template<class GraphT>
class InducedSubgraph
{
public:
	/**
	 * @param G target graph
	 * @param root id of vertex to start from
	 * @param p predicate which meet any vertex in subgraph must meet
	 */
	template<typename PredVertex>
	InducedSubgraph(const GraphT& G, const Vector2& site, index_t root, PredVertex pred)
		: site_(site), boundaryHalfedges_(), bisectors_(), circumcircles_()
	{
		std::hash_set<index_t> visitedVertices;
		std::vector<index_t> tmpBoundaryHalfEdges;
		dfs0(G, root, pred, visitedVertices, tmpBoundaryHalfEdges);
		construct(G, site, tmpBoundaryHalfEdges);
	}

	InducedSubgraph(InducedSubgraph&& rhs)
		: site_(rhs.site_)
		, n_(rhs.n_)
		, boundaryHalfedges_(rhs.boundaryHalfedges_)
		, bisectors_(rhs.bisectors_)
		, circumcircles_(rhs.circumcircles_)
	{
		rhs.circumcircles_ =  0;
		rhs.bisectors_ = 0;
		rhs.boundaryHalfedges_ = 0;
		rhs.n_ = 0;
	}

private:
	//avoid compiler bug(miss-matching)
	InducedSubgraph(const InducedSubgraph&);
	InducedSubgraph& operator=(const InducedSubgraph&);

public:
	~InducedSubgraph()
	{
		delete[] circumcircles_; circumcircles_ = 0;
		delete[] bisectors_; bisectors_ = 0;
		delete[] boundaryHalfedges_; boundaryHalfedges_ = 0;
	}
	
	const Vector2C site() const
	{
		return Vector2C(site_);
	}

	const std::size_t cycleLength() const
	{
		return n_;
	}
	
	index_t halfedgeAroundVertex(index_t v, index_t i) const
	{
		assert(i < 3);
		index_t hes[3];
		halfedgesAroundVertex(v, hes);
		return hes[i];
	}
	
	void halfedgesAroundVertex(index_t v, index_t hes[3]) const
	{
		hes[0] = boundaryHalfedges_[v];
		hes[1] = v;
		hes[2] = (v - 1) % n_;
	}

	/** @attention return as id of super-graph */
	index_t boundaryHalfedge(index_t vs) const
	{
		return boundaryHalfedges_[vs];
	}
	
	//for PrintableGraph
	template<typename FuncT>
	FuncT forEachRegion(FuncT f) const
	{
		f(static_cast<index_t>(0));
		return (f);
	}

	//for PrintableGraph
	const Vector2C site(index_t r) const
	{
		assert(r == 0);
		return site();
	}

	//for PrintableGraph	
	template<typename FuncT>
	FuncT forEachEdge(FuncT f) const
	{
		for (index_t e = 0; e < n_; ++e)	f(e);
		return (f);
	}

	//for PrintableGraph
	index_t startVertexOfEdge(index_t e) const
	{
		return e;
	}

	//for PrintableGraph
	index_t endVertexOfEdge(index_t e) const
	{
		return (e + 1) % n_;
	}

	//for PrintableGraph
	const Circumcircle& circumcircle(index_t v) const
	{
		return circumcircles_[v];
	}

	//for PrintableGraph
	const Line2& bisector(index_t e) const
	{
		return bisectors_[e];
	}

private:
	void construct(const GraphT& G, const Vector2& site, const std::vector<index_t>& tmpBoundaryHalfEdges)
	{
		n_ = tmpBoundaryHalfEdges.size();
		if (n_ > 0)
		{
			boundaryHalfedges_ = new index_t[n_];
			bisectors_ = new Line2[n_];
			circumcircles_ = new Circumcircle[n_];
			//TODO: error handling
			for (index_t i = 0; i < n_; ++i)
			{
				const auto he = tmpBoundaryHalfEdges[i];
				const auto rightSite = G.site(G.rightRegionOfHalfedge(he));
				const auto leftSite = G.site(G.leftRegionOfHalfedge(he));
				const auto b = ::bisector(site, leftSite.asEuclidPoint());	//TODO: Conformal‘Î‰ž
				const auto c = ::circumcircle(Vector2C(site), rightSite, leftSite);
				std::clog << "b="<<b<<'\n';
				std::clog << "q="<<c.center()<<'\n';
				boundaryHalfedges_[i] = he;
				bisectors_[i] = b;
				circumcircles_[i] = c;
			}
		}
	}

private:
	Vector2 site_;
	index_t n_;	//number of edges/vertices
	index_t* boundaryHalfedges_;
	Line2* bisectors_;
	Circumcircle* circumcircles_;

private:
	//depth first search
	template<typename PredVertex>
	static void dfs0(
		const GraphT& G,
		index_t root,
		PredVertex p,
		std::hash_set<index_t>& visitedEdges,
		std::vector<index_t>& boundaryHalfEdges
	)
	{
		const bool predResult = p(root);
		std::clog << IdString('V', root).c_str() << (predResult? 'o':'x') << '\n';

		G.forEachHalfedgeAroundVertex(root, [&](index_t he)
		{
			InducedSubgraph<GraphT>::dfs(
				G,
				he,
				predResult,
				p,
				visitedEdges,
				boundaryHalfEdges
			);
		});
	}

	//depth first search
	template<typename PredVertex>
	static void dfs(
		const GraphT& G,
		index_t currentHalfedge,
		bool prevPredResult,
		PredVertex p,
		std::hash_set<index_t>& visitedEdges,
		std::vector<index_t>& boundaryHalfEdges
	)
	{
		const auto currentEdge = G.parentEdge(currentHalfedge);
		std::clog << "|_" << IdString('E', currentEdge).c_str() << "=>\n";

		const auto range = visitedEdges.equal_range(currentEdge);
		if (range.first != range.second)
		{
			std::clog << "<=\n";
			return;
		}

		visitedEdges.insert(range.first, currentEdge);	//update

		const auto currentVertex = G.endVertexOfHalfedge(currentHalfedge);
		const bool predResult = p(currentVertex);

		std::clog << "   ." << IdString('V',currentVertex).c_str() << (predResult? 'o':'x') << '\n';

		if (prevPredResult != predResult)
		{
			if (predResult)
			{
				std::clog << "   BOUNDARY " << IdString('H', currentHalfedge) << '\n';
				std::clog << "<=\n";
				boundaryHalfEdges.push_back(currentHalfedge);
				return;	//end
			}
			else
			{
				const auto backHalfedge = G.mateHalfedge(currentHalfedge);
				std::clog << "   BOUNDARY " << IdString('H', backHalfedge) << '\n';
				boundaryHalfEdges.push_back(backHalfedge);
			}
		}

		G.forEachHalfedgeAroundVertex(currentVertex, [&G, predResult, &p, &visitedEdges, &boundaryHalfEdges](index_t he)
		{
			InducedSubgraph<GraphT>::dfs(G, he, predResult, p, visitedEdges, boundaryHalfEdges);
		});

	}
};

#include <iostream>

template<class GraphT>
class DifferenceGraph
{
public:
	DifferenceGraph(const GraphT& G, InducedSubgraph<GraphT>&& T)
		: G_(G)
		, T_(std::move(T))
	{
	}
	
	DifferenceGraph(DifferenceGraph&& rhs)
		: G_(rhs.G_)
		, T_(std::move(rhs.T_))
	{
	}
	
	~DifferenceGraph()
	{
	}

	//for PrintableGraph
	template<typename FuncT>
	FuncT forEachRegion(FuncT f) const
	{
		G_.forEachRegion(f);
		T_.forEachRegion([&](index_t r) { f(r + regionOffset()); });
		return f;
	}

	//for PrintableGraph
	const Vector2C site(index_t r) const
	{
		const auto dr = regionOffset();
		return (r < dr)? G_.site(r) : T_.site(r - dr);
	}

	template<typename FuncT>
	FuncT forEachVertexAroundRegion(index_t r, FuncT f) const
	{
		//when r is of sub-graph
		const auto n = T_.cycleLength();
		if (r == regionOffset())
		{
			const auto dv = vertexOffset();
			for (index_t v = 0; v < n; ++v)
				f(v + dv);
			return f;
		}

		//when r is boolean operated region of super-graph <=> r is left-region of each boundary-halfedge
		for (index_t v = 0; v < n; ++v)
		{
			const auto heFirst = T_.boundaryHalfedge(v);
			if (r == G_.leftRegionOfHalfedge(heFirst))
			{
				f(v);
				const auto heLast = G_.mateHalfedge(T_.boundaryHalfedge((v+1)%n));
				for (index_t he = heFirst; he != heLast; he = G_.halfedgeCCWAroundRegion(r, he))
				{
					f(G_.endVertexOfHalfedge(he));
				}
				return f;
			}
		}

		//when is r is an unmodified region of super-graph
		return G_.forEachVertexAroundRegion(r, f);
	}
	
	index_t forwardHalfedge(index_t e) const
	{
		return EdgeItem::forwardHalfedge(e);
	}
	
	index_t backwardHalfedge(index_t e) const
	{
		return EdgeItem::backwardHalfedge(e);
	}

	//for PrintableGraph	
	template<typename FuncT>
	FuncT forEachEdge(FuncT f) const
	{
		G_.forEachEdge(f);	//TODO: remove filter
		T_.forEachEdge([&](index_t e) { f(e + edgeOffset()); });
		return f;
	}

	//for PrintableGraph
	index_t startVertexOfEdge(index_t e) const
	{
		const auto de = edgeOffset();
		if (e >= de)
		{
			const auto dv = vertexOffset();
			return T_.startVertexOfEdge(e - de) + dv;
		}
		else
		{
			for (index_t v = 0; v < T_.cycleLength(); ++v)
			{
				const auto heSuper = T_.boundaryHalfedge(v);
				const auto e2 = EdgeItem::parentEdge(heSuper);
				if (e == e2 && G_.forwardHalfedge(e) == heSuper)
				{
					const auto dv = vertexOffset();
					return v + dv;
				}
			}
			return G_.startVertexOfEdge(e);
		}
	}
		
	index_t endVertexOfHalfedge(index_t he) const
	{
		const auto e = parentEdge(he);
		return (he == EdgeItem::forwardHalfedge(e))? endVertexOfEdge(e) : startVertexOfEdge(e);
		
	}

	//for PrintableGraph
	index_t endVertexOfEdge(index_t e) const
	{
		const auto de = edgeOffset();
		const auto dv = vertexOffset();
		if (e >= de)
		{
			return T_.endVertexOfEdge(e - de) + dv;
		}
		else
		{
			for (index_t v = 0; v < T_.cycleLength(); ++v)
			{
				const auto heSuper = T_.boundaryHalfedge(v);
				const auto e2 = EdgeItem::parentEdge(heSuper);
				if (e == e2 && G_.backwardHalfedge(e) == heSuper)
				{
					const auto dv = vertexOffset();
					return v + dv;
				}
			}
			return G_.endVertexOfEdge(e);
		}
	}

	index_t parentEdge(index_t he) const
	{
		return EdgeItem::parentEdge(he);
	}

	index_t leftRegionOfHalfedge(index_t he) const
	{
		const auto de = edgeOffset();
		const auto e = parentEdge(he);
		return (e < de)? G_.leftRegionOfHalfedge(he) : (0 + regionOffset());
	}

	index_t regionAroundVertex(index_t v, index_t i) const
	{
		const auto dv = vertexOffset();
		if (v < dv)
		{
			return G_.regionAroundVertex(v, i);
		}
		else
		{
			switch (i)
			{
			case 0:
				return 0 + regionOffset();
			case 1:
				return G_.rightRegionOfHalfedge(T_.boundaryHalfedge(v));
			case 2:
				return G_.leftRegionOfHalfedge(T_.boundaryHalfedge(v));
			default:
				assert(0);
				return NOTHING;
			}
		}
	}

	template<typename FuncT>
	FuncT forEachHalfedgeAroundVertex(index_t v, FuncT f) const
	{
		const auto dv = vertexOffset();
		return (v < dv)? forEachHalfedgeAroundVertexSuper(v, f) : forEachHalfedgeAroundVertexSub(v - dv, f);
	}

	//for PrintableGraph
	const Circumcircle& circumcircle(index_t v) const
	{
		const auto dv = vertexOffset();
		return (v < dv)? G_.circumcircle(v) : T_.circumcircle(v - dv);
	}

	//for PrintableGraph
	const Line2& bisector(index_t e) const
	{
		const auto de = edgeOffset();
		return (e < de)? G_.bisector(e) : T_.bisector(e - de);
	}
	
	index_t maxIdOfRegion() const
	{
		return G_.maxIdOfRegion() + 1;
	}
	
	index_t maxIdOfEdge() const
	{
		return G_.maxIdOfEdge() + T_.cycleLength();
	}

	index_t maxIdOfVertex() const
	{
		return G_.maxIdOfVertex() + T_.cycleLength();
	}

private:
	index_t edgeOffset() const
	{
		return G_.maxIdOfEdge();
	}
	
	index_t vertexOffset() const
	{
		return G_.maxIdOfVertex();
	}

	index_t regionOffset() const
	{
		return G_.maxIdOfRegion();
	}

	template<typename FuncT>
	FuncT forEachHalfedgeAroundVertexSuper(index_t v, FuncT f) const
	{
		G_.forEachHalfedgeAroundVertex(v, [this, &f](index_t he)
		{
			for (index_t v = 0; v < T_.cycleLength(); ++v)
			{
				if (he == T_.boundaryHalfedge(v))
				{
					f(v);	//boundaryHalfedge of v = v
					return;
				}
			}
			f(he);
		});
		return f;
	}

	template<typename FuncT>
	FuncT forEachHalfedgeAroundVertexSub(index_t v, FuncT f) const
	{
		index_t hes[3];
		T_.halfedgesAroundVertex(v, hes);
		f(hes[0]);
		f(hes[1]);
		f(hes[2]);
		return f;
	}

private:
	const GraphT& G_;
	InducedSubgraph<GraphT> T_;
};
