#include "stdafx.h"
#include "Graph.h"
#include <algorithm>
#include <cassert>


int VertexItem::indexOf(index_t he) const
{
	const auto it = std::find(halfedges, halfedges + DEGREE, he);
	assert(it != halfedges + DEGREE);
	return std::distance(halfedges, it);
}


//射影空間上で四面体上の位相構造を作成する。
RootGraph::RootGraph(const Vector2& p0, const Vector2& p1, const Vector2& p2)
{
	const Vector2 sites[] =
	{
		p0,
		p1,
		p2
	};

	//WARNING: 母点が反時計回りに並んでいることを仮定
	//TODO: 順序が反時計回りでなければ反転させる

	//初期の3点に対するボロノイ図 Phi_2 の構築
	//ボロノイ領域：無限遠点
	RInf_ = RegionItemC(EdgeItem::backwardHalfedge(3));

	//ボロノイ領域
	for (index_t i = 0; i < 3; ++i)
	{
		R_[i] = RegionItem(EdgeItem::backwardHalfedge(i), sites[i]);
	}
	
	//ボロノイ辺
	for (index_t i = 0; i < 3; ++i)
	{
		const auto left  = i;
		const auto right = (i+1) % 3;
		E[i] = EdgeItem(0, 1+i, 1+i, 1+ (i+1)%3, ::bisector(sites[left], sites[right]));
	}

	//ボロノイ辺：無限遠点
	const auto w_equal_0 = Line2(Vector2H(0, 0, 1));
	for (index_t i = 0; i < 3; ++i)
	{
		E[3+i] = EdgeItem(1+(i+2)%3, 1+i, 1+i, 0, w_equal_0);
	}

	//ボロノイ頂点
	V[0] = VertexItem(
		forwardHalfedge(0),
		forwardHalfedge(1),
		forwardHalfedge(2),
		::circumcircle(Vector2C(sites[0]), Vector2C(sites[1]), Vector2C(sites[2]))
	);

	//ボロノイ頂点：無限遠点
	for (index_t i = 0; i < 3; ++i)
	{
		V[1+i] = VertexItem(
			backwardHalfedge(i),
			backwardHalfedge(3 + i),
			forwardHalfedge(3 + (i+1)%3),
			::circumcircle(Vector2C(sites[i]), Vector2C(sites[(i+3)%3]), Vector2C::infinity())
		); 
	}
}

RootGraph::RootGraph(const RootGraph& rhs)
	: RInf_(rhs.RInf_), R_(rhs.R_), E(rhs.E), V(rhs.V)
{
}
