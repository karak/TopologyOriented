// Voronoi.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"
#include "IdString.h"
#include "Graph.h"
#include "io.h"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <malloc.h>	//_alloca
#include <vector>
#include <queue>


index_t findNearestSite(const Vector2 p[], index_t l)
{
	//とりあえず素朴な実装
	if (l == 0)
	{
		assert(!"input number larger than 0");
		return 0;
	}
	const Vector2& q = p[l];
	index_t result;
	float* d2 = new float[l];
	//float* d2 = static_cast<float*>(_alloca(l * sizeof (float)));
	{
		std::transform(p, p + l, d2, [&q](const Vector2& p_i)
		{
			return squaredDistance(p_i, q);
		});
		result = std::distance(d2, std::min_element(d2, d2 + l));
	}
	//_freea(d2);
	delete[] d2;
	return result;
}


//外接円に対する二乗距離（射影空間上では原点を加えた外接球）
//ボロノイ図の位相に関わる述語を構成する
template<typename GraphT>
float H(const GraphT& G, index_t v, const Vector2& p_l)
{
	const auto C = G.circumcircle(v);
	const auto region = G.regionAroundVertex(v, 0);
	std::clog << "H():\n";
	std::clog << "q_ijk = "<< C.center() <<'\n';
	return C.eval(Vector2C(p_l));
}

template<typename GraphT>
index_t findNearestCircum(const GraphT& G, index_t site, const Vector2& p)
{
	index_t v = NOTHING;
	float value = std::numeric_limits<float>::max();

	G.forEachVertexAroundRegion(
		site + 1,	//TODO: siteToRegion(site)
		[&G, &p, &v, &value](index_t v2)
		{
			auto value2 = H(G, v2, p);
			if (value2 < value)
				v = v2;
		}
	);
	return v;
}

#include <sstream>
#include <iomanip>

template<class GraphT>
InducedSubgraph<GraphT> embedRegion(const GraphT& G, const Vector2 sites[], index_t l)
{
	//新規母点 p_3 の追加
	const auto& p_l = sites[l];
	
	//既存の母点中で新規母点に対して最近のものを取得する
	const auto m = findNearestSite(sites, l);

	//v_m 内のボロノイ領域の境界上でもっともHが小さいものを探す
	index_t v_site_m = findNearestCircum(G, m, p_l);
	
	//新規母点のボロノイ領域に含まれる頂点集合を計算し、
	//その頂点集合によって誘導される部分グラフを埋め込む
	InducedSubgraph<GraphT> T(
		G,
		p_l,
		v_site_m,
		[&G, &p_l](index_t v) {
			return H(G, v, p_l) > 0; //except 0 for infinite point;
		}
	);
	std::stringstream ss;
	ss << "Voronoi-sub-" << std::setfill('0') << std::setw(6) << l << ".ps";
	WriteVoronoiToPsFile(ss.str().c_str(), T);
	return std::move(T);
}


template<class GraphT>
DifferenceGraph<GraphT> updatedVoronoiDiagram(const GraphT& G, const Vector2 sites[], index_t l)
{
	//埋め込まれた部分グラフを除去したグラフを構築し、新たなボロノイ図 Phi_3 を得る
	InducedSubgraph<GraphT> T1 = embedRegion(G, sites, l);

	//BUG! boundaryHalfedgeのbisector/Circumcenterがおかしい（頂点座標が全て原点になっている）
	DifferenceGraph<GraphT> G2(G, std::move(T1));
	return std::move(G2);
}



int _tmain(int argc, _TCHAR* argv[])
{
	Vector2 sites[] = {
		//WARNING: 最初の3つの母点は現在反時計回り前提
		Vector2(100, 500),
		Vector2(400, 500),
		Vector2(250, 500 + 150 * std::sqrtf(3)),
		//BUG: 4つ目元の三角形の内側以外はNG!
		//=>無限遠辺に対するHがおかしい。原点とのBISECTOR取らないように
		//  Vector3HのregionAroundVertexを変更し、Vector2 oneEuclidRegionAroundVertex とする。
		//Vector2(200, 550)
		//Vector2(350, 540)
		Vector2(220, 450)
	};

	RootGraph G0(sites[0], sites[1], sites[2]);
	std::size_t l = 3;		//処理済の母店数
	auto G1 = updatedVoronoiDiagram(G0, sites, l++);
	//BUG: 実行するとG1が壊れる！ moveコンストラクタの扱いのミス?
	//auto G2 = updatedVoronoiDiagram(G1, sites, l++);

	//出力
	WriteVoronoiToPsFile("voronoi0.ps", G0);
	WriteVoronoiToPsFile("voronoi1.ps", G1);
	//WriteVoronoiToPsFile("voronoi2.ps", G2);

	return 0;
}
