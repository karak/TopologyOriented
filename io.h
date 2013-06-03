#pragma once

#include "Vector2.h"
#include "PostscriptWriter.h"
#include <limits>
#include <cmath>


inline bool segmentOnEuclid(const Vector2H& s, const Vector2H& e, Vector2& sTrimmed, Vector2& eTrimmed)
{
	static const float eps = std::numeric_limits<float>::epsilon();
	static const float t = 0.8f;
	static const float one_minus_t = 1 - t;
	const bool isStartOnEuclid = std::abs(s.w) > eps;
	const bool isEndOnEuclid = std::abs(e.w) > eps;
	if (isStartOnEuclid)
	{
		sTrimmed = s.asEuclidPoint();
	}
	else
	{
		//トリム方法はとりあえずスケール揃えて間の点を適当に取るだけ
		sTrimmed = (one_minus_t * e + t * e.w * s).asEuclidPoint();
	}
	if (isEndOnEuclid)
	{
		eTrimmed = e.asEuclidPoint();
	}
	else
	{
		eTrimmed = (one_minus_t * s + t * s.w * e).asEuclidPoint();
	}

	return isStartOnEuclid || isEndOnEuclid;
};

template<class GraphT>
void WriteVoronoiToPsFile(const char* fileName, const GraphT& G)
{
	PostscriptWriter writer(fileName);
	writer.setFontBatch("Helvetica", 12);
	writer.comment("sites");
	G.forEachRegion([&G, &writer](index_t i)
	{
		const auto s = G.site(i);
		if(s.w < std::numeric_limits<float>::epsilon())
		{
			return;
		}

		const auto p = s.asEuclidPoint();
		writer
			.moveTo(p + Vector2(-5, -5))
			.lineTo(p + Vector2( 5, -5))
			.lineTo(p + Vector2( 5, 5))
			.lineTo(p + Vector2(-5, 5))
			.fill();
		IdString id('S', i);
		writer
			.moveTo(p + Vector2(10, 10))
			.show(id.c_str());
	});


	writer.comment("edges\n\texcept infinite");
	G.forEachEdge([&G, &writer](index_t e)
	{
		auto vs = G.startVertexOfEdge(e);
		auto ve = G.endVertexOfEdge(e);
#if 0
		std::clog
			<< IdString('E', e) << G.bisector(e)
			<< '='
			<< IdString('V', vs) << G.circumcenter(vs)
			<< ",\t"
			<< IdString('V', ve) << G.circumcenter(ve)
			<< '\n';
#endif
		Vector2 ps, pe;
		if (segmentOnEuclid(G.circumcircle(vs).center(), G.circumcircle(ve).center(), ps, pe))
		{
			writer
				.moveTo(ps)
				.lineTo(pe)
				.stroke();
			writer
				.moveTo((0.85f * ps + 0.15f * pe))
				.show(IdString('E', e).c_str());
		}
	});
}
