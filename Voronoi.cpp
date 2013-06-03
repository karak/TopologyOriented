// Voronoi.cpp : �R���\�[�� �A�v���P�[�V�����̃G���g�� �|�C���g���`���܂��B
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
	//�Ƃ肠�����f�p�Ȏ���
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


//�O�ډ~�ɑ΂����拗���i�ˉe��ԏ�ł͌��_���������O�ڋ��j
//�{���m�C�}�̈ʑ��Ɋւ��q����\������
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
	//�V�K��_ p_3 �̒ǉ�
	const auto& p_l = sites[l];
	
	//�����̕�_���ŐV�K��_�ɑ΂��čŋ߂̂��̂��擾����
	const auto m = findNearestSite(sites, l);

	//v_m ���̃{���m�C�̈�̋��E��ł����Ƃ�H�����������̂�T��
	index_t v_site_m = findNearestCircum(G, m, p_l);
	
	//�V�K��_�̃{���m�C�̈�Ɋ܂܂�钸�_�W�����v�Z���A
	//���̒��_�W���ɂ���ėU������镔���O���t�𖄂ߍ���
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
	//���ߍ��܂ꂽ�����O���t�����������O���t���\�z���A�V���ȃ{���m�C�} Phi_3 �𓾂�
	InducedSubgraph<GraphT> T1 = embedRegion(G, sites, l);

	//BUG! boundaryHalfedge��bisector/Circumcenter�����������i���_���W���S�Č��_�ɂȂ��Ă���j
	DifferenceGraph<GraphT> G2(G, std::move(T1));
	return std::move(G2);
}



int _tmain(int argc, _TCHAR* argv[])
{
	Vector2 sites[] = {
		//WARNING: �ŏ���3�̕�_�͌��ݔ����v���O��
		Vector2(100, 500),
		Vector2(400, 500),
		Vector2(250, 500 + 150 * std::sqrtf(3)),
		//BUG: 4�ڌ��̎O�p�`�̓����ȊO��NG!
		//=>�������ӂɑ΂���H�����������B���_�Ƃ�BISECTOR���Ȃ��悤��
		//  Vector3H��regionAroundVertex��ύX���AVector2 oneEuclidRegionAroundVertex �Ƃ���B
		//Vector2(200, 550)
		//Vector2(350, 540)
		Vector2(220, 450)
	};

	RootGraph G0(sites[0], sites[1], sites[2]);
	std::size_t l = 3;		//�����ς̕�X��
	auto G1 = updatedVoronoiDiagram(G0, sites, l++);
	//BUG: ���s�����G1������I move�R���X�g���N�^�̈����̃~�X?
	//auto G2 = updatedVoronoiDiagram(G1, sites, l++);

	//�o��
	WriteVoronoiToPsFile("voronoi0.ps", G0);
	WriteVoronoiToPsFile("voronoi1.ps", G1);
	//WriteVoronoiToPsFile("voronoi2.ps", G2);

	return 0;
}
