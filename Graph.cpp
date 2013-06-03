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


//�ˉe��ԏ�Ŏl�ʑ̏�̈ʑ��\�����쐬����B
RootGraph::RootGraph(const Vector2& p0, const Vector2& p1, const Vector2& p2)
{
	const Vector2 sites[] =
	{
		p0,
		p1,
		p2
	};

	//WARNING: ��_�������v���ɕ���ł��邱�Ƃ�����
	//TODO: �����������v���łȂ���Δ��]������

	//������3�_�ɑ΂���{���m�C�} Phi_2 �̍\�z
	//�{���m�C�̈�F�������_
	RInf_ = RegionItemC(EdgeItem::backwardHalfedge(3));

	//�{���m�C�̈�
	for (index_t i = 0; i < 3; ++i)
	{
		R_[i] = RegionItem(EdgeItem::backwardHalfedge(i), sites[i]);
	}
	
	//�{���m�C��
	for (index_t i = 0; i < 3; ++i)
	{
		const auto left  = i;
		const auto right = (i+1) % 3;
		E[i] = EdgeItem(0, 1+i, 1+i, 1+ (i+1)%3, ::bisector(sites[left], sites[right]));
	}

	//�{���m�C�ӁF�������_
	const auto w_equal_0 = Line2(Vector2H(0, 0, 1));
	for (index_t i = 0; i < 3; ++i)
	{
		E[3+i] = EdgeItem(1+(i+2)%3, 1+i, 1+i, 0, w_equal_0);
	}

	//�{���m�C���_
	V[0] = VertexItem(
		forwardHalfedge(0),
		forwardHalfedge(1),
		forwardHalfedge(2),
		::circumcircle(Vector2C(sites[0]), Vector2C(sites[1]), Vector2C(sites[2]))
	);

	//�{���m�C���_�F�������_
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
