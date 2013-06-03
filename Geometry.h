#pragma once

#include "Vector2.h"
#include <iosfwd>


class Line2
{
public:
	Line2() : equation_()
	{
	}
	
	Line2(float a, float b, float c) : equation_(a, b, c)
	{
	};

	explicit Line2(const Vector2H& equation) : equation_(equation)
	{
	};
	
	const Vector2H& equation() const
	{
		return equation_;
	}

	float eval(const Vector2& p) const
	{
		return dotProduct(equation_, p);
	}

	float eval(const Vector2H& p) const
	{
		return dotProduct(equation_, p);
	}

private:
	Vector2H equation_;
};

inline std::ostream& operator<<(std::ostream& os, const Line2& v)
{
	float a = v.equation().X;
	float b = v.equation().Y;
	float c = v.equation().w;
	return os << '(' << a << ' ' << b << ' ' << c << ')';
	//return os << '(' << a << ") * x + (" << b << ") * y + (" << c << ") = 0";
}

inline const Line2 line(const Vector2H& a, const Vector2H& b)
{
	return Line2(wedgeProduct(a, b));
}

inline const Line2 bisector(const Vector2& left, const Vector2& right)
{
	const Vector2H eq((right.x - left.x), (right.y - left.y), (left.norm2()- right.norm2()) / 2);
	return Line2(eq);
}

inline const float squaredDistance(const Vector2& a, const Vector2& b)
{
	return (b - a).norm2();
}

inline const float squaredDistance(const Vector2H& a, const Vector2H& b)
{
	return (b - a).norm2();
}

/*
外接円

ユークリッド平面上の3点の場合の行列（本来のconformal modelから導出すると二次の項に1/2がかかる）:
M = | 1 xi yi xi^2+yi^2 |
	| 1 xj yj xj^2+yj^2 |
	| 1 xk yk xk^2+yk^2 |
	| 1 x  y  x ^2+y ^2 |
うち1点が無限遠点の場合は直線に縮退する
M = | 0  0  0         1 |
	| 1 xj yj xj^2+yj^2 |
	| 1 xk yk xk^2+yk^2 |
	| 1 x  y  x ^2+y ^2 |

mii = 余因子展開で行列式を求める
f(x, y) = m44 * (x^2+y^2) + m33 * y + m22 * x + m11 
*/
class Circumcircle
{
public:
	Circumcircle() : m11(), m22(), m33(), m44()
	{
	}

	Circumcircle(float m11, float m22, float m33, float m44) : m11(m11), m22(m22), m33(m33), m44(m44)
	{
	}

	Vector2 center() const
	{
		const float minus_half_inv_m44 = -0.5f / m44;
		return Vector2(minus_half_inv_m44 * m22, minus_half_inv_m44 * m33);
	}

	float eval(const Vector2C& p) const
	{
		return m44 * p.i + m33 * p.Y + m22 * p.X + m11 * p.w;
	}

private:
	float m11;
	float m22;
	float m33;
	float m44;
};

inline float det3x3(float a11, float a12, float a13,
             float a21, float a22, float a23,
             float a31, float a32, float a33)
{
	return a11*a22*a33+a12*23*a31+a13*a21*a32-a11*a23*a32-a12*a21*a33-a13*a22*a31;
}

inline Circumcircle circumcircle(const Vector2C& a, const Vector2C& b, const Vector2C& c)
{
	return Circumcircle(
		det3x3(a.X, a.Y, a.i, b.X, b.Y, b.i, c.X, c.Y, c.i),
		det3x3(a.w, a.Y, a.i, b.w, b.Y, b.i, c.w, c.Y, c.i),
		det3x3(a.w, a.X, a.i, b.w, b.X, b.i, c.w, c.X, c.i),
		det3x3(a.w, a.X, a.Y, b.w, b.X, b.Y, c.w, c.X, c.Y)
	);
}
