#pragma once

#include <ostream>


class Vector2
{
public:
	Vector2() : x(), y()
	{
	}

	Vector2(float x, float y) : x(x), y(y)
	{
	}

	///squared norm
	float norm2() const
	{
		return x * x + y * y;
	}
	
	Vector2& operator+=(const Vector2& rhs)
	{
		x += rhs.x;
		y += rhs.y;
		return *this;
	}
	
	Vector2& operator-=(const Vector2& rhs)
	{
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	}

	Vector2& operator*=(float rhs)
	{
		x *= rhs;
		y *= rhs;
		return *this;
	}

	union
	{
		struct
		{
			float x;
			float y;
		};
		float v[2];
	};
};

inline const Vector2 operator+(const Vector2& lhs, const Vector2& rhs)
{
	Vector2 result(lhs);
	result += rhs;
	return (result);
}

inline const Vector2 operator-(const Vector2& lhs, const Vector2& rhs)
{
	Vector2 result(lhs);
	result -= rhs;
	return (result);
}

inline const Vector2 operator*(const Vector2& lhs, const float& rhs)
{
	Vector2 result(lhs);
	result *= rhs;
	return (result);
}

inline const Vector2 operator*(const float& lhs, const Vector2& rhs)
{
	return (rhs * lhs);
}

inline std::ostream& operator<<(std::ostream& os, const Vector2& v)
{
	os << '[' << v.x << ' ' << v.y << ']';
	return os;
}

/** 3D homogenius */
class Vector2H
{
public:
	Vector2H() : X(), Y(), w()
	{
	}

	Vector2H(float X, float Y, float w) : X(X), Y(Y), w(w)
	{
	}

	Vector2H(const Vector2& rhs) : X(rhs.x), Y(rhs.y), w(1)
	{
	}

	Vector2H& operator*=(float k)
	{
		X *= k;
		Y *= k;
		w *= k;
		return *this;
	}

	Vector2H& operator+=(const Vector2H& rhs)
	{
		X += rhs.X;
		Y += rhs.Y;
		w += rhs.w;
		return *this;
	}

	Vector2H& operator-=(const Vector2H& rhs)
	{
		X -= rhs.X;
		Y -= rhs.Y;
		w -= rhs.w;
		return *this;
	}

	const Vector2H operator-() const
	{
		return Vector2H(-X,-Y,-w);
	}
	
	///squared norm
	float norm2() const
	{
		return X * X + Y * Y + w * w;
	}

	const Vector2 asEuclidPoint() const
	{
		return Vector2(X/w, Y/w);
	}

	union
	{
		struct
		{
			float X;
			float Y;
			float w;
		};
		float v[3];
	};

	//make omega to be positive
	const Vector2H normalizedOmega() const
	{
		return (w < 0)? -(*this) : *this;
	}
};

inline std::ostream& operator<<(std::ostream& os, const Vector2H& v)
{
	os << '[' << v.X << ' ' << v.Y << ' ' << v.w << ']';
	return os;
}

inline const Vector2H wedgeProduct(const Vector2& a, const Vector2& b)
{
	return Vector2H(a.y - b.y, b.x - a.x, a.x * b.y - a.y * b.x);
}

inline const Vector2H wedgeProduct(const Vector2H& a, const Vector2H& b)
{
	return Vector2H(a.Y * b.w - a.w * b.Y, a.w * b.X - a.X * b.w, a.X * b.Y - a.Y * b.X);
}

inline float dotProduct(const Vector2H& a, const Vector2& b)
{
	return (a.X * b.x + a.Y * b.y + a.w);
}

inline float dotProduct(const Vector2H& a, const Vector2H& b)
{
	return (a.X * b.X + a.Y * b.Y + a.w * b.w);
}

inline const Vector2H operator+(const Vector2H& lhs, const Vector2H& rhs)
{
	Vector2H result(lhs);
	result += rhs;
	return result;
}

inline const Vector2H operator-(const Vector2H& lhs, const Vector2H& rhs)
{
	Vector2H result(lhs);
	result -= rhs;
	return result;
}

inline const Vector2H operator*(const Vector2H& lhs, float rhs)
{
	Vector2H result(lhs);
	result *= rhs;
	return result;
}

inline const Vector2H operator*(float lhs, const Vector2H& rhs)
{
	return rhs * lhs;
}

#include <cassert>

/** 4D conformal */
class Vector2C
{
public:
	float w;	//e
	float X;	//e1
	float Y;	//e2
	float i;	//e12

public:
	static Vector2C origin() { return Vector2C(1, 0, 0, 0); }
	static Vector2C infinity() { return Vector2C(0, 0, 0, 1); }
	
	Vector2C() : w(1), X(), Y(), i()
	{
	}

	explicit Vector2C(const Vector2& rhs) : w(1), X(rhs.x), Y(rhs.y), i(rhs.x*rhs.x+rhs.y*rhs.y)
	{
	}

	Vector2C(float w, float X, float Y, float v): w(w), X(X), Y(Y), i(i)
	{
	}

	Vector2 asEuclidPoint() const
	{
		assert(radiusSquared() <= std::numeric_limits<float>::epsilon());
		return Vector2(X/w, Y/w);
	}
private:
	float radiusSquared() const
	{
		//if negative meaning imaginary!
		return std::abs(X * X + Y * Y - i*i) / (w * w);
	}
};

/*
inline Vector2C geometricProduct(const Vector2C& lhs, const Vector2C& rhs)
{
}
*/

