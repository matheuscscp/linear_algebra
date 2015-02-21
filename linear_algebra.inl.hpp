/*
 * linear_algebra.inl.hpp
 *
 *  Created on: Feb 21, 2015
 *      Author: Pimenta
 */

#ifndef LINEAR_ALGEBRA_INL_HPP_
#define LINEAR_ALGEBRA_INL_HPP_

#include <cmath>

namespace linear_algebra {

template <typename T>
struct Vector3 {
  T x, y, z;
  
  Vector3(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) {
    
  }
  
  // arithmetic operators
  inline Vector3 operator+(const Vector3& other) {
    return Vector3(x + other.x, y + other.y, z + other.z);
  }
  inline Vector3 operator-(const Vector3& other) {
    return Vector3(x - other.x, y - other.y, z - other.z);
  }
  inline Vector3 operator-() {
    return Vector3(-x, -y, -z);
  }
  inline Vector3 operator*(T scalar) {
    return Vector3(x*scalar, y*scalar, z*scalar);
  }
  inline Vector3 operator/(T scalar) {
    return Vector3(x/scalar, y/scalar, z/scalar);
  }
  
  // compound assignment
  Vector3& operator+=(const Vector3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
  Vector3& operator-=(const Vector3& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }
  Vector3& operator*=(T scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
  }
  Vector3& operator/=(T scalar) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
  }
  
  // relational operators
  inline bool operator==(const Vector3& other) {
    return x == other.x && y == other.y && z == other.z;
  }
  inline bool operator!=(const Vector3& other) {
    return x != other.x || y != other.y || z != other.z;
  }
  
  // methods
  inline T length() {
    return sqrt(x*x + y*y + z*z);
  }
  inline Vector3 versor() {
    T len = length();
    if (len == 0) {
      return *this;
    }
    return operator/(len);
  }
  inline T dot(const Vector3& other) {
    return x*other.x + y*other.y + z*other.z;
  }
  inline Vector3 cross(const Vector3& v) {
    return Vector3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.x - y*v.x);
  }
  inline T angle(const Vector3& other) {
    T len = length()*other.length();
    if (len == 0) {
      return 0;
    }
    return acos(dot(other)/len);
  }
  inline Vector3 proj(const Vector3& other) { // other on this
    Vector3 v = versor();
    return v*v.dot(other);
  }
  inline Vector3 rej(const Vector3& other) { // other from this
    return other - proj(other);
  }
  inline T scalarproj(const Vector3& other) { // other on this
    return versor().dot(other);
  }
  inline Vector3 rotate(const Vector3& other, T angle) { // other around this
    Vector3 u = versor();
    T a2 = u.x*u.x, ab = u.x*u.y, ac = u.x*u.z;
    T b2 = u.y*u.y, bc = u.y*u.z, c2 = u.z*u.z;
    T sint = sin(angle), cost = cos(angle), _1mcost = 1 - cost;
    T asint = u.x*sint, bsint = u.y*sint, csint = u.z*sint;
    T x = other.x, y = other.y, z = other.z;
    return Vector3(
      x*(a2*_1mcost +  cost) + y*(ab*_1mcost - csint) + z*(ac*_1mcost + bsint),
      x*(ab*_1mcost + csint) + y*(b2*_1mcost +  cost) + z*(bc*_1mcost - asint),
      x*(ac*_1mcost - bsint) + y*(bc*_1mcost + asint) + z*(c2*_1mcost +  cost)
    );
  }
};

// determinant
template <typename T>
T det(const Vector3<T>& a, const Vector3<T>& b) {
  return a.x*b.y - a.y*b.x;
}
template <typename T>
T det(const Vector3<T>& a, const Vector3<T>& b, const Vector3<T>& c) {
  return
      a.x*b.y*c.z + b.x*c.y*a.z + c.x*a.y*b.z
    - a.z*b.y*c.x - b.z*c.y*a.x - c.z*a.y*b.x
  ;
}

// linear system
template <typename T>
bool solvesys(
  const Vector3<T>& a, const Vector3<T>& b,
  const Vector3<T>& c,
  T& x, T& y
) {
  T D = det(a, b);
  if (D == 0) {
    return false;
  }
  x = det(c, b)/D;
  y = det(a, c)/D;
  return true;
}
template <typename T>
bool solvesys(
  const Vector3<T>& a, const Vector3<T>& b, const Vector3<T>& c,
  const Vector3<T>& d,
  T& x, T& y, T& z
) {
  T D = det(a, b, c);
  if (D == 0) {
    return false;
  }
  x = det(d, b, c)/D;
  y = det(a, d, c)/D;
  z = det(a, b, d)/D;
  return true;
}

// basic types
typedef Vector3<float>  Vector3f;
typedef Vector3<double> Vector3d;

} // namespace linear_algebra

#endif /* LINEAR_ALGEBRA_INL_HPP_ */
