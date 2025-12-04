// vector.h
// Classes for (double) vectors
// James Liu, Rachel Cali, Graham Harper, ColoState; 2007/01--2018/02

#ifndef Vector_H
#define Vector_H

#include <iostream>

// namespace LinLite{
  class Vector {
  private: 
    int n;
    double *p;
  public: 
    Vector() {n=0; p=0;}  // Default empty constructor
    Vector(int numElts);  // Constructor for a zero vector
    Vector(int numElts, double *data);  // Constructor with given data
    Vector(const Vector &v);  // Copy constructor
    Vector &operator=(const Vector &v);  // Copy assignment
    ~Vector() {delete[] p;}  // Destructor
    // Copy constructor is implicitly used by any function returnig a vector
    // Destructor is also implicitly called by those functions
    int size() const {return n;}
    double l0norm() const;  // l_infty norm
    double l1norm() const;  // l_1 norm
    double l2norm() const;  // l_2 norm
    double &operator[](int i) const {return p[i];}  // 0<=i<n
    double &operator()(int i) const {return p[i-1];}  // 1<=i<=n
    void resize(int numElts);  // resize to a zero vector of a given size
    Vector operator-() const;  // Unary -
    Vector &operator+=(const Vector &w);  // v += w
    Vector &operator-=(const Vector &w);  // v -= w
    Vector getSubVector(int i, int length) const;  // 1<=i: starting position
    Vector &replaceSubVector(int i, const Vector &w);  // 1<=i: starting pos.
    Vector &setSubVector(int i, const Vector &w);  // 1<=i: starting position
    Vector &addSubVector(int i, const Vector &w);  // 1<=i: starting position
    void save2file(char *filename) const;  // save to a data file
  friend std::ostream &operator<<(std::ostream&strm, const Vector &v);
  friend Vector operator+(const Vector &v1, const Vector &v2);  // v = v1 + v2;
  friend Vector operator-(const Vector &v1, const Vector &v2);  // v = v1 + v2;
  friend Vector operator*(double a, const Vector &v);  // Scalar-vector multipl.
  friend double dotProd(const Vector &v1, const Vector &v2);  // Dot product
  };
// }

#endif  // Vector_H
