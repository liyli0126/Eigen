// vector.cpp
// Source code for (double) vectors
// James Liu, Rachel Cali, Graham Harper, ColoState; 2007/01--2018/02

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "vector.h"
using namespace std;
// using namespace LinLite;


// Vector: Constructor with given size and all elements 0

Vector::Vector(int numElts)
{
  n = numElts;
  p = new double[n];
  for (int i=0; i<n; ++i)  p[i]=0;
}


// Vector: Constructor with given size and data

Vector::Vector(int numElts, double *data)
{
  n = numElts;
  p = new double[n];
  for (int i=0; i<n; ++i)  p[i] = data[i];
}


// Vector: Copy constructor

Vector::Vector(const Vector &v)
{
  n = v.n;
  p = new double[n];
  for (int i=0; i<n; ++i)  p[i] = v.p[i];
}


// Vector: Copy assignment

Vector &Vector::operator=(const Vector &v)
{
  if (this!=&v) {
    delete[] p;
    n = v.n;
    p = new double[n];
    for (int i=0; i<n; ++i)  p[i] = v.p[i];
  }
  return *this;
}


// Vector: l_infty norm

double Vector::l0norm() const
{
  double norm = fabs(p[0]);
  for (int i=1; i<n; ++i) {
    if (norm<fabs(p[i]))  norm = fabs(p[i]);
  }
  return norm;
}


// Vector: l_1 norm

double Vector::l1norm() const
{
  double norm = 0.0;
  for (int i=0; i<n; ++i)  norm += fabs(p[i]);
  return norm;
}


// Vector: l_2 norm

double Vector::l2norm() const
{
  double norm = 0.0;
  for (int i=0; i<n; ++i) {
    double tmp = fabs(p[i]);
    norm += tmp*tmp;
  }
  return sqrt(norm);
}


// Vector: resize to a zero vector of a given size

void Vector::resize(int numElts)
{
  delete[] p;
  n = numElts;
  p = new double[n];
  for (int i=0; i<n; ++i)  p[i] = 0;
  return;
}


// Vector: Unary -

Vector Vector::operator-() const
{
  Vector u(n);
  for (int i=0; i<n; ++i)  u[i] = -p[i];
  return u;
}


// Vector: v += w

Vector &Vector::operator+=(const Vector &w)
{
  if (n!=w.n) {
    std::cout << "Bad vector sizes\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<n; ++i)  p[i] += w.p[i];
  return *this;
}


// Vector: v -= w

Vector &Vector::operator-=(const Vector &w)
{
  if (n!=w.n) {
    std::cout << "Bad vector sizes\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<n; ++i)  p[i] -= w.p[i];
  return *this;
}


// Vector: get a sub-vector
// 1<=i: starting position

Vector Vector::getSubVector(int i, int length) const
{
  if (n < i+length-1){
    std::cout << "Bad vector sizes\n";
    exit(EXIT_FAILURE);
  }
  Vector w(length);
  for (int j=0; j<length; ++j)  w.p[j] = p[(i-1)+j];
  return w;
}


// Vector: replace a sub-vector
// 1<=i: starting position

Vector &Vector::replaceSubVector(int i, const Vector &w)
{
  if (n < i+w.size()-1){
    std::cout << "Bad vector sizes\n";
    exit(EXIT_FAILURE);
  }
  for (int j=0; j<w.size(); ++j)  p[(i-1)+j] = w.p[j];
  return *this;
}


// Vector: set/replace a sub-vector
// 1<=i: starting position

Vector &Vector::setSubVector(int i, const Vector &w)
{
  if (n < i+w.size()-1){
    std::cout << "Bad vector sizes\n";
    exit(EXIT_FAILURE);
  }
  for (int j=0; j<w.size(); ++j)  p[(i-1)+j] = w.p[j];
  return *this;
}


// Vector: add a sub-vector
// 1<=i: starting position

Vector &Vector::addSubVector(int i, const Vector &w)
{
  if (n < i+w.size()-1){
    std::cout << "Bad vector sizes\n";
    exit(EXIT_FAILURE);
  }
  for (int j=0; j<w.size(); ++j)  p[(i-1)+j] += w.p[j];
  return *this;
}


// Vector: save to a data file

void Vector::save2file(char *filename) const
{
  FILE *fp;

  if (NULL==(fp=fopen(filename, "w"))) {
    std::cout << "Open data file failed.  Exit!\n";
    exit(-1);
  }

  fprintf(fp, "%8d\n", n);

  for (int i=0; i<n; ++i)
    fprintf(fp, "%8d  %21.9f\n", i+1, p[i]);

  fclose(fp);
  return;
}


// Vector: Output stream

std::ostream &operator<<(std::ostream &strm, const Vector&v)
{
  for (int i=0; i<v.n; i++) {
     strm << v[i] << "  ";
     if (i%5==4)  strm << "\n";
  }
  return strm;
}


// Vector: v = v1 + v2

Vector operator+(const Vector &v1, const Vector &v2) 
{
  if (v1.n!=v2.n)  {
    std::cout << "Bad vector sizes!\n";
    exit(EXIT_FAILURE);
  }
  Vector v(v1);
  v += v2;
  return v;
}


// Vector: v = v1 - v2
// Notice that the implementation is different than 
// that in the overloaded binary operator +

Vector operator-(const Vector &v1, const Vector &v2) 
{
  if (v1.n!=v2.n) {
    std::cout << "Bad vector sizes!\n";
    exit(EXIT_FAILURE);
  }
  Vector v(v1.n);
  for (int i=0; i<v1.n; ++i)  v[i] = v1[i] - v2[i];
  return v;
}


// Vector: scalar-vector multiplication

Vector operator*(double a, const Vector &v)
{
  Vector u(v.n);
  for (int i=0; i<v.n; ++i)  u[i] = a*v[i];
  return u;
}


// Vector: dot product

double dotProd(const Vector &v1, const Vector &v2)
{
  if (v1.n!=v2.n) {
    std::cout << "Bad vector sizes!\n";
    exit(EXIT_FAILURE);
  }
  double dp = 0;
  for (int i=0; i<v1.n; ++i)  dp += v1[i]*v2[i];
  return dp;
}
