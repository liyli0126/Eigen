// LinSys.cpp
// Solvers for dense/full and sparse linear systems
// James Liu, Rachel Cali, Graham Harper, ColoState; 2007/01--2018/02

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "LinSys.h"
#include "matrix.h"
#include "vector.h"

using namespace std;


// Solving a small full linear system by the 
// Gaussian Elimination with Partial Pivoting (GEPP) 

Vector slvFullLinSysGEPP(const FullMatrix &A, const Vector &b) 
{
  int i, im, j, k, m, n;
  double cmax, lier, tmp;
  double *p, *q;

  m = A.rowSize();
  n = A.colSize();

  if (m!=n) {
    std::cout << "Not a square matrix! EXit!\n";
    exit(EXIT_FAILURE);
  }

  p = new double[n*n];
  q = new double[n];

  for (i=1; i<=n; ++i) 
    for (j=1; j<=n; ++j) 
      p[(i-1)*n+(j-1)] = A(i,j);

  for (i=1; i<=n; ++i)  q[i-1] = b(i);

  for (j=1; j<=(n-1); ++j) {
    cmax = fabs(p[(j-1)*n+(j-1)]);
    im = j;

    for (i=j+1; i<=n; ++i) {
      if (fabs(p[(i-1)*n+(j-1)])>cmax)  im = i;
    }

    if (cmax==0) {
      std::cout << "Singular matrix! Exit!\n";
      std::cout << "Column " << j << "\n";
      exit(EXIT_FAILURE);
    }

    for (k=j; k<=n; ++k) {
      tmp = p[(j-1)*n+(k-1)];
      p[(j-1)*n+(k-1)] = p[(im-1)*n+(k-1)];
      p[(im-1)*n+(k-1)] = tmp;
    }

    tmp = q[j-1];
    q[j-1] = q[im-1];
    q[im-1] = q[j-1];

    for (i=j+1; i<=n; ++i) {
      lier = p[(i-1)*n+(j-1)]/p[(j-1)*n+(j-1)];
      q[i-1] -= lier*q[j-1];
      for (k=j+1; k<=n; ++k) 
        p[(i-1)*n+(k-1)] -= lier*p[(j-1)*n+(k-1)];
    }
  }

  for (i=1; i<=n; ++i) {
    if (fabs(p[(i-1)*n+(i-1)])==0) {
      std::cout << "Singular matrix! Exit!\n";
      std::cout << "trouble at diagonal" << i << "\n";
      exit(EXIT_FAILURE);
    }
  }

  Vector x(n);

  x(n) = q[n-1]/p[n*n-1];
  for (i=(n-1); i>=1; --i) {
    x(i) = q[i-1];
    for (k=i+1; k<=n; ++k) {
      x(i) -= p[(i-1)*n+(k-1)]*x(k);
    }
    x(i) /= p[(i-1)*n+(i-1)];
  }

  delete[] p;
  delete[] q;

  return x;
}


// Solving a small full lower triangular system

Vector slvFullLowerTrigSys(const FullMatrix &L, const Vector &b) 
{
  int i, j, m, n;
  double sum;

  m = L.rowSize();
  n = L.colSize();

  if (m!=n) {
    std::cout << "Not a square matrix!\n";
    exit(EXIT_FAILURE);
  }

  for (i=1; i<=(n-1); ++i) {
    for (j=i+1; j<=n; ++j) {
      if (L(i,j)!=0) {
        std::cout << "Not a lower triangular matrix!\n";
        exit(EXIT_FAILURE);
      }
    }
  }

  for (i=1; i<=n; ++i) {
    if (L(i,i)==0.0) {
      std::cout << "Singular lower triangular matrix!\n";
      exit(EXIT_FAILURE);
    }
  }

  Vector x(n);

  x(1) = b(1)/L(1,1);
  for (i=2; i<=n; ++i) {
    sum = b(i);
    for (j=1; j<=(i-1); ++j)  sum -= L(i,j)*x(j);
    x(i) = sum/L(i,i);
  }

  return x;
}


// Solving a small full upper triangular system

Vector slvFullUpperTrigSys(const FullMatrix &U, const Vector &b) 
{
  int i, j, m, n;
  double sum;

  m = U.rowSize();
  n = U.colSize();

  if (m!=n) {
    std::cout << "Not a square matrix!\n";
    exit(EXIT_FAILURE);
  }

  for (i=2; i<=n; ++i) {
    for (j=1; j<=(i-1); ++j) {
      if (U(i,j)!=0) {
        std::cout << "Not an upper triangular matrix!\n";
        exit(EXIT_FAILURE);
      }
    }
  }

  for (i=1; i<=n; ++i) {
    if (U(i,i)==0.0) {
      std::cout << "Singular upper triangular matrix!\n";
      exit(EXIT_FAILURE);
    }
  }

  Vector x(n);

  x(n) = b(n)/U(n,n);
  for (i=(n-1); i>=1; --i) {
    sum = b(i);
    for (j=i+1; j<=n; ++j)  sum -= U(i,j)*x(j);
    x(i) = sum/U(i,i);
  }

  return x;
}


// Solving a small full SPD system by the Cholesky factorization

Vector slvFullSpdSysCholesky(const FullMatrix &A, const Vector &b) 
{
  int i, j, k;
  double sum1, sum2;

  int m = A.rowSize();
  int n = A.colSize();

  if (m!=n) {
    std::cout << "Not a square matrix! Exit!\n";
    exit(EXIT_FAILURE);
  }

  FullMatrix L(n,n);  // The lower triangular ChoLesky
  FullMatrix U(n,n);  // U = L' (transpose)

  for (j=1; j<=n; ++j) {
    sum1 = 0.0;
    for (k=1; k<=(j-1); ++k)  sum1 += L(j,k)*L(j,k);
    L(j,j) = sqrt(A(j,j)-sum1);
    for (i=j+1; i<=n; ++i) {
      sum2 = 0.0;
      for (k=1; k<=(j-1); ++k)  sum2 += L(i,k)*L(j,k);
      L(i,j) = (A(i,j)-sum2)/L(j,j);
    }
  }

  U = transpose(L);
  Vector x(n), y(n);

  y = slvFullLowerTrigSys(L, b);
  x = slvFullUpperTrigSys(U, y);

  return x;
}


// Solving a small full SPD system by the Conjugate Gradient method
// x brings in an initial guess and returns the solution 

void slvFullSpdSysCG(Vector &x, const FullMatrix &A, const Vector &b, 
  int &itr, int maxitr, double threshold, double tol) 
{
  int n;
  double alpha, b2, mu, nu, r2;

  b2 = b.l2norm();

  n = x.size();
  Vector r(n);
  Vector u(n);
  Vector v(n);

  itr = 0;
  r = b-A*x;
  v = r;
  mu = dotProd(r,r);

  while (itr<maxitr) {
    if (v.l2norm()<threshold)  break;
    u = A*v;
    alpha = mu/dotProd(u,v);
    x += alpha*v;
    r -= alpha*u;
    nu = dotProd(r,r);
    r2 = sqrt(nu);
    if (r2<tol*b2)  break;
    v = r + (nu/mu)*v;
    mu = nu;
    itr++;
  }

  return;
}


// Solving a diagonal system directly

Vector slvDiagSys(const DiagMatrix &A, const Vector &b) 
{
  int n = b.size();
  Vector x(n);
  for (int i=1; i<=n; ++i)  x(i) = b(i)/A.entry(i);
  return x;
}


// Solving a nonsymmetric sparse linear system by 
// BiCGStab (Preconditioned Bi-Conjugate Gradient Stablized)
// Courtesy of Victor Ginting 
// cf: p.27 of the SIAM Templates book
// The return value indicates 
//   convergence within maxitr iterations (0), 
//   no convergence within maxitr iterations (1) 
// Upon successful return: 
//       x -- approximate solution to Ax = b
//  maxitr -- the number of iterations performed 
//            before the tolerance was reached
//     tol -- the residual after the final iteration

int slvSpaLinSysBiCGStab(Vector &x, 
  const SparseMatrix &A, const Vector &b, const DiagMatrix &B, 
  int &maxitr, double &tol, double atol, int printit)
{
  double alpha, beta, omega, res, rho1, rho2;

  int n = A.colSize();
  Vector p(n), phat(n), r(n), rtilde(n), s(n), shat(n), t(n), v(n);

  // JL20141204: Added for initialization
  alpha = 0;
  omega = 1;
  rho1 = 1;
  rho2 = 1;

  r = b - A*x;
  rtilde = r;

  res = dotProd(r,r);
  if (printit)  std::cout << "itr " << 0 << ", (r,r)=" << res << "\n";

  tol *= res;
  tol = (atol > tol) ? atol : tol;

  if (res<=tol) {
    maxitr = 0;
    tol = res;
    return 0;
  }

  for (int i=1; i<=maxitr; i++) {
    rho1 = dotProd(rtilde,r);

    if (rho1==0) {
      tol = res;
      if (printit)  std::cout << "itr " << i << ", (r,r)=" << res << "\n";
      return 2;
    }

    if (i==1)
      p = r;
    else {
      beta = (rho1/rho2)*(alpha/omega);
      p = p - omega*v;
      p = r + beta*p;
    }

    phat = slvDiagSys(B,p);
    v = A*phat;

    alpha = rho1/dotProd(rtilde,v);
    s = r - alpha*v;
    res = dotProd(v,v);

    if (res<tol) {
      x = x + alpha*phat;
      tol = res;
      if (printit)  std::cout << "itr " << i << ", (s,s)=" << res << "\n";
      return 0;
    }

    if (printit)  std::cout << "itr " << i << ", (s,s)=" << res << ";";

    shat = slvDiagSys(B,s);
    t = A*shat;
    omega = dotProd(s,t)/dotProd(t,t);

    x += alpha*phat;
    x += omega*shat;
    r = s - omega*t;

    rho2 = rho1;
    res = dotProd(r,r);

    if (printit)  std::cout << "(r,r)=" << res << "\n";

    if (res<tol) {
      maxitr = i;
      tol = res;
      return 0;
    }

    if (omega==0) {
      tol = res;
      return 3;
    }
  }

  tol = res;
  return 1;
}


// Solving a sparse SPD system by the Conjugate Gradient method
// x brings in an initial guess and returns the solution 

void slvSpaSpdSysCG(Vector &x, 
  const SparseMatrix &A, const Vector &b,
  int &itr, int maxitr, double threshold, double tol) 
{
  int n;
  double alpha, b2, mu, nu, r2;

  b2 = b.l2norm();

  n = x.size();
  Vector r(n);
  Vector u(n);
  Vector v(n);

  itr = 0;
  r = b-A*x;
  v = r;
  mu = dotProd(r,r);

  std::cout << "itr/100= ";
  while (itr<maxitr) {
    if (itr%100==0)  std::cout << itr/100 << "  " << std::flush;
    if (v.l2norm()<threshold)  break;
    u = A*v;
    alpha = mu/dotProd(u,v);
    x += alpha*v;
    r -= alpha*u;
    nu = dotProd(r,r);
    r2 = sqrt(nu);
    if (r2<tol*b2)  break;
    v = r + (nu/mu)*v;
    mu = nu;
    itr++;
  }
  std::cout << "done!\n";

  return;
}


// Solving a sparse block SPD system by the Conjugate Gradient method
// x brings in an initial guess and returns the solution 

void slvSpaBlkSpdSysCG(Vector &x, const SparseBlockMatrix &A, const Vector &b,
  int &itr, int maxitr, double threshold, double tol) 
{
  int n;
  double alpha, b2, mu, nu, r2;
  
  std::cout << "In slvSpaBlkSpdSysCG...\n" << std::flush;

  b2 = b.l2norm();

  n = x.size();
  Vector r(n);
  Vector u(n);
  Vector v(n);

  itr = 0;
  r = b-A*x;
  v = r;
  mu = dotProd(r,r);

  std::cout << "itr= ";
  while (itr<maxitr) {
    if (itr%100==0)  std::cout << itr << "  " << std::flush;
    if (v.l2norm()<threshold)  break;
    u = A*v;
    alpha = mu/dotProd(u,v);
    x += alpha*v;
    r -= alpha*u;
    nu = dotProd(r,r);
    r2 = sqrt(nu);
    if (r2<tol*b2)  break;
    v = r + (nu/mu)*v;
    mu = nu;
    itr++;
  }
  std::cout << "done!\n";

  return;
}


// Solving a block diagonal SPD system 
// by the conjugate gradient method
// assuming each diagonal block is small and SPD

void slvBlkDiagSpdSysCG(Vector &x, const BlockDiagMatrix &A, 
  const Vector &b) 
{
  int i, itr, j, maxitr;
  double threshold=1E-9, tol=1E-9;
  FullMatrix B;
  Vector c, y; 

  for (i=1; i<=A.numBlocks(); ++i) {
    B = A.block(i);
    c.resize(A.dimBlock(i));
    y.resize(A.dimBlock(i));

    for (j=1; j<=A.dimBlock(i); ++j) 
      c(j) = b(A.bgnPosBlock(i)+j);

    maxitr = c.size();
    slvFullSpdSysCG(y, B, c, itr, maxitr, threshold, tol);

    for (j=1; j<=A.dimBlock(i); ++j) 
      x(A.bgnPosBlock(i)+j) = y(j);
  }

  return ;
}


// Solving a block diagonal SPD system by Cholesky facorization 
// assuming each diagonal block is small and SPD

void slvBlkDiagSpdSysCholesky(Vector &x, 
  const BlockDiagMatrix &A, const Vector &b) 
{
  int i, j;
  FullMatrix B;
  Vector c, y; 

  for (i=1; i<=A.numBlocks(); ++i) {
    B = A.block(i);
    c.resize(A.dimBlock(i));
    y.resize(A.dimBlock(i));

    for (j=1; j<=A.dimBlock(i); ++j) 
      c(j) = b(A.bgnPosBlock(i)+j);

    y = slvFullSpdSysCholesky(B, c);

    for (j=1; j<=A.dimBlock(i); ++j) 
      x(A.bgnPosBlock(i)+j) = y(j);
  }

  return ;
}


// Solving 

void slvBlkDiagSpdSysCholeskyLower(Vector &x, 
  const BlockDiagMatrix &A, const Vector &b) 
{
  int i, j;
  FullMatrix L;
  Vector c, y; 

  for (i=1; i<=A.numBlocks(); ++i) {
    L = A.block(i);
    c.resize(A.dimBlock(i));
    y.resize(A.dimBlock(i));

    for (j=1; j<=A.dimBlock(i); ++j) 
      c(j) = b(A.bgnPosBlock(i)+j);

    y = slvFullLowerTrigSys(L, c);

    for (j=1; j<=A.dimBlock(i); ++j) 
      x(A.bgnPosBlock(i)+j) = y(j);
  }

  return;
}


// Solving 

void slvBlkDiagSpdSysCholeskyUpper(Vector &x, 
  const BlockDiagMatrix &A, const Vector &b) 
{
  int i, j;
  FullMatrix U;
  Vector c, y; 

  for (i=1; i<=A.numBlocks(); ++i) {
    U = transpose(A.block(i));
    c.resize(A.dimBlock(i));
    y.resize(A.dimBlock(i));

    for (j=1; j<=A.dimBlock(i); ++j) 
      c(j) = b(A.bgnPosBlock(i)+j);

    y = slvFullUpperTrigSys(U, c);

    for (j=1; j<=A.dimBlock(i); ++j) 
      x(A.bgnPosBlock(i)+j) = y(j);
  }

  return;
}


// Solving a sparse block SPD system 
// by the Preconditioned Conjugate Gradient method 
// using the block diagonal SPD matrix as a preconditioner 

void slvSpaBlkSpdSysPCG(Vector &x, const SparseBlockMatrix &A, 
  const Vector &b, const BlockDiagMatrix &B, 
  int &itr, int maxitr, double threshold, double tol) 
{
  int n;
  double alpha, b2, mu, nu, r2;
 
  b2 = b.l2norm();
 
  n = x.size();
  Vector r(n), s(n), t(n), u(n), v(n);

  BlockDiagMatrix  C = CholeskyBlockDiagSPD(B);

  itr = 0;
  r = b-A*x;
  // slvBlkDiagSpdSysCG(s, B, r);
  // slvBlkDiagSpdSysCholesky(s, B, r);
  slvBlkDiagSpdSysCholeskyLower(t, C, r);
  slvBlkDiagSpdSysCholeskyUpper(s, C, t);
  mu = dotProd(r,s);
  v = s;

  // cout << "itr= ";
  while (itr<maxitr) {
    // if (itr%50==0)  std::cout << itr << "  ";
    if (v.l2norm()<threshold)  break;
    u = A*v;
    alpha = mu/dotProd(u,v);
    x += alpha*v;
    r -= alpha*u;
    // slvBlkDiagSpdSysCG(s, B, r);
    // slvBlkDiagSpdSysCholesky(s, B, r);
    slvBlkDiagSpdSysCholeskyLower(t, C, r);
    slvBlkDiagSpdSysCholeskyUpper(s, C, t);
    nu = dotProd(r,s);
    r2 = r.l2norm();
    if (r2<tol*b2)  break;
    v = s + (nu/mu)*v;
    mu = nu;
    itr++;
  }
  std::cout << "done!\n";

  return;
}


// Solving a block diagonal system by GEPP 

void slvBlkDiagSysGEPP(Vector &x, 
  const BlockDiagMatrix &A, const Vector &b) 
{
  int i, j;
  FullMatrix B;
  Vector c, y; 

  for (i=1; i<=A.numBlocks(); ++i) {
    B = A.block(i);
    c.resize(A.dimBlock(i));
    y.resize(A.dimBlock(i));

    for (j=1; j<=A.dimBlock(i); ++j) 
      c(j) = b(A.bgnPosBlock(i)+j);

    y = slvFullLinSysGEPP(B,c);

    for (j=1; j<=A.dimBlock(i); ++j) 
      x(A.bgnPosBlock(i)+j) = y(j);
  }

  return ;
}


// Solving a nonsymmetric block sparse linear system by 
// the BiCGStab (Preconditioned Bi-Conjugate Gradient Stablized)

int slvSpaBlkLinSysBiCGStab(Vector &x, const SparseBlockMatrix &A, 
  const Vector &b, const BlockDiagMatrix &B, 
  int &maxitr, double &tol, double atol, int printit)
{
  double alpha, beta, omega, res, rho1, rho2;

  int n = A.colSize();
  Vector p(n), phat(n), r(n), rtilde(n), s(n), shat(n), t(n), v(n);

  r = b - A*x;
  rtilde = r;

  res = dotProd(r,r);
  if (printit)  std::cout << "itr " << 0 << ", (r,r)=" << res << "\n";

  tol *= res;
  tol = (atol > tol) ? atol : tol;

  if (res<=tol) {
    maxitr = 0;
    tol = res;
    return 0;
  }

  for (int i=1; i<=maxitr; i++) {
    if (i%50==0)  cout << i << "  ";
    rho1 = dotProd(rtilde,r);

    if (rho1==0) {
      tol = res;
      if (printit)  std::cout << "itr " << i << ", (r,r)=" << res << "\n";
      return 2;
    }

    if (i==1)
      p = r;
    else {
      beta = (rho1/rho2)*(alpha/omega);
      p = p - omega*v;
      p = r + beta*p;
    }

    slvBlkDiagSysGEPP(phat,B,p);
    v = A*phat;

    alpha = rho1/dotProd(rtilde,v);
    s = r - alpha*v;
    res = dotProd(v,v);

    if (res<tol) {
      x = x + alpha*phat;
      tol = res;
      if (printit)  std::cout << "itr " << i << ", (s,s)=" << res << "\n";
      return 0;
    }

    if (printit)  std::cout << "itr " << i << ", (s,s)=" << res << ";";

    slvBlkDiagSysGEPP(shat,B,s);
    t = A*shat;
    omega = dotProd(s,t)/dotProd(t,t);

    x += alpha*phat;
    x += omega*shat;
    r = s - omega*t;

    rho2 = rho1;
    res = dotProd(r,r);

    if (printit)  std::cout << "(r,r)=" << res << "\n";

    if (res<tol) {
      maxitr = i;
      tol = res;
      return 0;
    }

    if (omega==0) {
      tol = res;
      return 3;
    }
  }

  tol = res;
  return 1;
}


// Subroutines for GMRES

// Applying a Givens rotation
inline void ApplGivensRot(double &dx, double &dy, double cs, double sn)
{
  double tmp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = tmp;
}

// Generating a Givens rotation 
inline void GenGivensRot(double dx, double dy, double &cs, double &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    double tmp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + tmp*tmp );
    cs = tmp * sn;
  } else {
    double tmp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + tmp*tmp );
    sn = tmp * cs;
  }
}

// Updating in GMRES
inline void UpdateGMRES(Vector &x, const Vector &a, const Vector *q, 
  const FullMatrix &H, int k)
{
  int i, j;
  Vector y(k);
  // Back substitution for solving upper trig. lin. sys.
  y(k) = a(k)/H(k,k);
  for (i=k-1; i>=1; --i) {
    y(i) = a(i);
    for (j=k; j>=i+1; --j)  y(i) = y(i) - H(i,j)*y(j);
    y(i) = y(i)/H(i,i);
  }
  // Initial guess plus linear combination of orthogonal vectors 
  for (j=1; j<=k; ++j)  x = x + y(j)*(q[j]);
}

// PreCondGMRES: Preconditioned Generalized Minimal RESidual 
// for solving a nonsymmetric sparse linear system
// cf: p.20 of the SIAM Templates book
// The return value indicates 
//   convergence within maxitr iterations (0) 
//   no convergence within maxitr iterations (1) 
// Upon successful return: 
//       x -- approximate solution to Ax = b
//  maxitr -- the number of iterations performed 
//            before the tolerance was reached
//     tol -- the residual after the final iteration

int slvSpaBlkSysPreCondGMRES(Vector &x, const Vector &b, 
  const SparseBlockMatrix &A, const BlockDiagMatrix &B, 
  int m, int &maxitr, double &atol, double &tol, int printit)
{
  int i, j, k, n;
  double b2, beta, res;
  n = A.colSize();
  Vector a(m+1), cs(m), sn(m);  // cosine,sine for Givens rotations
  Vector p(n), r(n), s(n), u(n);
  Vector *q = new Vector[m+2];  // orth. vectors in Krylov subspace
  for (k=1; k<=m+1; ++k)  q[k].resize(n);  // q[0] wasted
  FullMatrix H(m+1,m);

  slvBlkDiagSysGEPP(s,B,b);  // solve Bs=b
  b2 = s.l2norm();
  if (b2==0.0)  b2 = 1.0;

  r = b - A*x;
  slvBlkDiagSysGEPP(s,B,r);  // solve Bs=r
  beta = s.l2norm();  // beta=||s||

  res = beta/b2;
  if (res<=tol) {
    tol = res;
    maxitr = 0;
    return 0;
  }

  // if (printit)  std::cout << "residual= " << beta << "\n";

  tol *= b2;
  tol = (atol>tol) ? atol : tol;

  k = 1;
  while (k<=maxitr) {
    a(1) = beta;
    for (i=2; i<=m+1; ++i)  a(i) = 0.0;
    q[1] = (1.0/beta)*s;  // q[1]=s/||s||

    for (j=1; j<=m && k<=maxitr; ++j) {
      p = A*(q[j]);
      slvBlkDiagSysGEPP(u,B,p);  // solve Bu=A(q[j])

      // Arnoldi iteration
      for (i=1; i<=j; ++i) {
        H(i,j) = dotProd(q[i],u);
        u = u - H(i,j)*(q[i]);
      }
      H(j+1,j) = u.l2norm();  // H(j+1,j)=||u||
      q[j+1] = (1.0/H(j+1,j))*u;

      // Givens rotation for QR-factorization 
      for (i=1; i<=(j-1); ++i) {
        ApplGivensRot(H(i,j), H(i+1,j), cs(i), sn(i));
      }
      GenGivensRot(H(j,j), H(j+1,j), cs(j), sn(j));
      ApplGivensRot(H(j,j), H(j+1,j), cs(j), sn(j));
      ApplGivensRot(a(j), a(j+1), cs(j), sn(j));

      res = fabs(a(j+1));
      if (printit) {
        std::cout << "j=" << j << "  ";
        std::cout << "k=" << k << "  ";
        std::cout << "Residual=" << res << "  ";
        std::cout << "H(j,j)= " << H(j,j) << "  ";
        std::cout << "H(j+1,j)= " << H(j+1,j) << "  ";
        std::cout << "\n";
      }
      if (res<tol) {
        UpdateGMRES(x, a, q, H, j);
        tol = res;
        maxitr = k;
	    delete[] q;
        return 0;
      }
      k++;
    }

    UpdateGMRES(x, a, q, H, m);
    r = b - A*x;
    slvBlkDiagSysGEPP(s,B,r);  // solve Bs=r
    beta = s.l2norm();  // beta=||s||
    if (res<tol) {
      tol = res;
      maxitr = k;
      delete[] q;
      return 0;
    }

    if (printit)  std::cout << "Restarting..." << "\n";
  }

  tol = res;
  delete[] q;
  return(1);
}


// Solving a linear system based on Block-Diagonal-Schur-complement & CG
// Assuming compatible patterns, mainly for weak Galerkin FEMs

int slvBDSchurCG(Vector &x0, Vector &x1,
                 const BlockDiagMatrix &A00, const SparseBlockMatrix &A01,
                 const SparseBlockMatrix &A10, const SparseBlockMatrix &A11,
                 const Vector &b0, const Vector &b1,
                 int &itr, int maxitr, double threshold, double tol)
{
  int n;
  double alpha, beta2, mu, nu, r2;

  std::cout << "In slvBDSchurCG ...\n" << std::flush;

  n = b1.size();
  Vector b(n), r(n), u(n), v(n), x(n);
  
  BlockDiagMatrix A00inv;
  A00inv = inverseGEPP(A00);

  b = b1 - A10 * (A00inv * b0);
  beta2 = b.l2norm();

  // Initializing iteration
  // x = zero vector;
  itr = 0;
  r = b - (A11*x - A10 * (A00inv * (A01*x)));
  v = r;
  mu = dotProd(r,r);

  std::cout << "itr/100= ";
  while (itr<=maxitr) {
    if (itr%100==0)  std::cout << itr/100 << "  " << std::flush;
    if (v.l2norm()<threshold)  break;
    u = A11 * v - A10 * (A00inv * (A01*v));
    alpha = mu/dotProd(u,v);
    x += alpha * v;
    r -= alpha * u;
    nu = dotProd(r,r);
    r2 = sqrt(nu);
    if (r2 < tol * beta2)  break;
    v = r + (nu/mu) * v;
    mu = nu;
    itr++;
  }
  x1 = x;
  std::cout << "\n";

  // "solving" for x0
  x0 = A00inv * (b0 - A01 * x1);

  return(0);  // If successful
}

// LinSys.cpp
