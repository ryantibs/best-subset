#include <R.h>
#include <Rmath.h>

// Matrices are stored as vectors, in column-major order

// Givens rotation of a and b, stored in c and s
void givens(double a, double b, double *c, double *s) {
  if (b==0) {
    *c = 1; 
    *s = 0; 
  }
  else {
    if (fabs(b)>fabs(a)) {
      double t = -a/b;
      *s = 1/sqrt(1+t*t);
      *c = (*s)*t;
    }
    else {
      double t = -b/a;
      *c = 1/sqrt(1+t*t);
      *s = (*c)*t;
    }
  }
}

// Givens rotation applied to rows i1 and i2 of the m x n
// matrix A, on the subset of columns j1 through j2
void rowrot(double *A, int i1, int i2, int m, int n, int j1, int j2, double c, double s) {
  int j;
  double t1,t2;
  for (j=j1; j<=j2; j++) {
    t1 = A[i1+j*m]; 
    t2 = A[i2+j*m];
    A[i1+j*m] = c*t1-s*t2;
    A[i2+j*m] = s*t1+c*t2;
  }
}

// Givens rotation applied to columns j1 and j2 of the m x n
// matrix A, on the subset of rows i1 through i2
void colrot(double *A, int j1, int j2, int m, int n, int i1, int i2, double c, double s) {
  int i;
  double t1,t2;
  for (i=i1; i<=i2; i++) {
    t1 = A[i+j1*m];
    t2 = A[i+j2*m];
    A[i+j1*m] = c*t1-s*t2;
    A[i+j2*m] = s*t1+c*t2;
  }
}

// Downdate the QR factorization after deleting column j0, 
// where Q1 is m x n and R is n x n. The other part of 
// the Q matrix, Q2 m x (m-n), isn't needed so it isn't 
// passed for efficiency
void downdate1(double *Q1, double *R, int *j0p, int *mp, int *np) {
  int j0,m,n,j;
  j0 = *j0p;
  m = *mp;
  n = *np;

  double c,s;
  for (j=j0+1; j<n; j++) {
    // Compute the appropriate c and s
    givens(R[j-1+j*n],R[j+j*n],&c,&s);

    // Pre-multiply R
    rowrot(R,j-1,j,n,n,j,n-1,c,s);

    // Post-mutiply Q1
    colrot(Q1,j-1,j,m,n,0,m-1,c,s);
  }
}

// Update the QR factorization after adding a column z.
// For convenience, we are given w=Q2'z. Here Q2 is m x 
// (m-n). The other part of the Q matrix, Q1 m x n, and 
// R are not needed, so they aren't passed for efficiency 
void update1(double *Q2, double *w, int *mp, int *kp) {
  int m,k,j;
  m = *mp;
  k = *kp;

  double c,s;
  for (j=k-1; j>=1; j--) {
    // Compute the appropriate c and s
    givens(w[j-1],w[j],&c,&s);

    // Pre-multiply w
    rowrot(w,j-1,j,k,1,0,0,c,s);

    // Post-multiply Q2
    colrot(Q2,j-1,j,m,k,0,m-1,c,s);
  }
}

// Downdate the QR factorization after deleting the first row, 
// where Q is m x m and R is m x n
void downdate2(double *Q, double *R, int *mp, int *np) {
  int m,n,i;
  m = *mp;
  n = *np;

  double c,s;
  for (i=m-1; i>=1; i--) {
    // Compute the appropriate c and s
    givens(Q[(i-1)*m],Q[i*m],&c,&s);

    // Post-mutiply Q
    colrot(Q,i-1,i,m,m,0,m-1,c,s);

    // Pre-multiply R
    if (i<=n) rowrot(R,i-1,i,m,n,i-1,n-1,c,s); 
  }
}

// Update the QR factorization after adding the last row,
// where Q is m x m and R is m x n. For efficiency, Q is not
// passed, and only the first row of R is passed. Not counting
// its first row, the first q columns of R are zero
void update2(double *y, double *D, double *r, int *mp, int *np, int *qp) {
  int m,n,q,j;
  m = *mp;
  n = *np;
  q = *qp;

  double c,s;
  for (j=0; j<q-1; j++) {
    // Compute the appropriate c and s
    givens(r[j+1],r[j],&c,&s);

    // Pre-multiply r
    rowrot(r,j+1,j,n,1,0,0,c,s);

    // Post-mutiply D
    colrot(D,j+1,j,m,n,0,m-1,c,s);

    // Pre-multiply y
    rowrot(y,j+1,j,n,1,0,0,c,s);
  }
}

// Make the R factor upper triangular, by Givens rotating its
// columns. Here A is m x n, and R is n x n with rank(R) = n-k
void maketri1(double *y, double *A, double *R, int *mp, int *np, int *kp) {
  int m,n,k,i,j;
  m = *mp;
  n = *np;
  k = *kp; 
  
  double c,s;
  for (i=n-k-1; i>=0; i--) {
    for (j=i; j<i+k; j++) {
      // Compute the appropriate c and s
      givens(R[i+(j+1)*n],R[i+j*n],&c,&s);
      
      // Post-multiply R
      colrot(R,j+1,j,n,n,0,i,c,s);

      // Post-multiply D
      colrot(A,j+1,j,m,n,0,m-1,c,s);
    
      // Pre-multiply y
      rowrot(y,j+1,j,n,1,0,0,c,s);
    }
  }
}

// Make the R factor upper triangular, by Givens rotating its 
// columns. Here A is m x n, and R is m x n with rank(R) = min(m,n)-k
void maketri2(double *y, double *A, double *R, int *mp, int *np, int *kp) {
  int m,n,k,r,d,i,j;
  m = *mp;
  n = *np;
  k = *kp;
 
  r = imin2(m,n);
  d = imax2(n-m,0);

  double c,s;
  for (i=r-k-1; i>=0; i--) {
    for (j=i; j<i+k+d; j++) {
      // Compute the appropriate c and s
      givens(R[i+(j+1)*m],R[i+j*m],&c,&s);
      
      // Post-multiply R
      colrot(R,j+1,j,m,n,0,i,c,s);

      // Post-multiply D
      colrot(A,j+1,j,m,n,0,m-1,c,s);
    
      // Pre-multiply y
      rowrot(y,j+1,j,n,1,0,0,c,s);
    }
  }
}

// Make the R factor upper triangular, by Givens rotating
// its columns. Here A is m1 x n, and R is m2 x n with 
// rank(R) = n-q-1. The first q columns of R are zero
void maketri3(double *y, double *A, double *R, int *m1p, int *m2p, int *np, int *qp) {
  int m1,m2,n,q,i;
  m1 = *m1p;
  m2 = *m2p;
  n = *np;
  q = *qp;
  
  double c,s;
  for (i=n-q-2; i>=0; i--) {
    // Compute the appropriate c and s
    givens(R[i+(i+q+1)*m2],R[i+(i+q)*m2],&c,&s);
      
    // Post-multiply R
    colrot(R,i+q+1,i+q,m2,n,0,i,c,s);

    // Post-multiply D
    colrot(A,i+q+1,i+q,m1,n,0,m1-1,c,s);
    
    // Pre-multiply y
    rowrot(y,i+q+1,i+q,n,1,0,0,c,s);
  }
}

// Make the R factor upper triangular, by Givens rotating
// its columns and rows, appropriately. Here A is m1 x n, 
// Q is m2 x m2, and R is m2 x n with rank(R) = n-q-1. The 
// first q columns of R are zero. The kth row of R is the 
// last row with a zero element on the diagonal    
void maketri4(double *y, double *A, double *Q, double *R, int *m1p, int *m2p, int *np, int *qp, int *kp) {
  int m1,m2,n,q,k,i,j;
  m1 = *m1p;
  m2 = *m2p;
  n = *np; 
  q = *qp;
  k = *kp;

  double c,s;

  // First rotate the columns
  for (i=k-1; i>=0; i--) {
    // Compute the appropriate c and s
    givens(R[i+(i+q+1)*m2],R[i+(i+q)*m2],&c,&s);
      
    // Post-multiply R
    colrot(R,i+q+1,i+q,m2,n,0,i,c,s);

    // Post-multiply D
    colrot(A,i+q+1,i+q,m1,n,0,m1-1,c,s);
    
    // Pre-multiply y
    rowrot(y,i+q+1,i+q,n,1,0,0,c,s);
  }

  // Next rotate the rows
  for (j=k+q+1; j<n; j++) {
    // Compute the appropriate c and s
    givens(R[j-q-1+j*m2],R[j-q+j*m2],&c,&s);

    // Pre-multiply R
    rowrot(R,j-q-1,j-q,m2,n,j,n-1,c,s);

    // Post-multiply Q
    colrot(Q,j-q-1,j-q,m2,m2,0,m2-1,c,s);
  }
}

