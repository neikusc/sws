#ifndef WESSEL_ULA_LAPACK
#define WESSEL_ULA_LAPACK

#include "matrixtypes.h"
#include <complex>

namespace ula {
  
  void _message( int info, int lwork ) {
    if( info > 0 ) 
      std::cerr << "Warning: diagonalization failed to converge." << std::endl;
    else 
      if( info < 0 ) 
	std::cerr << "Warning: illegal argument value in diagonalization." << std::endl;
#ifdef BL_DEBUG
      else std::cerr << "Optimal cache size is lwork = " << lwork << "." << std::endl;
#endif
    
  }
}

/*
  ---------------------------------------------------
  
  DIAGONALIZATION OF REAL NON-SYMMETRIC MATRIX
  
  H: n-by-n input matrix
  E: eigenvalues
  
  ---------------------------------------------------
*/

extern "C" void dgeev_(char* jobvl, char* jobvr, 
		       int* n, double* a, int* lda, 
		       double* er, double* ei, 
		       double* vl, int* ldvl, double* vr, int* ldvr,
                       double* work, int* lwork, int* info);

extern "C" void dgeevx_(char* balanc, 
                        char* jobvl, char* jobvr, 
                        char* sense,
		        int* n, double* a, int* lda, 
		        double* er, double* ei, 
		        double* vl, int* ldvl, double* vr, int* ldvr,
                        int* ilo, int* ihi, double* scale, double* abnrm, double* rconde, double* rcondv,
                        double* work, int* lwork, 
                        int* iwork,
                        int* info);

namespace ula {
  
  void ns_diag( RealMatrix& H, ComplexVector& E) {
    RealMatrix H_=H;
    int calc_evl=0;
    int calc_evr=0;
    char jobvl = ( calc_evl ? 'V' : 'N' );  // calculate left  eigenvectors?
    char jobvr = ( calc_evr ? 'V' : 'N' );  // calculate right eigenvectors?
    int n = H.size1();                      // squareness is not checked
    int lda = n;                            // leading dimension;
    double* er = new double[n];             // real part of eigenvalues
    double* ei = new double[n];             // imag part of eigenvalues
    int ldvl=1;
    double* vl = new double[ldvl];
    int ldvr=1;
    double* vr = new double[ldvr];
    int lwork=3*n;
    double* work = new double[lwork];
    int info;
    dgeev_( &jobvl, &jobvr, 
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            work, &lwork, &info );  
    _message( info, lwork );
    for (int i=0;i<n;++i)
      E(i)=std::complex<double>(er[i],ei[i]);
    delete[] er;
    delete[] ei;
    delete[] vl;
    delete[] vr;
    delete[] work;
  }

  void ns_diag( RealMatrix& H, ComplexVector& E, ComplexMatrix& V) {
    RealMatrix H_=H;
    int calc_evl=0;
    int calc_evr=1;
    char jobvl = ( calc_evl ? 'V' : 'N' );  // calculate left  eigenvectors?
    char jobvr = ( calc_evr ? 'V' : 'N' );  // calculate right eigenvectors?
    int n = H.size1();                      // squareness is not checked
    int lda = n;                            // leading dimension;
    double* er = new double[n];             // real part of eigenvalues
    double* ei = new double[n];             // imag part of eigenvalues
    int ldvl=1;
    double* vl = new double[ldvl*ldvl];
    int ldvr=n;
    double* vr = new double[ldvr*ldvr];
    int lwork=5*n;
    double* work = new double[lwork];
    int info;
    dgeev_( &jobvl, &jobvr, 
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            work, &lwork, &info );  
    _message( info, lwork );
    for (int i=0;i<n;++i) {
      if (ei[i]) {
        E(i)=std::complex<double>(er[i],ei[i]);
        E(i+1)=std::complex<double>(er[i+1],ei[i+1]);
        for (int j=0;j<n;++j) {
          V(j,i)=std::complex<double>(vr[n*i+j],vr[n*(i+1)+j]);
          V(j,i+1)=std::complex<double>(vr[n*i+j],-vr[n*(i+1)+j]);
        }
        ++i;
      }
      else {
        E(i)=std::complex<double>(er[i],ei[i]);
        for (int j=0;j<n;++j)
          V(j,i)=vr[n*i+j];
      }
    }   
    delete[] er;
    delete[] ei;
    delete[] vl;
    delete[] vr;
    delete[] work;
  }

  void ns_diagx( RealMatrix& H, ComplexVector& E, ComplexMatrix& V) {
    RealMatrix H_=H;
    int calc_evl=1;
    int calc_evr=1;
    char balanc= 'B';
    char jobvl = ( calc_evl ? 'V' : 'N' );  // calculate left  eigenvectors?
    char jobvr = ( calc_evr ? 'V' : 'N' );  // calculate right eigenvectors?
    char sense = 'B';
    int n = H.size1();                      // squareness is not checked
    int lda = n;                            // leading dimension;
    double* er = new double[n];             // real part of eigenvalues
    double* ei = new double[n];             // imag part of eigenvalues
    int ldvl=n;
    double* vl = new double[ldvl*ldvl];
    int ldvr=n;
    double* vr = new double[ldvr*ldvr];
    int ilo;
    int ihi;
    double* scale=new double[n];
    double abnrm;
    double* rconde=new double[n];
    double* rcondv=new double[n];
    int lwork=n*(n+6);
    double* work = new double[lwork];
    int* iwork = new int[2*n-2];
    int info;
    dgeevx_( &balanc,
            &jobvl, &jobvr, 
            &sense,
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            &ilo, &ihi , scale, &abnrm, rconde, rcondv,
            work, &lwork, 
            iwork,
            &info );  
    _message( info, lwork );
    for (int i=0;i<n;++i) {
      if (ei[i]) {
        E(i)=std::complex<double>(er[i],ei[i]);
        E(i+1)=std::complex<double>(er[i+1],ei[i+1]);
        for (int j=0;j<n;++j) {
          V(j,i)=std::complex<double>(vr[n*i+j],vr[n*(i+1)+j]);
          V(j,i+1)=std::complex<double>(vr[n*i+j],-vr[n*(i+1)+j]);
        }
        ++i;
      }
      else {
        E(i)=std::complex<double>(er[i],ei[i]);
        for (int j=0;j<n;++j)
          V(j,i)=vr[n*i+j];
      }
    }   
    delete[] er;
    delete[] ei;
    delete[] vl;
    delete[] vr;
    delete[] work;
  }
}

/*
  ---------------------------------------------------
  
  DIAGONALIZATION OF REAL SYMMETRIC MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  E: eigenvalues
  V: eigenvectors
  
  ---------------------------------------------------
*/

extern "C" void dsyev_(char* jobz, char* uplo, 
                       int* n, double* a, int* lda, 
                       double* w, 
                       double* work, int* lwork,
                       int* info);

namespace ula {
  
  void _diag( RealMatrix& H, RealVector& E, bool calc_ev ) {
    char jobz = ( calc_ev ? 'V' : 'N' );  // calculate eigenvectors?
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    int lda = n;                          // leading dimension;
    int lwork = 3*n-1;
    double* work = new double[lwork];
    int info;
    dsyev_( &jobz, &uplo, 
	    &n, &(H(0,0)), &lda, 
	    &(E(0)), 
	    work, &lwork, 
	    &info);  
    _message( info, lwork );
    delete[] work;
  }
  
  
  inline void diag( RealMatrix& H, RealVector& E ) {
    _diag( H, E, false );
  }
  
  inline void diag( RealMatrix& H, RealVector& E, RealMatrix& V ) {
    V=H;
    _diag( V, E, true );
  }
}

/*
  ---------------------------------------------------
  
  DIAGONALIZATION OF A HERMITIAN MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  E: eigenvalues
  V: eigenvectors
  
  ---------------------------------------------------
*/

extern "C" void zheev_( char* jobz, char* uplo, 
                        int* n, std::complex<double>* a, int* lda, 
                        double* w,
                        std::complex<double>* work, int* lwork, double* rwork,
                        int* info);

namespace ula {
  
  void _cdiag(ComplexMatrix& H, RealVector& E, bool calc_ev ) {
    char jobz = ( calc_ev? 'V' : 'N' );   // calculate eigenvectors?
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    int lda = n;                          // leading dimension
    int lwork = 3*n-1;                    // check optimal size in lwork after call
    int lrw = 3*n-2;                      // this size is hardwired
    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[lrw];
    int info;
    zheev_( &jobz, &uplo, 
	    &n, &(H(0,0)), &lda, 
	    &(E(0)), 
	    work, &lwork, rwork, 
	    &info);
    _message( info, lwork );
    delete[] work;
    delete[] rwork;
  }
  
  void diag(ComplexMatrix& H, RealVector& E) {
    _cdiag( H, E, false );
  }
  
  void diag(ComplexMatrix& H, RealVector& E, ComplexMatrix& V) {
    V=H;
    _cdiag( V, E, true );
  }
}

/*
  ---------------------------------------------------
  
  DETERMINANT OF A HERMITIAN MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  
  ----------------------------------------------------
*/

extern "C" void zgetrf_(int *m, int *n, 
                        std::complex<double> *a, int *lda, 
                        int* ipiv, 
                        int* info );

namespace ula {
  
  std::complex<double> determinant( ComplexMatrix& H ) {
    int n = H.size1();      // squareness is not checked!!!
    int lda = n;
    int* ipiv = new int[n];
    int info;
    zgetrf_( &n, &n, 
	     &(H(0,0)), &lda, 
	     ipiv, 
	     &info );
    int sign = 1;
    for( int i = 1; i <= n; i++ )
      if( ipiv[i-1] != i ) sign = -sign;
    std::complex<double> det = 1.0*sign;
    for( int i = 0; i < n; i++ ) det *= H(i,i*n);
    delete[] ipiv;
    return det;
  }
}


namespace ula {
  
  /*
    ---------------------------------------------------
    
    SORT ROUTINE FOR REAL EIGENVALUES
    
    ----------------------------------------------------
  */
  
  template <typename valuetype>
    void eigensort(RealVector& d, boost::numeric::ublas::matrix<valuetype, boost::numeric::ublas::column_major, std::vector<valuetype > >& v) {
    int k,j,i;
    int n=v.size1();
    for (i=0;i<n-1;++i) {
      double p=d(k=i);
      for (j=i+1;j<n;++j)
	if (d(j)>=p)
	  p=d(k=j);
      if (k!=i) {
	d(k)=d(i);
	d(i)=p;
	for (j=0;j<n;++j) {
	  valuetype v_=v(j,i);
	  v(j,i)=v(j,k);
	  v(j,k)=v_;
	};
      }
    };
  }
}


/*
  ---------------------------------------------------

  INVERSE OF REAL MATRIX

  A : n-by-n input matrix
  A_: A^{-1}

  ----------------------------------------------------
*/

extern "C" void dgesv_(int *n, int *nrhs,
              double *a, int *lda,
              int *ipiv,
              double *b, int *ldb,
              int* info );

namespace ula {

  void invr(RealMatrix& A, RealMatrix& A_) {
    int n = A.size1();
    RealMatrix B(n,n);
    int lda = n;
    int ldb = n;
    int info;
    int nrhs=n;
    int* ipiv = new int[n];
    A_=A;
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
       if(i==j)B(i,i)=1;
       else B(i,j)=0;
    dgesv_( &n, &nrhs ,
           &(A_(0,0)), &lda,
           ipiv,
           &(B(0,0)), &ldb,
           &info );
    if(info==0)A_=B;
    else std::cout<<"Info="<<info<<std::endl;
  }
}



#endif

