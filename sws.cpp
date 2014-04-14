#include <iostream>
#include <fstream>
#include <cmath>
#include "matrixtypes.h"
#include "lapack.h"
#include "sws.h"
using namespace std;

double g_scalar_product(ula::RealMatrix& vr, int l, int r) {
 double res=0;
    for (unsigned int i=0;i<vr.size1()/2;++i)
      res+=vr(i,l)*vr(i,r);
    for (unsigned int i=vr.size1()/2;i<vr.size1();++i)
      res-=vr(i,l)*vr(i,r);
 return res;
}

complex<double> g_scalar_product(ula::ComplexMatrix& v, int l, int r) {
 complex<double> res=0;
    for (unsigned int i=0;i<v.size1()/2;++i)
      res+=conj(v(i,l))*v(i,r);
    for (unsigned int i=v.size1()/2;i<v.size1();++i)
      res-=conj(v(i,l))*v(i,r);
 return res;
}

void g_ortho(ula::RealMatrix& vr, int l, int r,int sgn) {
  for(int i=l;i<=r;++i) {
    for (int j=l;j<i;++j) {
      double coeff=-sgn*g_scalar_product(vr,j,i);
      for (unsigned int k=0;k<vr.size1();++k)
        vr(k,i)+=coeff*vr(k,j);
    }
    double inorm=1./sqrt(sgn*g_scalar_product(vr,i,i));
    for (unsigned int k=0;k<vr.size1();++k)
      vr(k,i)*=inorm;
  }
}

void g_ortho(ula::ComplexMatrix& v, int l, int r,double sgn) {
  for(int i=l;i<=r;++i) {
    for (int j=l;j<i;++j) {
      complex<double> coeff=-sgn*g_scalar_product(v,j,i);
      for (unsigned int k=0;k<v.size1();++k)
        v(k,i)+=coeff*v(k,j);
    }
    double inorm=1./sqrt(sqrt(norm(g_scalar_product(v,i,i))));
    for (unsigned int k=0;k<v.size1();++k)
      v(k,i)*=inorm;
  }
}

int main()
{
  nu =0.0;
//  for(nu = 0.0; nu < 4.01; nu+=.5){
  uni_mag = 0.0;
  sta_mag = 0.0;
  ene = 0.0;
    
  IntConf();
  BuildMatrix();
  DialMatrix();
  CalcEne();
  CalcMag();
  PrintOut();
    
}

/*------------------------------------------------------------------------------
Read the coordinates of spins and count number of spins for each kind of lattice 
------------------------------------------------------------------------------*/
void IntConf()
{
  ifstream InPos("site.dat");
  int countA, countB, x, y, index;

  InPos >> NoSpin;                                    /* Read number of spins */
  for(int i = 0; i < NoSpin; ++i)                   /* Read position of spins */
    InPos >> PosSpin[i];
  
  countA = 0;
  countB = N/2; 
  for (int i = 0; i < Nx; ++ i)     /* Assigne each position of site with an index */
    for (int j = 0; j < Ny; ++ j)
      if (sign(i,j)) {
	ind[i][j] = countA;
	countA++;
      }
      else {
	ind[i][j] = countB;
	countB++;
      }

  NB = 0;                        /* Count number of spins in lattice A and B */
  for(int i = 0; i < NoSpin; ++i){
    x = PosSpin[i]/Ny;
    y = PosSpin[i]%Ny; 
    index = ind[x][y];
    eta(index) = 1.;
    if(index >= N/2 && eta(index) != 0)
      ++NB;
  }
}

void BuildMatrix()
{
  int i1, i2, j1, j2, tmp;
  Nb=0;
  for (i = 0; i < Nx; ++ i)
    for (j = 0; j < Ny; ++ j) {
      if (i < Nx-1) {
	i1 = ind[i][j];
	i2 = ind[i+1][j];
	if (eta(i1)*eta(i2) !=0 ) {
	  if (i1>i2) {
	    tmp = i1;
	    i1 = i2;
	    i2 = tmp;
	  }
	  m(i1,i1) += S;
	  m(i2,i2) += S;
	  m(i1,i2) += S*d;
	  m(i2,i1) += S*d;
	  ++Nb;
	}
      }
      if (j < Ny-1) {
	j1 = ind[i][j];
	j2 = ind[i][j+1];
	if (eta(j1)*eta(j2)!=0) { 
	  if (j1>j2) {
	    tmp = j1;
	    j1 = j2;
	    j2 = tmp;
	  }
	  m(j1,j1) += S;
	  m(j2,j2) += S;
	  m(j1,j2) += S*d;
	  m(j2,j1) += S*d;
	  ++Nb;
	}
      }
    }
  for (int i=0;i<N;++i)  {
    g(i,i)=(i < N/2) ? 1 : -1;
  }
  m = prod(g, m - nu*S*g);
}


void DialMatrix()
{
  ula::ns_diagx(m,e,v);
  er = real(e);
  ula::eigensort(er,v);
  vr = real(v);
}

void CalcEne()
{
  int left=0;
  double energy = er(left);
  for (i = 1; i<= N; ++i) {
    if (i==N || fabs(energy-er(i))>1E-10) {
      if (fabs(energy)>0) {
	g_ortho(v,left,i-1,(energy>0) ? 1 : -1 );
      }
      left = i;
      energy = er(left);
    }
  }
}

/*---------------------------------------------------------------------------
Calculate the manetiztion of each site as well as the uniform and staggered 
magnetization of lattice
----------------------------------------------------------------------------*/   
void CalcMag()
{
  double mA, mB, ms;
  int i, j, k, index;

  for(i = 0; i < Nx; ++ i)
  for(j = 0; j < Ny; ++ j)
    msH(i,j)=0.0;

  mA=0.0;
  mB=0.0;
  for (i = 0; i < Nx; ++ i)
    for (j = 0; j < Ny; ++ j){
      index = ind[i][j];
      ms = 0;
      if ((int)sign(i,j)) {
	for (k = 0; k < N; ++ k){
	  if (er(k)<-e_zero && fabs(er(k)+nu*S)>eps){   //eleminate the goldstone modes
	    ms += vr(index,k)*vr(index,k);
	  }
	}	
	msH(i,j) = (S-ms)*eta(index);
	mA += msH(i,j);
      }
      else {
	for (k = 0; k < N; ++k ){
	  if (er(k)>e_zero && fabs(er(k)+nu*S)>eps){
	    ms += vr(index,k)*vr(index,k);
	  }
	}
	msH(i,j) = (-S+ms)*eta(index);
	mB += msH(i,j);
      }
    }
  uni_mag = -(mA+mB)/(1.0*NoSpin);
  sta_mag = (mA-mB)/(1.0*NoSpin);
}

void PrintOut()
{
  ofstream out_sws("spinsite.dat");
  ofstream out_sta("sta_mag.dat");
  ofstream out_uni("uni_mag.dat");
  ofstream out_ene("energy.dat");

  if(nu==0) {
    for(i = 0; i < Nx; i++)
    for(j = 0; j < Ny; j++)
      if(msH(i,j)!=0){
	out_sws << i << "   ";
      }
    out_sws << "\n";
    
    for(i = 0; i < Nx; i++)
    for(j = 0; j < Ny; j++)
      if(msH(i,j)!=0){
	out_sws << j << "   ";
      }
    out_sws << "\n";
  }
  
  for(i = 0; i < Nx; i++)
  for(j = 0; j < Ny; j++)
    if(msH(i,j)!=0){
      out_sws << msH(i,j) << "   ";
    }
  out_sws<<endl;
  
  for(i = 0; i < N; ++i)
    out_ene << er(i) << endl;
  out_ene<<endl;

  out_uni << nu <<"    "<< uni_mag << "   " << endl;
  out_sta << nu <<"    "<< sta_mag << "   " << endl;
}

