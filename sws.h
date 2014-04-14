#define S 0.5                                                 /* Spin */
#define d 1.0                       /* The strength of XY interaction */
#define eps 1E-15                    /* Eleminate the goldstone modes */
#define e_zero 1E-10        /* Cut off energy levels which is below 0 */
#define Nx 8                         /* The dimensions of the lattice*/
#define Ny 8 
#define N Nx*Ny                                /* The number of sites */

ula::RealMatrix m(N,N);
ula::RealMatrix g(N,N);
ula::RealMatrix vr(N,N);
ula::RealVector er(N);
ula::RealVector eta(N); /* Each site is assignes to index 1 or 0 realatively to a spin or an impurity */
ula::RealMatrix msH(N,N);
ula::ComplexMatrix v(N,N);
ula::ComplexVector e(N);

int i, j, k;
int ind[N][N];                     /* The index of spin in the lattice */
int NoSpin;           /* The number of spins in the lattice */
int NB;                /* The real number of spins belong to lattice B */
int Nb;                /* The number of bonds */
int PosSpin[N];       /* The position of spin in one dimensional array */
double uni_mag, sta_mag, ene;
double nu;

void IntConf();                                /* Prototypes functions */
void BuildMatrix();
void DialMatrix();
void CalcEne();
void CalcMag();
void PrintOut();
int sign(int ix,int iy) { 
  return (ix+iy)%2; 
}

