#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);

void step(cmplx* const psi0,cmplx* const psi1,const double dx, const double dt, const double omega, const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40 ;
    const double xmax = 40;
	const double Tend = 10 * M_PI;
	const double dx = (2 * xmax)/(Nx) ;
	const double dt = dx / 1000.0 ;
    double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
    const double omega = 0.2;
    const double alpha = sqrt(omega);

    stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
    cmplx* psi1 = new cmplx[Nx];
    cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
            
            step(psi0, psi1, dx, dt, omega, Nx,xmin);
            h = psi0;
            psi0 = psi1;
            psi1 = h;
            
            
            t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;

    delete[] psi0;
    delete[] psi1;
    
  return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}

//-------------------------------------

void step(cmplx* const psi0, cmplx* const psi1, const double dx, const double dt, const double omega, const int Nx, const double xmin){
                        cmplx* d = new cmplx[Nx]; // initialise complex arrays for entries of matrix A
                        cmplx* a = new cmplx[Nx];
    
                        for (int i = 0; i<Nx; i++) { // filling matrix A with values for Schroedinger equation
                            d[i] = cmplx(1, dt / (2.0 * pow(dx,2)) + dt * pow(omega,2) * pow(xmin+i*dx,2) / 4.0);
                            a[i] = cmplx(0, -dt / 4 * pow(dx,2));
                        }
    
    
                        cmplx* psi2 = new cmplx[Nx];
    
    
                        psi2[0] = conj(d[0])*psi0[0] + conj(a[0])*psi0[1]; // calculate A^* psi_0
                        for( int i = 1; i< Nx-1; i++){
                            psi2[i] = conj(a[i])*psi0[i-1] + conj(d[i])*psi0[i] + conj(a[i])*psi0[i+1];
                        }
                        psi2[Nx-1] = conj(d[Nx-1])*psi0[Nx-1] + conj(a[Nx-1])*psi0[Nx-2];
    
    
                         for (int i = 1; i< Nx; i++ ){
                             d[i] = d[i] - a[i]/d[i-1] * a[i-1]; // upper triangular form of matrix A
                             psi2[i] = psi2[i] - a[i]/d[i-1] * psi2[i-1]; // apply same operations to right hand side of equation
                         }
    
    
    
                         // backward substitution for left hand side
                         psi1[Nx-1] = psi2[Nx-1] / d[Nx-1];
                         for(int i = Nx-2;i>=0;i--){
                             psi1[i] = (psi2[i] - a[i] * psi1[i+1]) / d[i];
                         }
    
                         delete[] d;
                         delete[] a;
                         delete[] psi2;
    
}
