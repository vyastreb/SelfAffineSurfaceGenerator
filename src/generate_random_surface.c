//===============================================//
//  Rough Surface Generator (RSG)		 //
//                                               //
//  V.A. Yastrebov, 2013-2015                    //
//  Centre des Materiaux, CNRS, MINES ParisTech  //
//===============================================//
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "normal.h"

using namespace::std;

double my_random() {
  return random()/double(RAND_MAX);
  //int r = random(); 
  //return i4_normal_ab(0,10000,r)/30000.; 
}
void logo() {
  cout << endl << "/=RSG Rough Surface Generator=========================================" << endl;
}
void message(string msg_) {
  cout << "/ " << msg_ << endl;
}

int main(int argc, char* argv[]) {
  logo();
// Read input
  if (argc != 9) { 
    cout << "I need 8 arguments: [1] lower cutoff wavenumber k1 (int)" << endl 
	 << "                    [2] upper cutoff wavenumber k2 (int)" << endl
	 << "                    [3] Hurst exponent, H (double) // 0 < H < 1" << endl
	 << "                    [4] number of points per size,  L (int) // powers of 2, i.e. 128, 256, 512, etc" << endl
	 << "                    [5] seed for the random number generator, s (int)" << endl
	 << "                    [6] standard deviation of heights, rms (double)" << endl
	 << "                    [7] number of bins in pdf data, Npdf (int)" << endl
	 << "                    [8] if the power spectral density has a plateau up to k1? (0 - non, 1 -yes)" << endl;
    abort();
  }
  int    k1        = atoi(argv[1]);
  int    k2        = atoi(argv[2]);
  double H         = atof(argv[3]);
  int    L         = atoi(argv[4]);
  int    seed      = atoi(argv[5]);
  double rms       = atof(argv[6]);
  int    Npdf      = atoi(argv[7]);
  bool has_plateau = atoi(argv[8]);
  double dx        = 1.;
  int    L2        = L*L;

// Predefine some variables
  ofstream out,out2;
  int k1_2=k1*k1;
  int k2_2=k2*k2;
  srand(seed);
  std::ostringstream complement;
  complement << "_k1_" << k1 << "_k2_" << k2 << "_H_" << H << "_seed_" << seed << "_L_" << L << "_rms_" << rms << ".data";

// Generate initial random distribution R and transform it to FFT
  message("Constructing initial white noise.");
  double zmin=1e15,zmax=-1e15;
  fftw_plan p;
  fftw_complex* r, *r_fft;
  r     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L2);
  r_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L2);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      r[i+j*L][0] = my_random(); 
      r[i+j*L][1] = 0.;
    }
  }
  message("Transforming the noise in Fourier space.");
  p = fftw_plan_dft_2d(L, L, r, r_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_free(r); 

// Filter H_FFT
  message("Constructing filter in Fourier space.");
  fftw_complex* h_fft;
  h_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L * L);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
	h_fft[i+L*j][0] = 0;
	h_fft[i+L*j][1] = 0;
    }
  }
      
  int nb_modes = 0; 
  if (k2>L/2) {
    cout << "Error: k2 must hot exceed L/2!" << endl;
    abort();
  }
  for (int i = 0; i < k2; i++) {
    int i_2=i*i;
    for (int j = 0; j < k2; j++) {
      int j_2=j*j;
      if (has_plateau) {
              double plateau = pow(k1,-1.*(H+1));
	      if ( i_2 + j_2 < k1_2 && i_2 + j_2 > 0 ) {
		if (i==0)  
		  h_fft[i+L*j][0] = h_fft[i+L*(L-j)][0] = plateau; 
		else if (j==0)  
		  h_fft[i+L*j][0] = h_fft[L-i+L*j][0] = plateau;
		else 
		  h_fft[i+L*j][0] = h_fft[i+L*(L-j)][0] = h_fft[L-i+L*j][0] = h_fft[L-i+L*(L-j)][0] = plateau;
	      }
      }
      if ( i_2 + j_2 >= k1_2 &&  i_2+j_2 <= k2_2 ) {
	nb_modes++;
	double A = pow(sqrt(i_2+j_2),-1.*(H+1));
	if (i==0 && j==0)
          h_fft[i+L*j][0] = 0; 
	else if (i==0)  
          h_fft[i+L*j][0] = h_fft[i+L*(L-j)][0] = A; 
	else if (j==0)  
          h_fft[i+L*j][0] = h_fft[L-i+L*j][0] = A; 
	else 
          h_fft[i+L*j][0] = h_fft[i+L*(L-j)][0] = h_fft[L-i+L*j][0] = h_fft[L-i+L*(L-j)][0] = A; 
      }
    }
  }

  message("Computing product in Fourier space.");
// Constructs the random surface <z_fft> in Fourier space
  fftw_complex* z_fft, *z;
  z_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L * L);
  for (int i = 0; i < L; i++) 
    for (int j = 0; j < L; j++)
      z_fft[i+j*L][0] = z_fft[i+j*L][1] = 0; 
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      z_fft[i+j*L][0] = h_fft[i+j*L][0] * r_fft[i+j*L][0] - h_fft[i+j*L][1] * r_fft[i+j*L][1];
      z_fft[i+j*L][1] = h_fft[i+j*L][0] * r_fft[i+j*L][1] + h_fft[i+j*L][1] * r_fft[i+j*L][0];
    }
  }
  fftw_free(r_fft);
  fftw_free(h_fft);

// Compute moments
  message("Compute spectral moments.");
  vector<vector<double> > moment;
  moment.resize(6);
  for (int i = 0; i < 6; i++)
    moment[i].resize(6);


  for (int i = 0; i < k2; i++) {	
    double i0 = pow(double(i),0), i2 = pow(double(i),2), i4 = pow(double(i),4), im0 = pow(-double(i),0), im2 = pow(-double(i),2), im4 = pow(-double(i),4);
    for (int j = 0; j < k2; j++) {
      double Phi = z_fft[i+j*L][0]*z_fft[i+j*L][0]           + z_fft[i+j*L][1]*z_fft[i+j*L][1],
	     Phi_10, Phi_01, Phi_11;
      if (i==0 && j==0) {
	Phi_10 = 0;
	Phi_11 = 0;
	Phi_01 = 0;
      }
      else if (i==0) {
	Phi_10 = 0;
	Phi_11 = 0;
	Phi_01 = z_fft[i+(L-j)*L][0]*z_fft[i+(L-j)*L][0]     + z_fft[i+(L-j)*L][1]*z_fft[i+(L-j)*L][1];
      }
      else if (j==0) {
	Phi_10 = z_fft[L-i+j*L][0]*z_fft[L-i+j*L][0]         + z_fft[L-i+j*L][1]*z_fft[L-i+j*L][1]; 
	Phi_11 = 0;
	Phi_01 = 0;
      }
      else { 
	Phi_10 = z_fft[L-i+j*L][0]*z_fft[L-i+j*L][0]         + z_fft[L-i+j*L][1]*z_fft[L-i+j*L][1];
	Phi_11 = z_fft[L-i+(L-j)*L][0]*z_fft[L-i+(L-j)*L][0] + z_fft[L-i+(L-j)*L][1]*z_fft[L-i+(L-j)*L][1];
	Phi_01 = z_fft[i+(L-j)*L][0]*z_fft[i+(L-j)*L][0]     + z_fft[i+(L-j)*L][1]*z_fft[i+(L-j)*L][1];
      }
      double j0 = 1 /* ~ pow(double(j),0)*/, j2 = pow(double(j),2), j4 = pow(double(j),4), jm0 = pow(-double(j),0), jm2 = pow(-double(j),2), jm4 = pow(-double(j),4); 
      moment[0][0] += i0*j0*Phi + im0*j0*Phi_10 + i0*jm0*Phi_01 + im0*jm0*Phi_11;
      moment[0][2] += i0*j2*Phi + im0*j2*Phi_10 + i0*jm2*Phi_01 + im0*jm2*Phi_11;
      moment[2][0] += i2*j0*Phi + im2*j0*Phi_10 + i2*jm0*Phi_01 + im2*jm0*Phi_11;
      moment[2][2] += i2*j2*Phi + im2*j2*Phi_10 + i2*jm2*Phi_01 + im2*jm2*Phi_11;
      moment[0][4] += i0*j4*Phi + im0*j4*Phi_10 + i0*jm4*Phi_01 + im0*jm4*Phi_11;
      moment[4][0] += i4*j0*Phi + im4*j0*Phi_10 + i4*jm0*Phi_01 + im4*jm0*Phi_11;
    }
  }
  double mean_m4 = (3*moment[2][2]+moment[0][4]+moment[4][0])/3., 
         mean_m2 = 0.5*(moment[0][2]+moment[2][0]); 

// Inverse the random surface <z>
  message("Inverse transform of the computed surface.");
  z     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L * L);
  p = fftw_plan_dft_2d(L, L, z_fft, z, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);

// Compute surface gradient / Compute RMS height/ Output PSD
  double* zr, *gradz;
  zr = new double[L*L];
  gradz = new double[L*L];
  double rms_h=0;
  zmax=-1e5;
  zmin=1e5;
  //out.open(("psd"+complement.str()).c_str(),ios::out);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      zr[i+L*j] = z[i+j*L][0]; // + z[i+j*L][1];
      if (zr[i+L*j] > zmax)
	zmax = zr[i+L*j];
      else if (zr[i+L*j] < zmin)
	zmin = zr[i+L*j];
      rms_h += zr[i+j*L]*zr[i+j*L]; 
    }
  }
  rms_h = sqrt(rms_h)/L;
// Scale z to match rms_h
  double scalef = rms/rms_h;
  for (int i = 0; i < L; i++) 
    for (int j = 0; j < L; j++) 
       zr[i+L*j] *= scalef;	
  zmin*=scalef;
  zmax*=scalef;
// PDF heights
  vector<int> pdf_surface;
  pdf_surface.resize(Npdf);
  zmin = -5*rms; 
  zmax = -zmin;
  int norm_pdf_h=0;
  message("Computing PDF(heights)");
  for (int i = 0; i < L2; i++) {
      int bin = int((Npdf-1)*(zr[i]-zmin)/(zmax-zmin));
      if (bin >= 0 && bin < Npdf) {
        pdf_surface[bin] += 1;
	norm_pdf_h++;
    }
  }
// Joint PDF(h,grad)
  out.open(("Pdf_heights"+complement.str()).c_str(),ios::out);
  out << "# [1] h/rms [2] P(h/rms)" << endl;
  out << "# rms = " << rms << ", sqrt(m_02+m_20) = " << sqrt(moment[0][2]+moment[2][0]) <<  endl; 
  for (int i = 0; i < Npdf; i++)
    out << (zmin + i*(zmax-zmin)/Npdf)/rms << " " << Npdf*rms*pdf_surface[i]/double(norm_pdf_h) << endl; 
  out.close(); 
  
// Output surface
  message("Writing surface.");
  out.open(("surface"+complement.str()).c_str(),ios::out);
  out << "# k1=" <<k1<<", k2="<<k2<<", H="<<H<<", L="<<L<<", seed="<<seed<<",rms="<<rms<<endl;
  out << "# mean_m2 = " << mean_m2 << ", mean_m4 = " << mean_m4 << endl; 
  out << "# [1] x [2] y [3] z [4] FFT(z).real [5] FFT(z).imaginary" << endl;
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      out << i*dx << " " << j*dx << " " << zr[i+j*L] << " " << z_fft[i+j*L][0] << " " << z_fft[i+j*L][1] << endl; 
    }
    out << endl;
  }
  out.close();

// Free memory
  fftw_destroy_plan(p);
  fftw_free(z); 
  fftw_free(z_fft); 
  delete [] zr;
  delete [] gradz;

  message("Successfully finished.");
  return 0;
}
