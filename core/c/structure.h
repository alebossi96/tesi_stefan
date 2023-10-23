#ifndef STRUCTURE_H
#define STRUCTURE_H

#define N_REC_MAX	     7
#define N_PUNTI_TPSF_MAX 400


struct DataInput {
	// Structure of data input for the program cyl_lay2_turbo C version
	// up dated on june 2003
	// Modification with Raman scattering, updated July 2020
	//
	double z0;     // thickness of the first layer (mm)
	double z1;     // thickness of the second layer (mm)
	double rec[N_REC_MAX];  // set of source-receiver distances
	double ro;     // point where the TPSF is evaluated (distance from the laser beam mm)
	double R;      // lateral dimention of the cylinder (Radius of the cylinder mm)
	double ua0;    // absorption of the first layer (mm^-1) (only background uab0)
	double ua1;    // absorption of the second layer (mm^-1) (only background uab1)
	double ud0;    // reduced scattering coefficient of the first layer (mm^-1) 
	double ud1;    // reduced scattering coefficient of the second layer (mm^-1)
	double usR0;    // Raman scattering coefficient of the first layer (mm^-1)
	double usR1;    // Raman scattering coefficient of the second layer (mm^-1)
	double n0;     // refractive index of the first layer 
	double n1;     // refractive index of the second layer 
	double ne;     // refractive index of the external medium 
	double c;      // speed of light in air (mm/ps)
	double v0;     // speed of light in the first layer (mm/ps)
	double v1;     // speed of light in the second layer (mm/ps) 
	double n0e;   // relative refractive index first layer - external medium
	double n1e;   // relative refractive index second layer - external medium
	double an01;   // relative refractive index first-second layer
	double A0e;    // A factor between first layer and external
	double A1e;	 // A factor between second layer and external
	double dt;     // time step of the TPSF  (ps)
	double tmin;	 // Startign time of the TPSF  (ps)
	double dr;
	double rmin;
	double xacc;   // accuracy of the roots calculated with the routine rtsafe
	double precisione; // accuracy
	double pi;     // 2*asin(1)
	int    i_type;   // reflectance (1) or transmittance (2)
	int    n_rec;   // number of receiver
	int    n_tpsf;  // number of steps in the TPSF
	int    n_cw;
	int    ni;     // number of roots kn0 and kn1 
	int    n_kj;   // number of term kj 
	int    nh;     // partition number of the interval for real roots 
	int    nhi;    // partition number of the interval for immaginary roots
};

struct DataOutput_raw {
	// Structure of data output for the program cyl_lay2_turbo C version
	// up dated on june 2003
	//

	double *kappa_z0;     // vettore dove vengono immagazzinate le radici kz0
	double *kappa_j;		// vettore dove vengono immagazzianate le radici kj
	int    *i_bound;      // vettore che contiene il numero di radici kz0 per ogni kj

};

struct Data_plot {
	// Structure of data output for the program cyl_lay2_turbo C version
	// up dated on june 2003
	//

	double r_tpsf[N_PUNTI_TPSF_MAX][N_REC_MAX];     // vettore dove vengono immagazzinate la riflettanza risolta in tempo
	double t_tpsf[N_PUNTI_TPSF_MAX][N_REC_MAX];     // vettore dove vengono immagazzinate la trasmittanza risolta in tempo
	//double r_tpsf_e[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the fluence (at Raman emission wavelength) vs time is stored for the reflectance configuration
	//double r_tpsf_dua1[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the fluence vs time is stored for the reflectance configuration, but with different absorption ua1
	//double r_tpsf_dua2[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the fluence vs time is stored for the reflectance configuration, but with different absorption ua2
	//double t1avg[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the average time spent by detected light in the first layer is stored
	//double t2avg[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the average time spent by detected light in the second layer is stored
  //  double r_cw[1][N_REC_MAX];     // vettore dove vengono immagazzinate la riflettanza cw
  //  double t_cw[1][N_REC_MAX];     // vettore dove vengono immagazzinate la trasmittanza cw  
	double dt;              // larghezza temporale di uno step
	double tmin;              // tempo iniziale
	double dr;              // larghezza spaziale di uno step per la continua
	double rmin;              // distanza sorgente- ricevitore minima
	int    n_tpsf;      // numero di punti delle tpsf
	int    n_cw;        // numero di punti della continua


};

struct Data_plot_Raman
{
	double r_tpsf_e[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the fluence (at Raman emission wavelength) vs time is stored for the reflectance configuration
	double t1avg[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the average time spent by detected light in the first layer is stored
	double t2avg[N_PUNTI_TPSF_MAX][N_REC_MAX]; //array where the average time spent by detected light in the second layer is stored
};

//struct DataCurrent {
// Structure of data current for the program cyl_lay2_turbo C version
// up dated on june 2003
//

//  double x1;     // smaller bound value for the hunting of kz0
//  double x2;	 // larger bound value for the hunting of kz0
//  double kz0;    // eigenvalue for the z direction
//  double kj;     // eigenvalue for the radial direction
//  double t;      // time
//  double square; // square parameter

//};
#endif
