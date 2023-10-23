// Programma sviluppato nel giugno 2003. Primo up date Marzo 2004
// Nel programma si fa uso di una FUNZIONE SENO IPERBOLICO fatta "in casa" perchè la
// funzione intrinseca sinh del compilatore CVI versione 6.0 che normalmente viene usato  
// risulta dare valori sbagliati quando l'argomento della funzione supera in modulo il valore di |36|. 
// IL valore risulta il doppio del valore vero
//#include <analysis.h>
#ifndef FUNCTION_STEFAN
#define FUNCTION_STEFAN
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <ansi_c.h>
#include "structure.h"
#include <string.h>

#define NBB 6000    
//#define NBO 100
#define N_D	2000     
//#define LOG 

extern struct DataInput di; //data input
extern struct DataInput did1; //data input with differential absorption added for layer 1
extern struct DataInput did2; //data input with differential absorption added for layer 2
extern struct DataInput didd1; //data input with second differential absorption added for layer 1
extern struct DataInput didd2; //data input with second differential absorption added for layer 2
extern struct DataInput dis1; //data input with differential scattering added for layer 1
extern struct DataInput dis2; //data input with differential scattering added for layer 2
extern struct DataInput did1d2; //data input with differential absorption added for layer 1 and layer 2
extern struct DataInput did1s1; //data input with differential absorption added for layer 1 and differential scattering added for layer 1
extern struct DataInput did2s1; //data input with differential absorption added for layer 2 and differential scattering added for layer 1
extern struct DataInput did1s2; //data input with differential absorption added for layer 1 and differential scattering added for layer 2
extern struct DataInput did2s2; //data input with differential absorption added for layer 2 and differential scattering added for layer 2


extern struct DataOutput_raw d_o;
extern struct DataOutput_raw d_o_d1; //data output for differential absorption added for layer 1
extern struct DataOutput_raw d_o_d2; //data output for differential absorption added for layer 2
extern struct DataOutput_raw d_o_dd1; //data output for second differential absorption added for layer 1
extern struct DataOutput_raw d_o_dd2; //data output for second differential absorption added for layer 2
extern struct DataOutput_raw d_o_s1; //data output for differential scattering added for layer 1
extern struct DataOutput_raw d_o_s2; //data output for differential scattering added for layer 2
extern struct DataOutput_raw d_o_d1d2; //data output for differential absorption added for layer 1 and layer 2
extern struct DataOutput_raw d_o_d1s1; //data output for differential absorption added for layer 1 and differential scattering added for layer 1
extern struct DataOutput_raw d_o_d2s1; //data output for differential absorption added for layer 2 and differential scattering added for layer 1
extern struct DataOutput_raw d_o_d1s2; //data output for differential absorption added for layer 1 and differential scattering added for layer 2
extern struct DataOutput_raw d_o_d2s2; //data output for differential absorption added for layer 2 and differential scattering added for layer 2


//extern struct DataCurrent dc;
extern struct Data_plot dp;
extern struct Data_plot dpd1; //data plot for differential absorption added for layer 1
extern struct Data_plot dpd2; //data plot for differential absorption added for layer 1
extern struct Data_plot dpdd1; //data plot for second differential absorption added for layer 1
extern struct Data_plot dpdd2; //data plot for second differential absorption added for layer 2
extern struct Data_plot dps1; //data plot for differential scattering added for layer 1
extern struct Data_plot dps2; //data plot for differential scattering added for layer 2
extern struct Data_plot dpd1d2; //data output for differential absorption added for layer 1 and layer 2
extern struct Data_plot dpd1s1; //data output for differential absorption added for layer 1 and differential scattering added for layer 1
extern struct Data_plot dpd2s1; //data output for differential absorption added for layer 2 and differential scattering added for layer 1
extern struct Data_plot dpd1s2; //data output for differential absorption added for layer 1 and differential scattering added for layer 2
extern struct Data_plot dpd2s2; //data plot for differential absorption added for layer 2 and differential scattering added for layer 2

extern struct Data_plot_Raman dpR;
extern struct Data_plot_Raman dpRd1;
extern struct Data_plot_Raman dpRd2;
extern struct Data_plot_Raman dpRds1;
extern struct Data_plot_Raman dpRds2;

void inte(struct DataInput);
void integ(struct DataInput *);  
double dAA_cyl (double );
static void ftd_cyl_2_r_r(double , double ,double *,double *,struct DataInput);
static void ftd_cyl_2_r_i(double , double ,double *,double *,struct DataInput);
static void ftd_cyl_2_i_r(double , double ,double *,double *,struct DataInput);
static double ft_cyl_2_r_i(double , double ,struct DataInput);
static double ft_cyl_2_r_r(double , double ,struct DataInput);
static double ft_cyl_2_i_r(double , double ,struct DataInput);
double r_cyl_2_r_r (double ,double ,double ,struct DataInput);
double r_cyl_2_r_i (double ,double ,double ,struct DataInput);
double r_cyl_2_i_r (double ,double ,double ,struct DataInput);
//double t_cyl_2_r_r (double ,double ,double ,struct DataInput);
//double t_cyl_2_r_i (double ,double ,double ,struct DataInput);
//double t_cyl_2_i_r (double ,double ,double ,struct DataInput);
double seno_h (double);

double rtsafe_cyl_2(void (*funcd)(double, double, double *, double *, struct DataInput), double , double ,
	double, double, struct DataInput );
//float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2,
//	float xacc);

void zbrak_cyl_2(double (*f)(double,double,struct DataInput), double , double , int , double *,
	double *, int *, struct DataInput,double);
//void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
//	float xb2[], int *nb);

void f_roots(struct DataInput,struct DataOutput_raw *); 
void f_plot(struct DataInput ,struct DataOutput_raw *, struct Data_plot *,int );
void f_plot_Raman(struct DataInput, struct DataInput, struct DataInput, struct Data_plot *, struct Data_plot *, struct Data_plot *, struct Data_plot_Raman *, int);
 
double bessj1(double);
double bessj0(double);  

int  lettura_riga_chiave(FILE * ,char *, char *); 
void shift_left(char *, int);
void spazi(char *str);
int  read_ascii_line(FILE * , char *,int); 

void selectionSort(int, double *);
void copyStruct(struct DataInput, struct DataInput *);
void Jacobian(int);

double *dfdu1, *dfds1, *dfdu2, *dfds2, *dFidu1, *dFidu2;
#endif
