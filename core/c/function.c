#include "function.h"

void f_roots(struct DataInput di,struct DataOutput_raw *d_o)
{
//  Funzione che calcola le radici dell'equazione trascendente
   double square,ref,kz0min,kz1min,d0_0,d1_1,*xb1,*xb2,x1,x2;
   double zb1_r_i[NBB],zb1_i_r[NBB], zb2_r_i[NBB],zb2_i_r[NBB];
   double fx1,fx2,dfx1,dfx2,pi_2,base;
   double discontinuita_1[N_D],discontinuita_2[N_D],*discontinuita_t; 
   int i,ii,j,jj,nzz_r_i,nzz_i_r,M0, nb=NBB, index, i_r_i, i_i_r, i_r_r;
   int iden,idisc_t,ibis,it,status;
   discontinuita_t = NULL; 

   xb1=calloc(N_D, sizeof(double));
   xb2=calloc(N_D, sizeof(double));     
   idisc_t=0;
   pi_2=di.pi/2.;
   d0_0=1./(3*di.ud0);
   d1_1=1./(3*di.ud1);   

   for (ii=0;ii<di.n_kj;ii++)
   {
     square=pow(d_o->kappa_j[ii],2.);  // kappa_j(ii)**2
	 ref=square*(d1_1/di.n1-d0_0/di.n0)+(di.ua1/di.n1-di.ua0/di.n0);  // (kappa_j(ii)**2)*(D1/n1-D0/n0)+ua1/n1-ua0/n0
//     Rispetto alla teoria del manoscritto il termine ref e' uguale a -C*D1/n1
//     in the following if control we are ckecking which are the smallest
//     values of the k coefficients 
	 if (ref < 0.)
	 {
	   kz1min=sqrt(-ref/d1_1*di.n1);
	   kz0min=di.precisione;
	 }
     else if (ref > 0.)
     {
	   kz0min=sqrt(ref/d0_0*di.n0);
	   kz1min=di.precisione;
	 }
     else
     {
	   kz0min=di.precisione;
	   kz1min=di.precisione;
	 }
	 nzz_r_i=0;
	 nzz_i_r=0; 	  
//-----------------------------------------------------------*
//     Bracketing of immaginary roots	
     if (ref >= 0.)
     {                            // kz0 real, kz1 immaginary
        M0=(int)((sqrt(ref/d0_0*di.n0))*(di.z0+2*di.A0e*d0_0)/di.pi);
        if (M0 > 0 )
        {
           x1=M0*di.pi/((di.z0+2*di.A0e*d0_0))+di.precisione;
	       x2=sqrt(ref/d0_0*di.n0)-di.precisione;
           nb=6000;    
           zbrak_cyl_2(ft_cyl_2_r_i,x1,x2,di.nh,xb1,xb2,&nb,di,square);
           if (nb > 0)
           {
             for(j=0; j< nb; j++)
             {
               nzz_r_i=nzz_r_i+1; 
               zb1_r_i[nzz_r_i-1]=xb1[j];
               zb2_r_i[nzz_r_i-1]=xb2[j];
             }
           }  
	       for (i=1; i<=M0; i++)
	       {
              x1=(2*i-1)*di.pi/(2*(di.z0+2*di.A0e*d0_0));
	          x2=(2*i)*di.pi/(2*(di.z0+2*di.A0e*d0_0));
              nzz_r_i=nzz_r_i+1;	          
              zb1_r_i[nzz_r_i-1]=x1+di.precisione;
              zb2_r_i[nzz_r_i-1]=x2-di.precisione;
	       }
	    }
        else if (M0 == 0 )
        {
           if (kz0min*(di.z0+2*di.A0e*d0_0) > pi_2)
           {
              x1=(pi_2)/(di.z0+2*di.A0e*d0_0)+di.precisione ;
              x2=sqrt(ref/d0_0*di.n0)-di.precisione;
              nb=6000;    
              zbrak_cyl_2(ft_cyl_2_r_i,x1,x2,di.nh,xb1,xb2,&nb,di,square);
              if (nb > 0)
              {
                 for (j=0; j<nb; j++)
                 {
                    nzz_r_i=nzz_r_i+1; 
                    zb1_r_i[nzz_r_i-1]=xb1[j];
                    zb2_r_i[nzz_r_i-1]=xb2[j];
                 }
              }
           }
        }
      }
      else if (ref < 0.)    // kz0 immaginary, kz1 real
      {
	    M0=(int)((sqrt(-ref/d1_1*di.n1))*(di.z1+2*di.A1e*d1_1)/di.pi);
	    if (M0 > 0)
	    {
          x1=di.precisione;
	      x2=sqrt(-pow((M0*di.pi/((di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*(di.n0/di.n1)
            -ref/d0_0*di.n0)-di.precisione;
          nb=6000;    
          zbrak_cyl_2(ft_cyl_2_i_r,x1,x2,di.nh,xb1,xb2,&nb,di,square);     
          if (nb > 0)
          {
            for (j=0; j<nb; j++)
		    {
		      nzz_i_r=nzz_i_r+1; 		    
              zb1_i_r[nzz_i_r-1]=xb1[j];
              zb2_i_r[nzz_i_r-1]=xb2[j];
            }
          }
          for (i=1; i<=M0; i++)
          {
            x1=sqrt(-pow(((2*i)*di.pi/(2*(di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*di.n0/di.n1
               -ref/d0_0*di.n0)+di.precisione;
	        x2=sqrt(-pow(((2*i-1)*di.pi/(2*(di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*di.n0/di.n1
               -ref/d0_0*di.n0)-di.precisione;
	    	nzz_i_r=nzz_i_r+1;
           zb1_i_r[nzz_i_r-1]=x1+di.precisione;
           zb2_i_r[nzz_i_r-1]=x2-di.precisione;
		  }
        } 
	    else if (M0 == 0)
	    {
	      if (kz1min*(di.z1+2*di.A1e*d1_1) > pi_2)
	      {
            x1=di.precisione;
            x2=sqrt(-pow((di.pi/(2*(di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*di.n0/di.n1
               -ref/d0_0*di.n0)-di.precisione;
            nb=6000;    
            zbrak_cyl_2(ft_cyl_2_i_r,x1,x2,di.nh,xb1,xb2,&nb,di,square);
            if (nb > 0)
            {
              for (j=0; j<nb; j++)
              {
		        nzz_i_r=nzz_i_r+1;               
                zb1_i_r[nzz_i_r-1]=xb1[j];
                zb2_i_r[nzz_i_r-1]=xb2[j];
              }
            }
	      }			 			 	
	    }
	  }
//-----------------------------------------------------------*   
//-----------------------------------------------------------
//     Calcolo delle radici immaginarie Caso ref>0 ! kz0 real, kz1 immaginary
      if (ref > 0.)
      {
        i_r_i=0;
        for(i=0; i<nzz_r_i; i++)
        {
          x1=zb1_r_i[i];
          x2=zb2_r_i[i];
          ftd_cyl_2_r_i(square,x1,&fx1,&dfx1,di);
          if ( fx1 > 0. && dfx1 < 0.)
          {
	        index=(di.n_kj*3)*i_r_i+3*ii+0;
            d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_i,x1,x2,di.xacc,square,di);
     	    i_r_i=i_r_i+1;
          }
          else if (fx1 < 0. && dfx1 > 0.)
          {
	        index=(di.n_kj*3)*i_r_i+3*ii+0; 
            d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_i,x1,x2,di.xacc,square,di);
	        i_r_i=i_r_i+1;      
          }      
        }
        index=3*ii+0;
	    d_o->i_bound[index]=i_r_i;
	  }
//----------------------------------------------------------------
//     Calcolo delle radici immaginarie Caso ref<0 ! kz0 immaginary, kz1 real
	  else if (ref < 0.)
	  {
        i_i_r=0;
        for (i=0; i<nzz_i_r; i++)
        {
          x1=zb1_i_r[i];
          x2=zb2_i_r[i];        
          ftd_cyl_2_i_r(square,x1,&fx1,&dfx1,di);
          if (fx1 > 0. && dfx1 < 0.)
          {
            index=(di.n_kj*3)*i_i_r+3*ii+1;
            d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_i_r,x1,x2,di.xacc,square,di);            
            i_i_r=i_i_r+1;
          }
          else if (fx1 < 0. && dfx1 > 0.)
          {
            index=(di.n_kj*3)*i_i_r+3*ii+1;
            d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_i_r,x1,x2,di.xacc,square,di);
            i_i_r=i_i_r+1;
          }
        }
        index=3*ii+1;
	    d_o->i_bound[index]=i_i_r;        
	  }
//---------------------------------------------------------------
      i_r_r=0; 
//-----------------------------------------------------------*
//     classificazione delle discontinuita' per radici reali
	  jj=0;
	  for (i=1; i <= di.ni; i++)
	  {
	    if (((2*i-1)*pi_2/(di.z0+2*di.A0e/(3*di.ud0))) > kz0min)
	    {
	      jj=jj+1;  
	      discontinuita_1[jj-1]=(2*i-1)*pi_2/(di.z0+2*di.A0e/(3*di.ud0));
	    }
	  }
      iden=0;
      i=1;
	  base=0.;  
	  while (sqrt(fabs(base)) < discontinuita_1[jj-1])
		  {
		    discontinuita_2[i-1]=(2*i-1)*pi_2/(di.z1+2*di.A1e/(3*di.ud1));
            base=((ref*di.n0+pow(discontinuita_2[i-1],2.)*d1_1*di.n0/di.n1)/d0_0);
	        if (base > 0.)
	        {
	          if (sqrt(base) > kz0min)
	          {
	          	iden=iden+1;
	            discontinuita_2[iden-1]=sqrt(base);
	          }
	        }
	        i=i+1;
		  }
      it=jj+iden;	  
//      if (idisc_t == 1)
//      {
//        free(discontinuita_t);    
//	  }
//	  discontinuita_t=calloc(it,sizeof(double));
//      idisc_t=1;
      
	  if (idisc_t == 0)
		  {
	         discontinuita_t=calloc(it,sizeof(double));
             idisc_t=1;		  
		  }
	  else if(idisc_t == 1)
		  {
             free(discontinuita_t);
             discontinuita_t=calloc(it,sizeof(double));             
		  }
      
      		  
	  for (i=0; i<jj; i++)
		  {
		    discontinuita_t[i]=discontinuita_1[i];
		  }
	  for (i=0; i<iden; i++)
		  {
		    discontinuita_t[jj+i]=discontinuita_2[i];
		  }      	  
	  selectionSort(it, discontinuita_t);
// Questa routine può essere sostituita dalla routine Sort del numerical recipes in C
// che ho riportato in coda al programma. Prima di usarla c'è da aggiornare i dimentionamenti
// delle variabili da float a double
//-----------------------------------------------------------------	  
//     Calcolo delle radici reali kz0 e kz1
      nb=6000; 
      x1=kz0min+di.precisione;
      x2=discontinuita_t[0]-di.precisione;	 
      zbrak_cyl_2(ft_cyl_2_r_r,x1,x2,di.nh,xb1,xb2,&nb,di,square); 
	  if (nb > 0)
	  {
        for (i=0; i<nb; i++)
        {
          x1=xb1[i];
          x2=xb2[i];
          ftd_cyl_2_r_r(square,x1,&fx1,&dfx1,di);
		  if (fx1 > 0. && dfx1 < 0.)
			  {
                index=(di.n_kj*3)*i_r_r+3*ii+2;
                d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di);               
                i_r_r=i_r_r+1; 			  
			  }
		  else if (fx1 < 0. && dfx1 > 0.)
			  {
                index=(di.n_kj*3)*i_r_r+3*ii+2;
                d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di);               
                i_r_r=i_r_r+1;			  
			  }
        }  
	  }
	  i=1; 
	  while (i_r_r < di.ni && i <= it-1)
		  {
	        x1=discontinuita_t[i-1]+di.precisione;  
            x2=discontinuita_t[i]-di.precisione;
            ftd_cyl_2_r_r(square,x1,&fx1,&dfx1,di);
            ftd_cyl_2_r_r(square,x2,&fx2,&dfx2,di);
			if (fx1*fx2 < 0.)
				{
                index=(di.n_kj*3)*i_r_r+3*ii+2;
                d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di);               
                i_r_r=i_r_r+1;				
				}
			else if (fx1*fx2 > 0.)
				{
                  nb=6000;
                  zbrak_cyl_2(ft_cyl_2_r_r,x1,x2,di.nh,xb1,xb2,&nb,di,square);
				  if (nb > 0)
					  {
					   for (ibis=0; ibis<nb; ibis++)
						   {
                            x1=xb1[ibis];
                            x2=xb2[ibis];
                            ftd_cyl_2_r_r(square,x1,&fx1,&dfx1,di);
							if (fx1 > 0. && dfx1 < 0.)
								{
                                  index=(di.n_kj*3)*i_r_r+3*ii+2;
                                  d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di);               
                                  i_r_r=i_r_r+1;								
								}
							else if (fx1 < 0. && dfx1 > 0.)
								{
                                  index=(di.n_kj*3)*i_r_r+3*ii+2;
                                  d_o->kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di);               
                                  i_r_r=i_r_r+1;								
								}
						   }
					  }
				}
			i=i+1;	
		  }
        index=3*ii+2;
	    d_o->i_bound[index]=i_r_r; 	  
   }
   if (idisc_t == 0)
	   {
         free(discontinuita_t); 	   
	   }
}			 
// fine funzione f_roots
//-----------------------------------------------------------*
void f_plot(struct DataInput di, struct DataOutput_raw *d_o, struct Data_plot *dp,int i_r)
{
 int   i,ii,j,jj,index_a,index_b;
 double t,kj,kz0;
// questa funzione costruisce il vettore del profilo temporale uttilizzando le radici calcolate dalla funzione f_roots

	for (i=0; i<di.n_tpsf; i++)
		{
          t=di.tmin+i*di.dt;
		  for (j=0; j<di.n_kj; j++)  
			  {
			    kj=d_o->kappa_j[j];
			    index_b=j*3+0;
				for (ii=0; ii<d_o->i_bound[index_b]; ii++)
					{
					  index_a=(di.n_kj*3)*ii+3*j+0;
					  kz0=d_o->kappa_z0[index_a];
					  
					  dp->r_tpsf[i][i_r]=dp->r_tpsf[i][i_r]+r_cyl_2_r_i(t,kz0,kj,di);
					  //dp->t_tpsf[i][i_r]=dp->t_tpsf[i][i_r]+t_cyl_2_r_i(t,kz0,kj,di);
					}
			    index_b=j*3+1;
				for (ii=0; ii<d_o->i_bound[index_b]; ii++)
					{
					  index_a=(di.n_kj*3)*ii+3*j+1;
					  kz0=d_o->kappa_z0[index_a];
					  dp->r_tpsf[i][i_r]=dp->r_tpsf[i][i_r]+r_cyl_2_i_r(t,kz0,kj,di);
					  //dp->t_tpsf[i][i_r]=dp->t_tpsf[i][i_r]+t_cyl_2_i_r(t,kz0,kj,di);
					}
				index_b=j*3+2;
				for (ii=0; ii<d_o->i_bound[index_b]; ii++)
					{
					  index_a=(di.n_kj*3)*ii+3*j+2;
					  kz0=d_o->kappa_z0[index_a];
					  dp->r_tpsf[i][i_r]=dp->r_tpsf[i][i_r]+r_cyl_2_r_r(t,kz0,kj,di);
					  //dp->t_tpsf[i][i_r]=dp->t_tpsf[i][i_r]+t_cyl_2_r_r(t,kz0,kj,di);
					}	
			  }
		}
}
// fine funzione f_plot
//-----------------------------------------------------------*

void f_plot_Raman(struct DataInput di, struct DataInput did1, struct DataInput did2, struct Data_plot *dp, struct Data_plot *dpd1, struct Data_plot *dpd2, struct Data_plot_Raman *dpR, int i_r)
{
	for (int i = 0; i < di.n_tpsf; i++)
	{
		//dpR->t1avg[i][i_r] = (dp->r_tpsf[i][i_r] - dpd1->r_tpsf[i][i_r]) / (di.v0*(did1.ua0 - di.ua0)) / dp->r_tpsf[i][i_r];
		//dpR->t2avg[i][i_r] = (dp->r_tpsf[i][i_r] - dpd2->r_tpsf[i][i_r]) / (di.v1*(did2.ua1 - di.ua1)) / dp->r_tpsf[i][i_r];
		dpR->t1avg[i][i_r] = (1 - dpd1->r_tpsf[i][i_r] / dp->r_tpsf[i][i_r]) / (di.v0*(did1.ua0 - di.ua0));
		dpR->t2avg[i][i_r] = (1 - dpd2->r_tpsf[i][i_r] / dp->r_tpsf[i][i_r]) / (di.v1*(did2.ua1 - di.ua1));
		dpR->r_tpsf_e[i][i_r] = dp->r_tpsf[i][i_r] * (di.usR0*di.v0*dpR->t1avg[i][i_r] + di.usR1*di.v1*dpR->t2avg[i][i_r]);
	}
	
}

void Jacobian(int i_r)
{
	for (int i = 0; i < di.n_tpsf; i++)
	{
		dfdu1[i] = -di.usR0*(dpdd1.r_tpsf[i][i_r]-2*dpd1.r_tpsf[i][i_r]+dp.r_tpsf[i][i_r]) / ((did1.ua0 - di.ua0)*(did1.ua0 - di.ua0))
			-di.usR1*(dpd1d2.r_tpsf[i][i_r]-dpd1.r_tpsf[i][i_r]-dpd2.r_tpsf[i][i_r]+dp.r_tpsf[i][i_r])/((did1.ua0 - di.ua0)*(did2.ua1-di.ua1));

		dfds1[i] = -di.usR0*(dpd1s1.r_tpsf[i][i_r] - dpd1.r_tpsf[i][i_r] - dps1.r_tpsf[i][i_r] + dp.r_tpsf[i][i_r]) / ((did1.ua0 - di.ua0)*(dis1.ud0 - di.ud0))
			- di.usR1*(dpd2s1.r_tpsf[i][i_r] - dpd2.r_tpsf[i][i_r] - dps1.r_tpsf[i][i_r] + dp.r_tpsf[i][i_r]) / ((did2.ua1-di.ua1)*(dis1.ud0-di.ud0));

		dfdu2[i] = -di.usR0*(dpd1d2.r_tpsf[i][i_r] - dpd1.r_tpsf[i][i_r] - dpd2.r_tpsf[i][i_r] + dp.r_tpsf[i][i_r]) / ((did1.ua0 - di.ua0)*(did2.ua1 - di.ua1))
			- di.usR1*(dpdd2.r_tpsf[i][i_r] - 2*dpd2.r_tpsf[i][i_r] + dp.r_tpsf[i][i_r]) / ((did2.ua1 - di.ua1)*(did2.ua1 - di.ua1));

		dfds2[i] = -di.usR0*(dpd1s2.r_tpsf[i][i_r] - dpd1.r_tpsf[i][i_r] - dps2.r_tpsf[i][i_r] + dp.r_tpsf[i][i_r]) / ((did1.ua0 - di.ua0)*(dis2.ud1 - di.ud1))
			- di.usR1*(dpd2s2.r_tpsf[i][i_r] - dpd2.r_tpsf[i][i_r] - dps2.r_tpsf[i][i_r] + dp.r_tpsf[i][i_r]) / ((did2.ua1 - di.ua1)*(dis2.ud1 - di.ud1));

		dFidu1[i] = (dp.r_tpsf[i][i_r] - dpd1.r_tpsf[i][i_r]) / (did1.ua0 - di.ua0);

		dFidu1[i] = (dp.r_tpsf[i][i_r] - dpd2.r_tpsf[i][i_r]) / (did2.ua1 - di.ua1);
	}
}


      double dAA_cyl (double an12)
{      
      double A;
      if (an12 > 1.)
        A=504.332889-2641.00214*an12+5923.699064*an12*an12
        +-7376.355814*an12*an12*an12+5507.5304*pow(an12,4)-2463.357945
     *pow(an12,5.)+610.956547*pow(an12,6.)-64.8047*pow(an12,7.);
       else if (an12 == 1.)
        A=1.;
       else 
        A=3.084635-6.531194*an12+8.357854*an12*an12-5.082751
     *an12*an12*an12+1.171382*pow(an12,4.);
      return A;
}      
//-----------------------------------------------------------*
static void ftd_cyl_2_r_r(double square, double kz0,double *f,double *df,struct DataInput di)
// Funzione che calcola la funzione e la derivata dell'equazione trascendente
// Caso considerato = Reale - Reale
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
{
      double ga1,dga1,kz1,dkz1,ga0,dga0,d0,d1;      
      kz1=sqrt((pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+
       (3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.)));
      dkz1=kz0/kz1*(di.ud1/di.ud0)*(di.n1/di.n0);	
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);
      ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*d1*kz1);
      ga0=(2*di.A0e*d0*kz0);
      dga0=2*di.A0e*d0; 
      dga1=-dkz1*(di.z0+di.z1)-2*di.A1e*d1*dkz1;
      *f=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-
         (3*di.ud1/kz1)*tan(kz1*di.z0+ga1)*pow((di.n0/di.n1),2.);     
      *df=-(3*di.ud0/pow(kz0,2.))*tan(kz0*di.z0+ga0)+
          (3*di.ud1/pow(kz1,2.))*tan(kz1*di.z0+ga1)*dkz1*pow((di.n0/di.n1),2.)+
          (3*di.ud0/kz0)/pow((cos(kz0*di.z0+ga0)),2.)*
          (di.z0+dga0)-(3*di.ud1/kz1)/pow((cos(kz1*di.z0+ga1)),2.)* 
          (dkz1*di.z0+dga1)*pow((di.n0/di.n1),2.);     
}      
//-----------------------------------------------------------*
static void ftd_cyl_2_r_i(double square, double kz0,double *f,double *df,struct DataInput di)
// Funzione che calcola la funzione e la derivata dell'equazione trascendente
// Caso considerato = Reale - Immaginario
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
//     kz1 is immaginary
//     f and df are the values of the function and of his derivative with respect kz0
{
      double ga1,dga1,mod_kz1,dmod_kz1,kz1_2,ga0,dga0,d0,d1,arg,darg;
      kz1_2=(pow(kz0,2)*(di.ud1/di.ud0)*(di.n1/di.n0)+
      (3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.));
      mod_kz1=sqrt(fabs(-kz1_2));
      dmod_kz1=-kz0/mod_kz1*(di.ud1/di.ud0)*(di.n1/di.n0);      
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);      
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1);
      ga0=(2*di.A0e*d0*kz0);
      dga0=2*di.A0e*d0; 
      dga1=-dmod_kz1*(di.z0+di.z1)-2*di.A1e*d1*dmod_kz1;
      arg=mod_kz1*di.z0+ga1;
      darg=dmod_kz1*di.z0+dga1;	   
      *f=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-
      (3*di.ud1/mod_kz1)*tanh(arg)*pow((di.n0/di.n1),2.);      
      *df=-(3*di.ud0/pow(kz0,2.))*tan(kz0*di.z0+ga0)+
      (3*di.ud1/pow(mod_kz1,2.))*dmod_kz1*tanh(arg)*pow((di.n0/di.n1),2.)+
      (3*di.ud0/kz0)/pow((cos(kz0*di.z0+ga0)),2.)*(di.z0+dga0)-
      (3*di.ud1/mod_kz1)/pow(cosh(arg),2.)*darg*pow((di.n0/di.n1),2.);
}     
//-----------------------------------------------------------*
static void ftd_cyl_2_i_r(double square, double kz0,double *f,double *df,struct DataInput di)
// Funzione che calcola la funzione e la derivata dell'equazione trascendente
// Caso considerato = Reale - Immaginario
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
//     kz0 is immaginary
//     La variabile kz0 rappresenta in questo caso il modulo |kn0| della variabile kn0
//     f and df are the values of the function and of his derivative with respect |kn0|
{
      double ga1,dga1,mod_kz1,dmod_kz1,kz1_2,ga0,dga0,d0,d1,arg,darg; 
      kz1_2=(-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+
      (3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.));
      mod_kz1=sqrt(fabs(kz1_2));
      dmod_kz1=-kz0/mod_kz1*(di.ud1/di.ud0)*(di.n1/di.n0);	      
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);      
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1);
      ga0=(2*di.A0e*d0*kz0);
      dga0=2*di.A0e*d0; 
      dga1=-dmod_kz1*(di.z0+di.z1)-2*di.A1e*d1*dmod_kz1;
      arg=mod_kz1*di.z0+ga1;
      darg=dmod_kz1*di.z0+dga1;	   
      *f=(3*di.ud0/kz0)*tanh(kz0*di.z0+ga0)-
      (3*di.ud1/mod_kz1)*tan(arg)*pow((di.n0/di.n1),2.);     
      *df=-(3*di.ud0/pow(kz0,2.))*tanh(kz0*di.z0+ga0)+
      (3*di.ud1/pow(mod_kz1,2.))*dmod_kz1*tan(arg)*pow((di.n0/di.n1),2.)+
      (3*di.ud0/kz0)/pow((cosh(kz0*di.z0+ga0)),2.)*
      (di.z0+dga0)-
      (3*di.ud1/mod_kz1)/pow(cos(arg),2.)*darg*pow((di.n0/di.n1),2.);     
}
//-----------------------------------------------------------*
static  double ft_cyl_2_r_i (double square, double kz0,struct DataInput di)
// Funzione che calcola la funzione e la derivata dell'equazione trascendente
// Caso considerato = Reale - Immaginario
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
//     kz1 is immaginary
//     f and df are the value of the function with respect kz0
{
      double ga1,dga1,mod_kz1,dmod_kz1,kz1_2,ga0,dga0,d0,d1,arg,f_out;
      kz1_2=(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+
      (3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.));
	  mod_kz1=sqrt(fabs(-kz1_2));   //   valore assoluto inserito per precauzione
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);       
      ga0=(2*di.A0e*d0*kz0);
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1);
      arg=mod_kz1*di.z0+ga1;   
      f_out=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-
      (3*di.ud1/mod_kz1)*tanh(arg)*pow((di.n0/di.n1),2.);
      return f_out;
}       
//-----------------------------------------------------------*
static double ft_cyl_2_r_r (double square, double kz0,struct DataInput di)
// Funzione che calcola la funzione dell'equazione trascendente
// Caso considerato = Reale - Reale
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
{      
      double ga1,dga1,kz1,dkz1,ga0,dga0,d0,d1,f_out;
      kz1=sqrt(fabs(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+
      (3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.)));
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);      
      ga0=(2*di.A0e*d0*kz0);
      ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*d1*kz1);
      f_out=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-
      (3*di.ud1/kz1)*tan(kz1*di.z0+ga1)*pow((di.n0/di.n1),2.);
      return f_out;
}     
//-----------------------------------------------------------*
static double ft_cyl_2_i_r (double square, double kz0,struct DataInput di)
// Funzione che calcola la funzione e la derivata dell'equazione trascendente
// Caso considerato = Reale - Immaginario
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
//     kz0 is immaginary
//     La variabile kz0 rappresenta in questo caso il modulo |kn0| della variabile kn0
//     f and df are the values of the function and of his derivative with respect |kn0|
{
      double ga1,dga1,mod_kz1,dmod_kz1,kz1_2,ga0,dga0,d0,d1,arg,f_out;
      kz1_2=(-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+
      (3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.));
	  mod_kz1=sqrt(fabs(kz1_2));   //   valore assoluto dentro la radice inserito per precauzione
      d0=1./(3*di.ud0);                // per problemi di precisione si possono avere valori piccolissi e negativi di kz1_2
      d1=1./(3*di.ud1);      
      ga0=(2*di.A0e*d0*kz0);
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1);
      arg=mod_kz1*di.z0+ga1;   
      f_out=(3*di.ud0/kz0)*tanh(kz0*di.z0+ga0)-
      (3*di.ud1/mod_kz1)*tan(arg)*pow((di.n0/di.n1),2.);
      return f_out;
}
//-----------------------------------------------------------*
      double r_cyl_2_r_r (double t,double kz0,double kj,struct DataInput di)      
// Funzione che calcola la riflettanza da un mezzo a due strati
// Caso considerato = Reale - Reale
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
{
      double d0,d1,f_out,ga0,j1,j0;
      double kz1,b1,ga1,arg1,Nl;
	  j1=bessj1(di.R*kj);
	  j0=bessj0(di.ro*kj);
///      pi=2.d0*dasin(1.d0)
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);      
      kz1=sqrt(fabs(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+
         (pow(kj,2.))*(di.ud1/di.ud0*di.n1/di.n0-1.)));
      ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*d1*kz1);
      ga0=(2*di.A0e*d0*kz0);       
      arg1=sin(kz1*di.z0+ga1);   
      b1=sin(kz0*di.z0+ga0)/arg1*pow((di.n1/di.n0),2.);
      Nl=di.z0/2.+ga0/(2.*kz0)-pow(b1,2.)/(2*kz1)*(kz1*di.z0+ga1)*(di.n0/di.n1) // Il contributo del termine radiale
     -(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))            // e' inserito nella fil2_cyl
     +(pow(b1,2.)/(4*kz1))*(sin(2*(kz1*di.z0+ga1)))*(di.n0/di.n1);
      f_out=d0*di.v0*kz0*cos(ga0)*sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0 
     /Nl*exp(-di.ua0*di.v0*t)*exp(-(pow(kj,2.)+pow(kz0,2.))*d0*di.v0*t);   
      return f_out;
}
//-----------------------------------------------------------*
      double r_cyl_2_r_i (double t,double kz0,double kj,struct DataInput di)
// Funzione che calcola la riflettanza da un mezzo a due strati
// Caso considerato = Reale - Immaginario
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab      

{
      double d0,d1,f_out,ga0,j1,j0;
      double kz1_2,mod_kz1,b1,ga1,arg1,Nl,vedo,vedo1,vedo2;
	  j1=bessj1(di.R*kj);
	  j0=bessj0(di.ro*kj);
	  
//      pi=2.d0*dasin(1.d0)
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);
      kz1_2=((pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+
      (pow(kj,2.))*(di.ud1/di.ud0*di.n1/di.n0-1.)));
      mod_kz1=sqrt(-kz1_2);
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1);
      ga0=(2*di.A0e*d0*kz0);       
      arg1=(mod_kz1*di.z0+ga1);
//	  if (fabs(2*arg1) > 300.)
//		  {
//		  b1=d0;
//		  }
      
      b1=sin(kz0*di.z0+ga0)/seno_h(arg1)*pow((di.n1/di.n0),2.);
//      b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
      Nl=di.z0/2.+ga0/(2.*kz0)+pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1)   // Il contributo del termine radiale 
     -(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))            // e' inserito direttamente nella fil2_r_i
     -(pow(b1,2.)/(4*mod_kz1))*(seno_h(2*arg1))*(di.n0/di.n1);
//     +-(pow(b1,2.)/(4*mod_kz1))*((exp(-2*arg1)-exp(2*arg1))/2)
      f_out=d0*di.v0*kz0*cos(ga0)*sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))
      *j0/Nl*exp(-di.ua0*di.v0*t)*exp(-(pow(kj,2.)+pow(kz0,2.))*d0*di.v0*t);   
      return   f_out;
}      
//-----------------------------------------------------------*
      double r_cyl_2_i_r (double t,double kz0,double kj,struct DataInput di)      
// Funzione che calcola la riflettanza da un mezzo a due strati
// Caso considerato = Reale - Immaginario
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
{
      double d0,d1,f_out,ga0,j1,j0;
      double kz1_2,mod_kz1,b1,ga1,arg1,Nl;
      
	  j1=bessj1(di.R*kj);
	  j0=bessj0(di.ro*kj);
//      pi=2.d0*dasin(1.d0)
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);      
      kz1_2=((-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+
           (pow(kj,2.))*(di.ud1/di.ud0*di.n1/di.n0-1.)));
      mod_kz1=sqrt(fabs(kz1_2));
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1);
      ga0=(2*di.A0e*d0*kz0);       
      arg1=(mod_kz1*di.z0+ga1);   
      b1=seno_h(kz0*di.z0+ga0)/sin(arg1)*pow((di.n1/di.n0),2.);
//      b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
      Nl=-di.z0/2.-ga0/(2.*kz0)-pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1) // IL contributo del termine radiale e' inserito
     +(1./(4.*kz0))*(seno_h(2*(kz0*di.z0+ga0)))          // direttamente nella fil2_i_r
     +(pow(b1,2.)/(4*mod_kz1))*(sin(2*arg1))*(di.n0/di.n1);
//     +-(b1**2/(4*mod_kz1))*((dexp(-2*arg1)-dexp(2*arg1))/2)
      f_out=d0*di.v0*kz0*cosh(ga0)*seno_h(kz0*3*d0+ga0)/
     (di.pi*pow(di.R,2)*pow(j1,2.))*j0
     /Nl*exp(-di.ua0*di.v0*t)*exp((-pow(kj,2.)+pow(kz0,2.))*d0*di.v0*t);   
      return f_out;
}      
//-----------------------------------------------------------*
//-----------------------------------------------------------*
/*double t_cyl_2_r_r (double t,double kz0,double kj,struct DataInput di)      
// Funzione che calcola la trasmittanza attraverso un mezzo a due strati
// Caso considerato = Reale - Reale
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
{
      double d0,d1,f_out,ga0,j1,j0;
      double kz1,b1,ga1,arg1,Nl;
	  j1=bessj1(di.R*kj);
	  j0=bessj0(di.ro*kj);
//      pi=2.d0*dasin(1.d0)
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);       
      kz1=sqrt(fabs(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+
         (pow(kj,2.))*(di.ud1/di.ud0*(di.n1/di.n0)-1.)));
      ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*(1./(3*di.ud1))*kz1);
      ga0=(2*di.A0e*d0*kz0);       
      arg1=sin(kz1*di.z0+ga1);   
      b1=sin(kz0*di.z0+ga0)/arg1*pow((di.n1/di.n0),2.);
      Nl=di.z0/2.+ga0/(2.*kz0)-pow(b1,2.)/(2*kz1)*(kz1*di.z0+ga1)*(di.n0/di.n1) // Il contributo del termine radiale
      -(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))            // e' inserito nella fil2_cyl
      +(pow(b1,2.)/(4*kz1))*(sin(2*(kz1*di.z0+ga1)))*(di.n0/di.n1);
      f_out=-d1*di.v1*di.n1/di.n0*b1*kz1*cos(kz1*(di.z0+di.z1)+ga1)*
      sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua1*di.v1*t)*
      exp(-(pow(kj,2.)+pow(kz1,2.))*d1*di.v1*t);   
      return f_out;
}*/
//-----------------------------------------------------------*
//-----------------------------------------------------------*
      /*double t_cyl_2_r_i (double t,double kz0,double kj,struct DataInput di)      
// Funzione che calcola la trasmittanza attraverso un mezzo a due strati
// Caso considerato = Reale - Immaginario
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
{
      double d0,d1,f_out,ga0,j1,j0;
      double kz1_2,mod_kz1,b1,ga1,arg1,Nl;
	  j1=bessj1(di.R*kj);
	  j0=bessj0(di.ro*kj);
//      pi=2.d0*dasin(1.d0)
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);      
      kz1_2=((pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+
           (pow(kj,2.))*(di.ud1/di.ud0*(di.n1/di.n0)-1.)));
      mod_kz1=sqrt(-kz1_2);
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*(1./(3*di.ud1))*mod_kz1);
      ga0=(2*di.A0e*d0*kz0);       
      arg1=(mod_kz1*di.z0+ga1);   
      b1=sin(kz0*di.z0+ga0)/seno_h(arg1)*pow((di.n1/di.n0),2.);
//      b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
      Nl=di.z0/2.+ga0/(2.*kz0)+pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1)   // Il contributo del termine radiale 
         -(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))            // e' inserito direttamente nella fil2_r_i
         -(pow(b1,2.)/(4*mod_kz1))*(seno_h(2*arg1))*(di.n0/di.n1);
//     +-(b1**2/(4*mod_kz1))*((dexp(-2*arg1)-dexp(2*arg1))/2)
      f_out=-d1*di.v1*di.n1/di.n0*mod_kz1*b1*cosh(mod_kz1*(di.z0+di.z1)+ga1)*
            sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua1*di.v1*t)*
            exp(-(pow(kj,2.)-pow(mod_kz1,2.))*d1*di.v1*t);
      return f_out;
}*/
//-----------------------------------------------------------*
//-----------------------------------------------------------*
      /*double t_cyl_2_i_r (double t,double kz0,double kj,struct DataInput di)      
// Funzione che calcola la trasmittanza attraverso un mezzo a due strati
// Caso considerato = Immaginario - Reale
// Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
{
      double d0,d1,f_out,ga0,j1,j0;
      double kz1_2,mod_kz1,b1,ga1,arg1,Nl;
	  j1=bessj1(di.R*kj);
	  j0=bessj0(di.ro*kj);
//      pi=2.d0*dasin(1.d0)
      d0=1./(3*di.ud0);
      d1=1./(3*di.ud1);    
      kz1_2=((-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+
           (pow(kj,2.))*(di.ud1/di.ud0*(di.n1/di.n0)-1.)));
      mod_kz1=sqrt(fabs(kz1_2));
      ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*(1./(3*di.ud1))*mod_kz1);
      ga0=(2*di.A0e*d0*kz0);       
      arg1=(mod_kz1*di.z0+ga1);   
      b1=seno_h(kz0*di.z0+ga0)/sin(arg1)*pow((di.n1/di.n0),2.);
//      b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
      Nl=-di.z0/2.-ga0/(2.*kz0)-pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1) // IL contributo del termine radiale e' inserito
         +(1./(4.*kz0))*(seno_h(2*(kz0*di.z0+ga0)))          // direttamente nella fil2_i_r
         +(pow(b1,2.)/(4*mod_kz1))*(sin(2*arg1))*(di.n0/di.n1);
//     +-(b1**2/(4*mod_kz1))*((dexp(-2*arg1)-dexp(2*arg1))/2)
      f_out=-d1*di.v1*di.n1/di.n0*mod_kz1*b1*cos(mod_kz1*(di.z0+di.z1)+ga1)*
            seno_h(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua1*di.v1*t)*
            exp((-pow(kj,2.)-pow(mod_kz1,2.))*d1*di.v1*t);   
      return f_out;
}*/
//-----------------------------------------------------------*
//-----------------------------------------------------------*
      double seno_h (double x)      
// Funzione seno iperbolico
{    double f_out;
//    f_out=tanh(x)*cosh(x);
      f_out=(exp(x)-exp(-x))/2.;
      return f_out;
}
//-----------------------------------------------------------*
double bessj0(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;
   
	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
               +y*(-0.6911147651e-5+y*(0.7621095161e-6
			   -y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}
//--------------------------------------------------------
double bessj1(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;
   
	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
             +y*(0.8449199096e-5+y*(-0.88228987e-6
			 +y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}
//--------------------------------------------------------------
#define MAXIT 3000
double rtsafe_cyl_2(void (*funcd)(double, double, double *, double *, struct DataInput), double x1, double x2,
	double xacc, double square, struct DataInput di)
{
//	void nrerror(char error_text[]);
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;
   
	(*funcd)(square,x1,&fl,&df,di);
	(*funcd)(square,x2,&fh,&df,di);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0));
//		nrerror("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(square,rts,&f,&df,di);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)

|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(square,rts,&f,&df,di);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
//	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}
#undef MAXIT
//--------------------------------------------------------------------
void zbrak_cyl_2(double (*fx)(double , double ,struct DataInput), double x1, double x2, int n, double xb1[],
	double xb2[], int *nb,struct DataInput di, double square)
{
	int nbb,i;
	double x,fp,fc,dx;
   
	nbb=0;
	dx=(x2-x1)/n;
	fp=(*fx)(square,x=x1,di);
	for (i=1;i<=n;i++) {
		fc=(*fx)(square,x +=dx,di);
//	fc=(*fx)(square,x =x1+dx*i,di); 
		if (fc*fp <= 0.0) {
			xb1[++nbb-1]=x-dx;
			xb2[nbb-1]=x;
			if(*nb == nbb) return;
   
		}
		fp=fc;
	}
	*nb = nbb;
}
//-------------------------------------------------------------------------

int lettura_riga_chiave(FILE *file ,char *chiave, char *str)
{															  
	int i,ii,iii,status=0;
	char p[200],str1[200],str2[200];
	
	rewind (file);
	for(;;){
      status= read_ascii_line(file , p, 200);
	  if(status < 0){
	    strcpy(str,"");
	    break;
	  }
	  ii=0;
	  for(i=0;i<strlen(p);i++)if(p[i]== 61)ii=i; // cerca ch "="
	  if(ii > 0){
	     for(i=0;i<ii;i++)str1[i]=p[i];
	     str1[ii]=0;
	     spazi(str1);
	     iii=0;
	     for(i=ii+1;i<strlen(p);i++)str2[iii++]=p[i];
	     str2[iii]=0;
	     spazi(str2);
	     if(!strcmp(str1,chiave)){
	       strcpy(str,str2);
	       return 0;
	     }  
	  }
	}  
    return -1;
}

void spazi(char *str)
{
   int ii=0,i=0;
   if(!strlen(str))return;
   while(str[ii] == 32)ii++;
   if(ii > 0)for(i=0;i<strlen(str);i++)str[i]=str[i+ii];
   ii=strlen(str)-1;
   while(str[ii] == 32)ii--;
   str[ii+1]=0;
}

/*--------------------------------------------------------------------*/

int substr_da_a(char *str, int c1, int c2, char *sstr)
{
   int i,ret;
   
   sstr[0]=0;
   ret=0;
   if(c1 > -1 && c2 > 0 && c2 > c1){
     for(i=c1; i<=c2; i++)sstr[ret++] = str[i];
   }
   sstr[ret]=0;
   return ret;
}

void shift_left(char *str, int c1)
{
   int i,ret,c2;
   
   ret=0;
   c2=strlen(str);
   if(c1 > 0 ){
     for(i=c1 ; i < c2; i++)str[ret++] = str[i];
   }
   str[ret]=0;
   return ;
}

/* ---------------------------------------------------------------------- */
int  read_ascii_line(FILE *fd , char *str, int maxch)
/*
   Legge da file dati riga trascurando :
      righe vuote  
      righe che iniziano con #
      la parte di riga a sinistra del carattere | 
      errore = status -1
*/

{
      char     *pch ,*p ;
      int      status ,nb=0 ;
      int      i,ii,trovata,n;
	  
      trovata = 1;
      while(trovata) {
         p=fgets(str, maxch, fd);
         if(p == NULL ) return -1;
         spazi(str);
         if(str[0] != '\n' && str[0] != '#' && (strlen(str)>0)){
	       if((pch = strchr(str,'\n')) != NULL) *pch='\0';
           if((pch = strchr(str,'|'))  != NULL) *pch='\0';
           i=strlen(str)-1;
           while(str[i] == 32 & i > 0)--i;
           str[i+1]='\0';
           trovata=0;
         }
      }
      return strlen(str);
}


//Simple selection sort
void selectionSort(int n, double *niz)
{
	double temp;
	for (int i = 0; i < n-1; i++)
	{
		int min = i;
		for (int j = i+1; j < n; j++)
		{
			if (niz[j] < niz[min]) min = j;
		}
		if (min != i)
		{
			temp = niz[i];
			niz[i] = niz[min];
			niz[min] = temp;
		}
	}
}

void copyStruct(struct DataInput source, struct DataInput *dest)
{
	dest->A0e = source.A0e;
	dest->A1e = source.A1e;
	dest->an01 = source.an01;
	dest->c = source.c;
	dest->dr = source.dr;
	dest->dt = source.dt;
	dest->i_type = source.i_type;
	dest->n0 = source.n0;
	dest->n0e = source.n0e;
	dest->n1 = source.n1;
	dest->n1e = source.n1e;
	dest->ne = source.ne;
	dest->nh = source.nh;
	dest->nhi = source.nhi;
	dest->ni = source.ni;
	dest->n_cw = source.n_cw;
	dest->n_kj = source.n_kj;
	dest->n_rec = source.n_rec;
	dest->n_tpsf = source.n_tpsf;
	dest->pi = source.pi;
	dest->precisione = source.precisione;
	dest->R = source.R;
	for (int i = 0; i < N_REC_MAX; i++)
	{
		dest->rec[i] = source.rec[i];
	}
	dest->rmin = source.rmin;
	dest->ro = source.ro;
	dest->tmin = source.tmin;
	dest->ua0 = source.ua0;
	dest->ua1 = source.ua1;
	dest->ud0 = source.ud0;
	dest->ud1 = source.ud1;
	dest->usR0 = source.usR0;
	dest->usR1 = source.usR1;
	dest->v0 = source.v0;
	dest->v1 = source.v1;
	dest->xacc = source.xacc;
	dest->z0 = source.z0;
	dest->z1 = source.z1;
}

