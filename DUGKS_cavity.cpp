//DUGKS  cavity flow  bounce back  boundary condition  in t+h time step
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <time.h>
const int M=128;
const int N=128;
const int Q=9;
const int M1=M+1;
const int N1=N+1;
const int M2=M+2;
const int N2=N+2;
const int IMIN=1;
const int JMIN=1;
const int IMAX=N;
const int JMAX=M;
const double PI=3.1415926;
const double Lx=1.0;
const double Ly=1.0;
const double U0=0.1;
const double Re=1000.0;
const double RT=1.0/3.0;
const double CFL=0.95;
double f[M2][N2][Q],fb[M2][N2][Q],w[Q],e[Q][2];//f=f~   fb=f-+  cell ceter
double fin[M1][N1][Q],fbinx[M2][N1][Q],fbiny[M1][N2][Q];//fin=f- interface   fbin=f-+ interface
double fx[M1][N1][Q],fy[M1][N1][Q];
double rho[M2][N2],u[M2][N2],v[M2][N2],ua[M2][N2],va[M2][N2];
double dx,dy,dt,h,tau,vis,CS,wx,wy,err;
void initial(void),flux_x(void),flux_y(void);
void macro(void),datadeal(void),boundary(void),evol(void);
double feq(int k, double RHO, double U, double V);
int main()
{
 int m,mmax,TEND,goon,M3,N3,j,i;
 double u0,d,duration,finish,start;
 double Erra,Errb;
 m=0;M3=M/2;N3=N/2;
 TEND=0;
 err=1.0;
 dx=Lx/N;
 dy=Ly/M;
 vis=Lx*U0/Re;
 CS=sqrt(RT);
 
 dt=CFL*dx/(sqrt(6*RT));

 tau=vis/RT;
 h=dt/2;
 wx=dt/dx;
 wy=dt/dy;
 
 	 d=sqrt(3*RT);
 	e[0][0]=0;  e[0][1]=0;//0-x,1-y
	e[1][0]=d;  e[1][1]=0;
	e[2][0]=0;  e[2][1]=d;
	e[3][0]=-d; e[3][1]=0;
	e[4][0]=0;  e[4][1]=-d;
	e[5][0]=d;  e[5][1]=d;
	e[6][0]=-d; e[6][1]=d;
	e[7][0]=-d; e[7][1]=-d;
	e[8][0]=d;  e[8][1]=-d;
	
		
	w[0]=4.0/9;
	w[1]=w[2]=w[3]=w[4]=1.0/9;
	w[5]=w[6]=w[7]=w[8]=1.0/36;
	
	initial();
	
	AA:
    printf("dt=%f,tau=%f\n",dt,tau);
	printf("input mmax:\n");
	scanf("%d",&mmax);
	TEND+=mmax ;
    start=clock();
	while((m<TEND)&&err>1.0e-6)//convergence?
    {
		m++;
		boundary();
	    flux_x();
	    flux_y();
	    evol();
		if(m%1000==0)
		{
		 Erra=0.0;Errb=0.0;
		 for(i=IMIN;i<=IMAX;i++)
           for(j=JMIN;j<=JMAX;j++)
		   {
		   Erra+=(fabs(u[j][i]-ua[j][i])+fabs(v[j][i]-va[j][i]));
		   }
		 for(i=IMIN;i<=IMAX;i++)
            for(j=JMIN;j<=JMAX;j++)
		   {
		   Errb+=fabs(u[j][i])+fabs(v[j][i]);
		   ua[j][i]=u[j][i];va[j][i]=v[j][i];//上个时间步的速度
		   }
		   err=Erra/Errb;//全局相对误差
		   printf("err=%e m=%d\n",err,m);
		 datadeal();
		}
		
    } 
	
	finish=clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf( "CPU time=%2.2f seconds\n", duration );
	printf("go on?\n");
	scanf("%d", &goon);
	if(goon) goto AA;
	return 0;
}
void initial()
{
  int j,i,k;
  
  for(i=IMIN;i<=IMAX;i++)
   for(j=JMIN;j<=JMAX;j++)
     {
	 rho[j][i]=1.0;
	 u[j][i]=0.0;
	 v[j][i]=0.0;
	 }
	 for(i=IMIN;i<=IMAX;i++)
      for(j=JMIN;j<=JMAX;j++)
	   for(k=0;k<=8;k++)
	      {
		  f[j][i][k]=feq(k,rho[j][i],u[j][i],v[j][i]);
		  }
}
double feq(int k,double RHO,double U,double V)
 {
  double re,eu,uv;
  eu=e[k][0]*U+e[k][1]*V;
  uv=U*U+V*V;
  re=w[k]*RHO*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
  return(re);
 }
 void boundary()
 {  int j,i,k;
    double w1,w2,RW;
     w1=(2*tau-h)/(2*tau+dt);
     w2=3*h/(2*tau+dt);
   for(j=JMIN;j<=JMAX;j++)
      for(i=IMIN;i<=IMAX;i++)
	  {
	   for(k=0;k<=8;k++)
	   {
	    fb[j][i][k]=w1*f[j][i][k]+w2*feq(k,rho[j][i],u[j][i],v[j][i]);
	   }
	  }
	 
	  for(i=IMIN;i<=IMAX;i++)for(k=0;k<=8;k++)
	{
	fb[JMIN-1][i][k]=2*fb[JMIN][i][k]-fb[JMIN+1][i][k];
	fb[JMAX+1][i][k]=2*fb[JMAX][i][k]-fb[JMAX-1][i][k];
	}
	
	
	for(j=JMIN-1;j<=JMAX+1;j++)for(k=0;k<=8;k++)
	   {
	   fb[j][IMIN-1][k]=2*fb[j][IMIN][k]-fb[j][IMIN+1][k];
	   fb[j][IMAX+1][k]=2*fb[j][IMAX][k]-fb[j][IMAX-1][k];
	   }
	  
	
  for(j=JMIN-1;j<=JMAX+1;j++)
      for(i=IMIN-1;i<=IMAX;i++)
	    for(k=0;k<=8;k++)
		{
		fbinx[j][i][k]=0.5*(fb[j][i][k]+fb[j][i+1][k]);
		}
	for(i=IMIN-1;i<=IMAX+1;i++)
      for(j=JMIN-1;j<=JMAX;j++)
	    for(k=0;k<=8;k++)
		{
		fbiny[j][i][k]=0.5*(fb[j][i][k]+fb[j+1][i][k]);
		}
 }
 
 
 void flux_x()
 {
 int j,i,k;
 double w1,w2;
 double grad;
 double rh,uh,vh;
 w1=(2*tau)/(2*tau+h);
 w2=h/(2*tau+h);
     for(j=JMIN;j<=JMAX;j++)
       for(i=IMIN-1;i<=IMAX;i++)
	    for(k=0;k<=8;k++)
	   {  
	      grad=e[k][0]*(fb[j][i+1][k]-fb[j][i][k])/dx+e[k][1]*(fbinx[j+1][i][k]-fbinx[j-1][i][k])/dy/2;
	      fin[j][i][k]=fbinx[j][i][k]-h*grad;
	   }
	   //bounce back
	   i=IMIN-1;
	for(j=JMIN;j<=JMAX;j++)
	{
	fin[j][IMIN-1][1]=fin[j][IMIN-1][3];
	fin[j][IMIN-1][5]=fin[j][IMIN-1][7];
	fin[j][IMIN-1][8]=fin[j][IMIN-1][6];
	
	fin[j][IMAX][3]=fin[j][IMAX][1];
	fin[j][IMAX][6]=fin[j][IMAX][8];
	fin[j][IMAX][7]=fin[j][IMAX][5];
	
	}
	 for(j=JMIN;j<=JMAX;j++)
      for(i=IMIN-1;i<=IMAX;i++)
      {   
	      rh=0.0;uh=0.0;vh=0.0;
		  for(k=0;k<=8;k++)
             {
			 rh+=fin[j][i][k];
			 uh+=fin[j][i][k]*e[k][0];
			 vh+=fin[j][i][k]*e[k][1];
			 }
	   uh=uh/rh;
	   vh=vh/rh;
	   for(k=0;k<=8;k++)
	      fx[j][i][k]=w1*fin[j][i][k]+w2*feq(k,rh,uh,vh);
	  }
	  
 }

void flux_y()
 {
 int j,i,k;
 double w1,w2;
 double grad;
 double rh,uh,vh,RW;
 w1=(2*tau)/(2*tau+h);
 w2=h/(2*tau+h);
   
      for(i=IMIN;i<=IMAX;i++)for(j=JMIN-1;j<=JMAX;j++)
	    for(k=0;k<=8;k++)
	   {  
	      grad=e[k][0]*(fbiny[j][i+1][k]-fbiny[j][i-1][k])/dx/2+e[k][1]*(fb[j+1][i][k]-fb[j][i][k])/dy;
	      fin[j][i][k]=fbiny[j][i][k]-h*grad;
	   }
	 //下板
	for(i=IMIN;i<=IMAX;i++)
	{
	RW=rho[j-1][i];
	//RW=0.5*(fin[j][i][0]+fin[j-1][i][0]+fin[j][i][1]+fin[j-1][i][1]+fin[j][i][3]+fin[j-1][i][3]);
   // RW=RW+fin[j][i][2]+fin[j-1][i][2]+fin[j][i][5]+fin[j-1][i][5]+fin[j][i][6]+fin[j-1][i][6];
	//RW=RW/(1-2*(w[5]*e[5][0]*U0+w[6]*e[6][0]*U0)/RT);
	fin[JMIN-1][i][2]=fin[JMIN-1][i][4];
    fin[JMIN-1][i][5]=fin[JMIN-1][i][7];
	fin[JMIN-1][i][6]=fin[JMIN-1][i][8];
	
	fin[JMAX][i][4]=fin[JMAX][i][2];
	fin[JMAX][i][7]=fin[JMAX][i][5]+2*w[7]*RW*U0*e[7][0]/RT;
	fin[JMAX][i][8]=fin[JMAX][i][6]+2*w[8]*RW*U0*e[8][0]/RT;
	}
	
      for(i=IMIN;i<=IMAX;i++)for(j=JMIN-1;j<=JMAX;j++)
      {
	    rh=0.0;uh=0.0;vh=0.0;
		  for(k=0;k<=8;k++)
             {
			 rh+=fin[j][i][k];
			 uh+=fin[j][i][k]*e[k][0];
			 vh+=fin[j][i][k]*e[k][1];
			 }
	   uh=uh/rh;
	   vh=vh/rh;
	   for(k=0;k<=8;k++)
	      fy[j][i][k]=w1*fin[j][i][k]+w2*feq(k,rh,uh,vh);
	  } 
 
 } 
 void evol()
 {
  int j,i,k;
  double fa;
    for(i=IMIN;i<=IMAX;i++)
      for(j=JMIN;j<=JMAX;j++)
	   for(k=0;k<=8;k++)
         {
		 fa=fb[j][i][k]*4/3-f[j][i][k]/3;
		 f[j][i][k]=fa-wx*e[k][0]*(fx[j][i][k]-fx[j][i-1][k])-wy*e[k][1]*(fy[j][i][k]-fy[j-1][i][k]);		 
		 }
	for(i=IMIN;i<=IMAX;i++)
      for(j=JMIN;j<=JMAX;j++)
		   {
		   rho[j][i]=0.0;u[j][i]=0.0;v[j][i]=0.0;
		   for(k=0;k<=8;k++)
             {
			 rho[j][i]+=f[j][i][k];
			 u[j][i]+=f[j][i][k]*e[k][0];
			 v[j][i]+=f[j][i][k]*e[k][1];
			 }
	   
	      u[j][i]=u[j][i]/rho[j][i];
	      v[j][i]=v[j][i]/rho[j][i];
		 }
 }

void datadeal()
 {
 int j,i;
 FILE *fp;
 if((fp=fopen("error.dat","a"))==NULL)
	{
		printf(" open file error\n");
		exit(1);
	}
	   fprintf(fp,"%e ",err);  
	fprintf(fp,"\n");
	fclose(fp);
 if((fp=fopen("rho.dat","w"))==NULL)
	{
		printf(" open file error\n");
		exit(1);
	}
	for(j=JMIN;j<=JMAX;j++)
      {
	  for(i=IMIN;i<=IMAX;i++)
	   fprintf(fp,"%e ",rho[j][i]*RT);  
	fprintf(fp,"\n");
	  }
	
	fclose(fp);
	
	 if((fp=fopen("u.dat","w"))==NULL)
	{
		printf(" open file error\n");
		exit(1);
	}
	for(j=JMIN;j<=JMAX;j++)
	
      {for(i=IMIN;i<=IMAX;i++)
	  		fprintf(fp,"%e ",u[j][i]/U0);  
	fprintf(fp,"\n");
	  }
	
	fclose(fp);
	
	 if((fp=fopen("v.dat","w"))==NULL)
	{printf(" open file error\n");
		exit(1);
	}
	for(j=JMIN;j<=JMAX;j++)
	
      {for(i=IMIN;i<=IMAX;i++)
		fprintf(fp,"%e ",v[j][i]/U0);  
	fprintf(fp,"\n");
	  }
	
	fclose(fp);
	
 }
