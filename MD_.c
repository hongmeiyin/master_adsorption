# include <unistd.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "mpi.h"
 
 
# define N     30                 /* Number of Residues(amino acids)                   */
# define Nwin  40
# define NC    250
double T;
double bn=4;  // bond length
 
# define MDCYC 80000000           /* Number of running cycles of MD                    */
# define MDCON 50                 /* Number of steps per cycle                         */
 
# define max(A,B) ((A)>(B)?(A):(B)) 
# define min(A,B) ((A)<(B)?(A):(B))
 
# define NE   1000
# define NBIN 1000               /* Number of  bins                                   */
# define NOBS 7                 /* Number of observable qualities                    */
# define IRES 100000            /* the inverval for saving in formation for restart  */
# define NTHERM 10000000           /* Number of discarded cycles                      */
# define IRT 100000               /* Number of cycles for output energy et. al       */
# define RT "LD_RT.out"          /* Enerfy, Rg, RMSD, Q, general qualities            */
# define NMAX    5000
# define MAXREP  5000
 
# define zdn  0
# define zup  240.0
 
////////////desolvation//////////////
int f_desolv=0;
double ksolv=0.1;                  
double ksolv1=0.2; 
////////////desolvation//////////////
 
long int con_t=0;
 
int ProgID, ProcNum;
void CloseMPI(void);
int NumRep;
void RunREM(void);
double xxArray[Nwin][N], yyArray[Nwin][N], zzArray[Nwin][N];
double vxArray[Nwin][N], vyArray[Nwin][N], vzArray[Nwin][N];
double fxArray[Nwin][N], fyArray[Nwin][N], fzArray[Nwin][N];
double EArray[Nwin];
double NumofExchange[Nwin], NumofExAccept[Nwin];
double delta_B;
double Tm[Nwin];
 
/************* geometry ***************************************************************/
double x[N], y[N], z[N];           /* atom coordinates                                */
double b[N];                       /* pseudo bond lenghts (0...N-1)                   */
double th[N];                      /* pseudo bond angles (1...N-2)                    */
double bx[N],by[N],bz[N];          /* pseudo bond vectors                             */
/************* energies and forces and velocities *************************************/
double Ekin,Epot,Eint,Eben,Ebon,Erep,Econ,Evol, Egft=0;   /* energy terms                                    */
double fx[N],fy[N],fz[N];          /* conformational force                            */
double fxo[N],fyo[N],fzo[N];       /* conformational force old                        */
double frdx[N],frdy[N],frdz[N];    /* random force                                    */
double frdxo[N],frdyo[N],frdzo[N]; /* random force old                                */
double vx[N],vy[N],vz[N];          /* velocities                                      */
/************* MD parameters **********************************************************/
const double tau=4.0;              /* time scale sqrt(m*l*l/eps0)                     */
double gam;                        /* friction coefficient                            */
double dt;                         /* time step                                       */
long imd;
/************* interactional parameters ***********************************************/
double epsilon=1.0;                /* energy unit                                     */
double kbon=100.0;                 /* coefficient of bond stretching                  */
double kth=8.0;                      /* coefficient of bond bending                     */
double krep=1;                     /* coefficient of unbonded non-native term         */
double kcont=1.0;
double kwp=0;
double kvol=100;
double sigsa=4.0;                  /* atom radius (-hard core potential in AA)        */
double cut=9.0;                    /* cutoff self avoidance                           */
/************* native structure *******************************************************/
/************* miscellaneous **********************************************************/
double c1,c2,c3;                   /* coefficients for getting new velocities         */
double pi,pi2,pid2;                /* pai,2*pai,pai/2                                 */
double tconst;                     /* random force magnitude                          */
double kcont60,krep12,kbon2;       /* 60*kcont, 12*krep, 2*kbon                       */
double cthmin,sthmin;
int carterr=0,therr=0;             /* for unreasonable conformation                   */
int flag=-1;                       /* judge teh conformation reasonability (1 or 0)   */
 
//long int iseed=-11;                 
//double gcos;                       /* for generating random number                    */
//long int ma[56];         /* for generating random number                    */
//int inext,inextp,IFF,ISET;         /* for generating random number                    */
 
long int ir[98];
long int iy;
int jjj;
long int iseed=-11, iseed1=-13;
 
int phase=0;
double V1, V2;
/************* functions **************************************************************/
int double_to_int(double x);        
void init(void);
void printhead(); 
void mdstep();
double bond();                  
double bend();                  
double rep(); // the fuction same as sac(), for the E_rep term
double gasdev2();               
double irand(long int *idum);
void histoe(int iflag,double x);
void historg(int iflag,double x);
int cart2dof(); 
double gyr2(); 
void runtime(char *fn,long it,double o[],int n);
double pore_rep();
double cont(int iflag);
double graft();
 
 
double pbin[NE];
double Elow=-200, Ehigh=300;
double dbin;
 
int ico=0;
long int step0=0;
 
int ip1[MAXREP], ip2[MAXREP];
int npair;
 
double xc0, yc0, zc0;
 
double Rgz2, Rgxy2;
///////////////////////////////////////////////////////////////////
int main (int argc,char *argv[])  // main function
{
  int i, j, m, k, nnc;
  double o[NOBS],so[NOBS],po[NOBS];
  int num_c=0,model=0;
  double rg, abc;
  char fn[100], Istring[60], fpname[60], filename[60], coorname[60];
  FILE *fp, *fp1, *fp2, *fp3, *fp4;
 
  double opt[10000][7];
  double E1=0, E2=0, istep=0, Cv;
  double U1=0, U2=0;
 
  double A1=0, A2=0;
 
  double rx, ry, rz, r2, rij;
  int histq[NMAX], nt=0;
  double rg1=5, rg2=205, drg=0.2, rgsq;
  double hisrg[NBIN]; 
  int rge[NBIN][NBIN];
  double pben[NE];
  double pcon[NE]; 
  double irg1=0, irg2=0, crg;
  double zbin=0.4, xm, ym, zm;  // w-p distance
  double Emin=-900, Emax=100, dE;
  double Eblow=0, Ebhigh=400,dF=0.4;
  double Eclow=-300, Echigh=100,dG=0.4;
  double dcm;
  double pscE=0, psc=0, Eave=0;
  long int icount=0;
  double t1, t2, t3;
 
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProgID);
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
 
  NumRep=Nwin;
  MPI_Bcast(&NumRep,1,MPI_INT,0,MPI_COMM_WORLD);
 
 
  for(i=0;i<NBIN;i++)
  {
     hisrg[i]=0;
  }
 
 
  for(i=0;i<NBIN;i++)
  {
      for(j=0;j<NBIN;j++)
      {
          rge[i][j]=0;
      }
  }
 
 
  fp=fopen("TList","r");
  for(i=0;i<Nwin;i++)
  {
     fscanf(fp,"%lf",&Tm[i]);
  }
  fclose(fp);
 
  MPI_Bcast(Tm,Nwin,MPI_DOUBLE,0,MPI_COMM_WORLD);
  T=Tm[ProgID];
 
  srand((unsigned)getpid());  // generate random seed
  abc=rand()/(double)RAND_MAX;
  iseed=-10000*abc;
 
  dbin=(Ehigh-Elow)/NBIN;
 
  for(i=0;i<NE;i++)
  {
     pbin[i]=0;
 
     pben[i]=0;
     pcon[i]=0;
  }
 
  for(i=0;i<NMAX;i++)
  {
     histq[i]=0;
  }
 
  sprintf(Istring,"LD_RT%d.out",ProgID);
 
  for(i=0;i<NOBS;i++) o[i]=po[i]=so[i]=0;
 
  init();     //Initiation 
 
  printf("Initial conformation\n");
  printf("T: %f\tepsilon: %f\te/T: %f\n", T, epsilon, epsilon/T);
  printf("Ekin %f  Epot %f\n",Ekin,Epot);
  printf("Ebon %g  Eben %g  Econ  %g\n",Ebon,Eben,Econ);
  printhead();
 
  for(m=0;m<NTHERM;m++)
  {
      for(j=0;j<MDCON;j++) mdstep();
 
      if((m+1) % 2000 == 0)
      {
          RunREM();
          for(i=0;i<N;i++)
          {
             fx[i]=0; fy[i]=0; fz[i]=0;
          }
          Ebon=bond(); Eben=bend(); Econ=cont(0); 
          Epot=Ebon+Eben+Econ;
      }
  }
  for(imd=con_t/MDCON;imd<MDCYC;imd++) 
  {
      for(j=0;j<MDCON;j++)
      {
          mdstep();
          E1+=Epot;
          E2+=Epot*Epot;
          istep+=1;
          i=(int)((Epot-Elow)/dbin);
          if(i>=0 && i<NE)
          {
             pbin[i]+=1;
          }
          
      }
 
      icount+=1;
 
      i=(int)((Eben-Eblow)/dF);
      if(i>=0 && i<NE)
      {
          pben[i]+=1;
      }
 
      i=(int)((Econ-Eclow)/dG);
      if(i>=0 && i<NE)
      {
         pcon[i]+=1;
      }
 
      rgsq=gyr2();
      rg=sqrt(rgsq);
 
      irg1+=rg;
      irg2+=(rg*Epot);
 
      Eave+=Epot;
       
 
      nnc=0;
      for(m=0;m<npair;m++)
      {
          i=ip1[m]; j=ip2[m];
          rx=x[j]-x[i]; 
          ry=y[j]-y[i]; 
          rz=z[j]-z[i];
          r2=rx*rx+ry*ry+rz*rz;
          rij=sqrt(r2);
          if(rij < 1.2*sigsa)
          {
             nnc+=1;
          }
      }
 
      histq[nnc]+=1;
 
      k=(rg-rg1)/drg;
      if(k >=0 && k<NBIN)
      {
          hisrg[k]+=1;          
      }
 
      i=(Epot-Elow)/dbin;
      rge[i][k]+=1;
 
 
  
      o[0]=Ekin; o[1]=Epot; o[2]=Ebon; o[3]=Eben; o[4]=Econ; o[5]=rg; o[6]=nnc;    
 
 
      if((imd+1)%IRT==0) runtime(Istring,(imd+1),o,NOBS);
 
      if((imd+1)% 200000==0)
      {
          sprintf(fpname,"conf%d",ProgID);
          fp2=fopen(fpname,"a");
          for(i=0;i<N;i++)
          {
              fprintf(fp2,"%g  %g  %g\n",x[i],y[i],z[i]);
          }
          fprintf(fp2,"%g  %g  %g  %g  %g  %d\n",Epot,Ebon,Eben,Econ,rg,nnc);
          fprintf(fp2,"\n");
          fclose(fp2);
      }
 
      if((imd+1) % 2000 == 0)
      {
          RunREM();
          for(i=0;i<N;i++)
          {
             fx[i]=0; fy[i]=0; fz[i]=0;
          }
          Ebon=bond(); Eben=bend(); Econ=cont(0);
          Epot=Ebon+Eben+Econ;
      }
 
      if(imd>NTHERM && (imd+1)%2000000==0)
      {
          if(ProgID==0)
          {
             fp1=fopen("Ratio.dat","a");
             fprintf(fp1,"Step = %ld\n",imd);
             for(i=0; i<Nwin-1; i++)
             {
            fprintf(fp1,"Accepted ratio: [%8.5lf, %8.5lf]\t%7.3lf\n",Tm[i], Tm[i+1], (NumofExAccept[i]/NumofExchange[i]));
         }
             fclose(fp1);
          }        
      }
 
      if((imd+1)% 1000000==0)
      {
          sprintf(fpname,"PE%d",ProgID);
          fp2=fopen(fpname,"w");
          for(i=0;i<NE;i++)
          {
             fprintf(fp2,"%d  %g\n",i,pbin[i]);
          }
          fclose(fp2); 
 
          sprintf(fpname,"Ebend%d",ProgID);
          fp2=fopen(fpname,"w");
          for(i=0;i<NE;i++)
          {
             fprintf(fp2,"%g  %g\n",Eblow+(i+0.5)*dF,pben[i]);
          }
          fclose(fp2);
 
          sprintf(fpname,"Econ%d",ProgID);
          fp2=fopen(fpname,"w");
          for(i=0;i<NE;i++)
          {
             fprintf(fp2,"%g  %g\n",Eclow+(i+0.5)*dG,pcon[i]);
          }
          fclose(fp2);
  
           
          Cv=(E2/istep)-(E1/istep)*(E1/istep);
          Cv=Cv/T/T;
 
          sprintf(filename,"CV%d",ProgID);
          fp3=fopen(filename,"a"); 
          fprintf(fp3,"%g  %g\n",T,Cv);   
          fclose(fp3);
 
          crg = (irg2/icount)-(Eave/icount)*(irg1/icount);   
          crg=crg/T/T;
 
          sprintf(filename,"CRG%d",ProgID);
          fp3=fopen(filename,"a"); 
          fprintf(fp3,"%g  %g\n",T,crg);   
          fclose(fp3);
           
 
          nt=0;
          for(i=0;i<=npair;i++)
          {
             nt+=histq[i];
          }   
 
          sprintf(fpname,"FE%d",ProgID);
          fp2=fopen(fpname,"w");
          for(i=0;i<NC;i++)
          {
             fprintf(fp2,"%d  %g\n",i,-log(1.0*histq[i]/nt));
          }
          fclose(fp2); 
 
          sprintf(fpname,"RG%d",ProgID);
          fp2=fopen(fpname,"w");
          for(i=0;i<NBIN;i++)
          {
             fprintf(fp2,"%g  %g\n",rg1+(i+0.5)*drg,hisrg[i]);
          }
          fclose(fp2);  
 
/*
          sprintf(fpname,"Ep_RG%d",ProgID);
          fp2=fopen(fpname,"w");
          for(i=0;i<NBIN;i++)
          {
             for(j=0;j<600;j++)
             {
                fprintf(fp2,"%d  %d  %d\n",i, j, rge[i][j]);
             }
          }
          fclose(fp2);*/
      }
     
  }
     
  CloseMPI();
  return 0;
}
/****************************************************************************/
/***** MEASUREMENTS *********************************************************/
/****************************************************************************/
 
void RunREM(void)
{
    int i, j;
    double Etmp;
    double P, h;
    double xp[N], yp[N], zp[N];
    double vxsum, vysum, vzsum;
    double Vmi, Vmj, Vni, Vnj;
    if(ProcNum==1)
    {
       return;
    }
    MPI_Gather(x, N, MPI_DOUBLE, xxArray[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(y, N, MPI_DOUBLE, yyArray[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(z, N, MPI_DOUBLE, zzArray[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Epot, 1, MPI_DOUBLE, EArray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    MPI_Gather(vx, N, MPI_DOUBLE, vxArray[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(vy, N, MPI_DOUBLE, vyArray[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(vz, N, MPI_DOUBLE, vzArray[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    MPI_Barrier(MPI_COMM_WORLD);
 
    if(ProgID==0)
    {
        for(i=Nwin-1;i>0;i--)
        {
            delta_B=(1.0/Tm[i]-1.0/Tm[i-1])*(EArray[i-1]-EArray[i]);        
            NumofExchange[i-1]+=1;
 
            if(delta_B>0)
            {
                P=exp(-delta_B);
                h=irand(&iseed1);
                if(h>=P)
                {
                   continue;
                }
            }
            /////////////////////////////////////////////////
            for(j=0;j<N;j++)
            {
                xp[j]=xxArray[i][j];
                yp[j]=yyArray[i][j];
                zp[j]=zzArray[i][j];                
            }
            for(j=0;j<N;j++)
            {
                xxArray[i][j]=xxArray[i-1][j];
                yyArray[i][j]=yyArray[i-1][j];
                zzArray[i][j]=zzArray[i-1][j];
            }
            for(j=0;j<N;j++)
            {
                xxArray[i-1][j]=xp[j];
                yyArray[i-1][j]=yp[j];
                zzArray[i-1][j]=zp[j];
            }
            Etmp=EArray[i]; EArray[i]=EArray[i-1]; EArray[i-1]=Etmp;
             
            /////////////////////////////////////////////////
 
            vxsum=vysum=vzsum=0;  //velocity
            for(j=0;j<N;j++) 
            {
                vxArray[i][j]=sqrt(Tm[i])*gasdev2();
                vyArray[i][j]=sqrt(Tm[i])*gasdev2();
                vzArray[i][j]=sqrt(Tm[i])*gasdev2();
 
                vxsum+=vxArray[i][j]; 
                vysum+=vyArray[i][j];
                vzsum+=vzArray[i][j];
            }
            for(j=0;j<N;j++) 
            {
                vxArray[i][j]-=vxsum/N;
                vyArray[i][j]-=vysum/N;
                vzArray[i][j]-=vzsum/N;
            }
 
            vxsum=vysum=vzsum=0;
            for(j=0;j<N;j++) 
            {
                vxArray[i-1][j]=sqrt(Tm[i-1])*gasdev2();
                vyArray[i-1][j]=sqrt(Tm[i-1])*gasdev2();
                vzArray[i-1][j]=sqrt(Tm[i-1])*gasdev2();
 
                vxsum+=vxArray[i-1][j]; 
                vysum+=vyArray[i-1][j];
                vzsum+=vzArray[i-1][j];
            }
            for(j=0;j<N;j++) 
            {
                vxArray[i-1][j]-=vxsum/N;
                vyArray[i-1][j]-=vysum/N;
                vzArray[i-1][j]-=vzsum/N;
            }
            NumofExAccept[i-1]+=1.0;            
        }
    }
 
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(xxArray[0], N, MPI_DOUBLE, x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(yyArray[0], N, MPI_DOUBLE, y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(zzArray[0], N, MPI_DOUBLE, z, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    MPI_Scatter(vxArray[0], N, MPI_DOUBLE, vx, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(vyArray[0], N, MPI_DOUBLE, vy, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(vzArray[0], N, MPI_DOUBLE, vz, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    MPI_Scatter(EArray, 1, MPI_DOUBLE, &Epot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
}
 
void CloseMPI(void)
{
     MPI_Finalize();
     return;
}
int double_to_int(double x)
{
     int bb,xx;
     bb=(int)x;
     if ((x-bb)*(x-bb)>=0.25) {xx=bb+1;}
     else xx=bb;
     return xx;
}
/****************************************************************************/
 
void histoe(int iflag,double x)
{
  static double low=0.0,high=600.0;
  static double his[NBIN],eps;
  static int out=0;
  double tmp1, tmp2;
  int j,k;
  FILE *fp;
 
  if(iflag<0) 
  {
    for(j=0;j<NBIN;j++)
    {
       his[j]=0;
    }
    eps=(high-low)/NBIN;
    printf("histoe: low %f high %f\n",low,high);
  }
   
  if(iflag==0) 
  {
    if(x>low && x<high) 
    {
      k=(x-low)/eps;
      his[k]++;
    }
    else
    {
      out++;
    }
 
  }
 
  if(iflag>0) 
  {
     fp=fopen("LD_hise","w");
     j=0; 
     for(k=0;k<NBIN;k++) j+=his[k];
     for(k=0;k<NBIN;k++)
     {
       if(j!=0)
       {
         fprintf(fp,"%f %g\n",low+(k+.5)*eps,his[k]/j/eps);
       }
     }
     fclose(fp);
  }
}
 
void historg(int iflag,double x) 
{
  static double low=5.0,high=205.0;  //  5.0,  40.0
  static double his[NBIN],eps;
  static int out=0;
  double tmp1, tmp2;
  int j, k;
  FILE *fp;
 
  if(iflag<0) 
  {
     for(j=0;j<NBIN;j++) his[j]=0;
     eps=(high-low)/NBIN;    
  }
   
  if(iflag==0) 
  {
    if(x>low && x<high) 
    {
      k=(x-low)/eps;
      his[k]++;
    }
    else out++;
     
  }
 
  if(iflag>0) 
  {
    fp=fopen("LD_hisrg","w");
    j=0; 
    for (k=0;k<NBIN;k++) j+=his[k];
    for (k=0;k<NBIN;k++)
      if (j!=0) fprintf(fp,"%f  %g\n",low+(k+.5)*eps,his[k]/j/eps);
    fclose(fp);    
  }
 
}
 
double bond() 
{
  int i,j;
  double fbx,fby,fbz,fb;
  double e=0,db;
 
  for(i=0;i<N-1;i++) 
  {
    j=i+1;
    db=b[i]-bn;
    e+=db*db;
    fb=-kbon2*db;
 
    fbx=fb*bx[i];   //revised by CHEN on Dec. 8 2008
    fby=fb*by[i];
    fbz=fb*bz[i];
 
    fx[i]-=fbx;
    fy[i]-=fby;
    fz[i]-=fbz;
 
    fx[j]+=fbx;
    fy[j]+=fby;
    fz[j]+=fbz;
  }
  return kbon*e;
}
/****************************************************************************/
double bend() 
{
  int i,j,k;
  double b1x,b1y,b1z,b1;
  double b2x,b2y,b2z,b2;
  double dix,diy,diz;
  double dkx,dky,dkz;
  double cth,sth,dth;
  double e=0,fben;
 
  for(i=0;i<N-2;i++) 
  {
    j=i+1;
    k=i+2;
    dth=th[j]-pi;  // 180 degree
    cth=max(cos(th[j]),-cthmin);
    sth=max(sin(th[j]),sthmin);
    if(sin(th[j])<sthmin) therr++;
     
    b1x=bx[i];
    b1y=by[i];
    b1z=bz[i];
    b1=b[i];
 
    b2x=bx[j];
    b2y=by[j];
    b2z=bz[j];
    b2=b[j];
 
    e+=dth*dth;
    fben=-2*kth*dth;
 
    dix=-(b2x+cth*b1x)/sth/b1;
    diy=-(b2y+cth*b1y)/sth/b1;
    diz=-(b2z+cth*b1z)/sth/b1;
 
    dkx=(b1x+cth*b2x)/sth/b2;
    dky=(b1y+cth*b2y)/sth/b2;
    dkz=(b1z+cth*b2z)/sth/b2;
 
    fx[i]+=fben*dix;
    fy[i]+=fben*diy;
    fz[i]+=fben*diz;
 
    fx[j]+=fben*(-dix-dkx);
    fy[j]+=fben*(-diy-dky);
    fz[j]+=fben*(-diz-dkz);
 
    fx[k]+=fben*dkx;
    fy[k]+=fben*dky;
    fz[k]+=fben*dkz;
  }
 
  return kth*e;
}
 
double rep()   //no-native interaction
{
  int i,j;
  double r2,r6,r12,rx,ry,rz,fr,e=0;
  double cut2,sigsa2;
  //double fs=0.8;
 
  cut2=cut*cut;
  sigsa2=sigsa*sigsa;   //sigsa is r_rep
 
/*
  for(i=0;i<N-2;i++)
  {
     for(j=i+2;j<N;j++)
     {
         rx=x[j]-x[i];
         ry=y[j]-y[i];
         rz=z[j]-z[i];
 
         if((r2=rx*rx+ry*ry+rz*rz)>cut2) continue;
         r6=sigsa2/r2; r6=r6*r6*r6; r12=r6*r6;
         e+=r12;
         fr=krep12*r12/r2;
 
         fx[i]-=fr*rx;
         fy[i]-=fr*ry;
         fz[i]-=fr*rz;
 
         fx[j]+=fr*rx;
         fy[j]+=fr*ry;
         fz[j]+=fr*rz;
     }
  }*/
 
  for(i=0;i<N;i++)
  {
      rz=z[i]-zup;
      if((r2=rz*rz)>cut2) continue;
      r6=sigsa2/r2; r6=r6*r6*r6; r12=r6*r6;
      e+=r12;
      fr=krep12*r12/r2;
 
      fz[i]+=fr*rz;
  }
 
  for(i=0;i<N;i++)
  {
      rz=z[i]-zdn;
      if((r2=rz*rz)>cut2) continue;
      r6=sigsa2/r2; r6=r6*r6*r6; r12=r6*r6;
      e+=r12;
      fr=krep12*r12/r2;
 
      fz[i]+=fr*rz;
  }
 
 
  return krep*e;
}
 
double cont(int iflag)   //potential of native contacts
{
  int i,j,m,n,k, tmp=0, md=0;
  double r2,r4,r6,rx,ry,rz;
  double sig,sig2,fr,e=0; 
 
  if(iflag<0)
  { 
     for(i=0;i<N;i++)
     {
        for(j=i+2;j<N;j++)
        {
             ip1[tmp]=i;  ip2[tmp]=j;
             tmp++;
        }
     }
    
     npair=tmp;
     printf("npair: %d\n",tmp);
     return 0;
  }
 
  for(k=0;k<npair;k++)
  {
      i=ip1[k];  j=ip2[k];
      rx=x[j]-x[i]; 
      ry=y[j]-y[i]; 
      rz=z[j]-z[i];
 
      r2=rx*rx+ry*ry+rz*rz;
      sig2=sigsa*sigsa;
      if((r2>9*sig2)||(r2<1e-8)) continue;   // cutoff
      r6=sig2/r2; r4=r6*r6; r6=r4*r6;
      e+=kcont*r6*(r6-2.0);
      fr=12.0*kcont*r6*(r6-1.0)/r2;
 
      fx[i]-=fr*rx;
      fy[i]-=fr*ry;
      fz[i]-=fr*rz;
 
      fx[j]+=fr*rx;
      fy[j]+=fr*ry;
      fz[j]+=fr*rz;  
  } 
  
  return e;
}
 
 
double graft()
{
    int i,j;
    double rx, ry, rz, dr, drc, fr;
    double e=0,db;
 
    rx=x[0]-0;
    ry=y[0]-0;
    rz=z[0]-0;
 
    drc=sqrt(rx*rx+ry*ry+rz*rz);
    dr=drc-sigsa;
    e+=kbon*dr*dr;
 
    fr=-2*kbon*dr;
 
    fx[0]+=fr*rx/drc;
    fy[0]+=fr*ry/drc;
    fz[0]+=fr*rz/drc;
 
    return e;
}
 
double gyr2()     //mean-square radius of gyration
{
  int i;
  double rmx,rmy,rmz,rg2;
 
  Rgz2=0; Rgxy2=0;
 
  rmx=rmy=rmz=rg2=0;
  for(i=0;i<N;i++) 
  {
    rmx+=x[i];
    rmy+=y[i];
    rmz+=z[i];
  }
  rmx/=N; rmy/=N; rmz/=N;
  for(i=0;i<N;i++) 
  {
    rg2+=(x[i]-rmx)*(x[i]-rmx)+
         (y[i]-rmy)*(y[i]-rmy)+
         (z[i]-rmz)*(z[i]-rmz);
 
    Rgxy2+=(x[i]-rmx)*(x[i]-rmx)+
           (y[i]-rmy)*(y[i]-rmy);
    Rgz2+=(z[i]-rmz)*(z[i]-rmz);
  }
  return rg2/N;
}
 
 
int cart2dof()   //Caculation of bond, bond angles, torsion angles.
{
  int i,j,ok=1;
  double b1x,b1y,b1z,b1;
  double b2x,b2y,b2z;
  double ux,uy,uz,u;
  double tmp1;
 
  /* bond vectors */
  for(i=0;i<N-1;i++) 
  {
    j=i+1;
    b1x=x[j]-x[i];
    b1y=y[j]-y[i];
    b1z=z[j]-z[i];
    b1=sqrt(b1x*b1x+b1y*b1y+b1z*b1z);
     
    bx[i]=b1x/b1;
    by[i]=b1y/b1;
    bz[i]=b1z/b1;
    b[i]=b1;
    if(b1>5.0 || b1<2.0)
    {
       ok=0; 
       //printf("cart2dof 1: %i b1 %f\n",i,b1);
    }
  }    
 
  /* bond angles */
  for(i=1;i<N-1;i++) 
  {
    b1x=bx[i-1];
    b1y=by[i-1];
    b1z=bz[i-1];
 
    b2x=bx[i];
    b2y=by[i];
    b2z=bz[i];
     
    tmp1=b1x*b2x+b1y*b2y+b1z*b2z;
    if (tmp1<=-1.0) {tmp1=-1+(1e-8); }
    if (tmp1>=1) {tmp1=1-(1e-8);}
    if (acos(tmp1)==0.0) {ok=0;}// printf("cart2dof 4: tmp1 %f\n",tmp1);}
    th[i]=pi-acos(tmp1);
  }
  return ok;
}
 
void mdstep() 
{
  int i;
 
  for(i=0;i<N;i++)
  {
     x[i]=x[i]+dt*vx[i]+c1*(fx[i]-gam*vx[i]+frdx[i]);
     y[i]=y[i]+dt*vy[i]+c1*(fy[i]-gam*vy[i]+frdy[i]);
     z[i]=z[i]+dt*vz[i]+c1*(fz[i]-gam*vz[i]+frdz[i]);
  }
  
  
  if(1!=cart2dof()) {carterr++; flag=0;} else flag=1;
 
  for(i=0;i<N;i++) 
  {
    frdxo[i]=frdx[i];
    frdyo[i]=frdy[i];
    frdzo[i]=frdz[i];
 
    fxo[i]=fx[i];
    fyo[i]=fy[i];
    fzo[i]=fz[i];
 
    frdx[i]=gasdev2()*tconst;
    frdy[i]=gasdev2()*tconst;
    frdz[i]=gasdev2()*tconst;
  }
 
  for(i=0;i<N;i++) fx[i]=fy[i]=fz[i]=0;
  Epot=(Ebon=bond())+(Eben=bend())+(Econ=cont(0));
 
  Ekin=0;
  for(i=0;i<N;i++) 
  {
     vx[i]=c2*vx[i]+c3*(fxo[i]+fx[i]+frdxo[i]+frdx[i]);
     vy[i]=c2*vy[i]+c3*(fyo[i]+fy[i]+frdyo[i]+frdy[i]);
     vz[i]=c2*vz[i]+c3*(fzo[i]+fz[i]+frdzo[i]+frdz[i]);
     Ekin+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  Ekin*=0.5;
}
 
void init(void) 
{
  int i,j,m;
  double xn[N],yn[N],zn[N];
  double vxsum,vysum,vzsum;
  double r2,r2min,r2max,c0;
  double x0, y0, z0;
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double E;
  FILE *fp;
  FILE *fp1;
 
  pi=acos(-1.);
  pi2=2.*pi;
  pid2=pi/2.;
 
  cthmin=cos(pi/90.);
  sthmin=sin(pi/90.);
 
  /* energy parameters */
  kbon*=epsilon;
  kth*=epsilon;
  krep*=epsilon;
  krep12=krep*12;
  kbon2=kbon*2;
   
  /* MD parameters */
  dt=0.005*tau;
  gam=0.05/tau;
  tconst=sqrt(2*gam*T/dt);
  c0=gam*dt/2.;
  c1=dt*dt/2.;
  c2=(1-c0)*(1-c0+c0*c0);
  c3=dt*(1-c0+c0*c0)/2.;
 
  printf("seed %ld\n",iseed);
///////////////////////////////////////////////////  
/*
  x[0]=0; y[0]=0; z[0]=sigsa;
  for(i=1;i<N;i++)
  {
     x[i]=x[i-1];
     y[i]=y[i-1];
     z[i]=z[i-1]+4;
  }   */
    
 
  fp=fopen("Initconf","r");
  for(i=0;i<N;i++)
  {
     fscanf(fp,"%lf %lf %lf",&x[i],&y[i],&z[i]);
  }
  fclose(fp); 
 
 
///////////////////////////////////////////////////
  if(1!=cart2dof()) printf("Error initial configuration");
 
  /* initialize velocities */
  vxsum=vysum=vzsum=0;
  for(i=0;i<N;i++) 
  {
    vx[i]=sqrt(T)*gasdev2();
    vy[i]=sqrt(T)*gasdev2();
    vz[i]=sqrt(T)*gasdev2();
    vxsum+=vx[i]; 
    vysum+=vy[i];
    vzsum+=vz[i];
  }
  for(i=0;i<N;i++) 
  {
    vx[i]-=vxsum/N;
    vy[i]-=vysum/N;
    vz[i]-=vzsum/N;
  }
 
//  histoe(-1,0);
//  historg(-1,0);
  cont(-1);
  //////////////////////////////////////////////
  Ekin=0;
  for(i=0;i<N;i++) Ekin+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  Ekin*=0.5;
 
  for(i=0;i<N;i++)
  {
    fx[i]=fy[i]=fz[i]=0; 
  }
 
  Epot=(Ebon=bond())+(Eben=bend())+(Econ=cont(0));
 
  for(i=0;i<N;i++) 
  {
    frdx[i]=gasdev2()*tconst;
    frdy[i]=gasdev2()*tconst;
    frdz[i]=gasdev2()*tconst;
  }
}
 
/****************************************************************************/
double gasdev2() 
{
  double U1, U2, Gauss;
 
  if(phase==0) 
  { 
     do
     {
        U1=irand(&iseed);
     }while(U1 <= 0);     
     U2=irand(&iseed);
 
     V1=sqrt(-2*log(U1));
     V2=pi2*U2;
     Gauss=V1*sin(V2);
  }
  else
  {
     Gauss=V1*cos(V2);    
  }
 
  phase=1-phase;
  return Gauss;
 
}
/****************************************************************************/
#define mm 714025
#define r_ia 1366
#define r_ic 150889
#define rmm 1.0/mm
 
double irand(long int *iseed)
{
        int j;
    double ran;
    if(*iseed<0)
    {   
        *iseed=fmod(r_ic-(*iseed),mm);
        for(j=1;j<=97; j++)
        {
            *iseed=fmod(r_ia*(*iseed)+r_ic,mm);
            ir[j]=*iseed;
        }    
        *iseed=fmod(r_ia*(*iseed)+r_ic,mm);
        iy=*iseed;
    }
    jjj=1+(97*iy)/mm;
    iy = ir[jjj];
    ran=iy*rmm;
    *iseed=fmod(r_ia*(*iseed)+r_ic,mm);
    ir[jjj] =*iseed;
    return ran;
}
 
#undef mm
#undef r_ia
#undef r_ic
#undef rmm
 
void printhead() 
{
  printf("\n***************************************************\n");
  printf("N %i MDCYC %i MDCON %i NTHERM %i\n",N,MDCYC,MDCON,NTHERM);
  printf("T %f\n",T);
  printf("Interaction parameters\n");
  printf("kbon %f kth %f krep %f \n",kbon,kth,krep);
  printf("MD parameters\n");
  printf("tau %f dt %f gam %f\n",tau,dt,gam);
  printf("c1 %f c2 %f c3 %f\n",c1,c2,c3);
  printf("***************************************************\n");
  fflush(0);
  return;
}
 
void runtime(char *fn,long it,double o[],int n) 
{
  int i;
  FILE *fp;
   
  fp=fopen(fn,"a");
  fprintf(fp,"%li  ",it);
  for(i=0;i<n;i++) fprintf(fp,"%.4f  ",o[i]);
  fprintf(fp,"\n");
  fclose(fp);
  return;
} 
