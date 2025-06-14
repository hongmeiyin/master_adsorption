# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
 
# define N   20
int main(void)
{
  int i, j,n=0; 
  double  k=0,l=0,m=6;
 
  FILE *fp;
 
  fp=fopen("nsur6","w");
  for(i=-N;i<=N;i++) for (j=-N;j<N;j++)
  { 
    n=n+1;
    k=i*(7.5); l=j*(-7.5);
    fprintf(fp,"%8.3lf %8.3lf %8.3lf\n",k,l,m); 
  }
  fclose(fp);
  return 0;
}
