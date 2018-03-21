#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


double dot(int N,double *vec1,double *vec2);
void davidson(int N,double **v,double a,double b,int nev,int first_ev);
double **GramSchmidt(int N,double **set0,int nvec0,double *set1);

#define NCHAR 1024

#define DELIM   " \n" 
#define MAXWORD 50
#define MAXLEN  512
typedef struct {
  char line[NCHAR];
  int ntokens;
  char words[MAXWORD][MAXLEN];
} Line;
Line *read_line(char *line){
  Line *sline;sline=malloc(sizeof(Line));   sprintf(sline->line,"%s",line);
  char *sep=NULL;  sep=strtok(line,DELIM);

  sline->ntokens=0;
  while (sep != NULL){
    strcpy(sline->words[sline->ntokens++], sep);
    sep = strtok(NULL, DELIM);
  }
  free(sep);
  return sline;
}


//#include <time.h>
double  *set_random_vector(int N,int init){
  double *vec; vec=malloc(N*sizeof(double));
  double norm;
  int i;
  srand(init); // initialisation de rand
  for(i=0;i<N;i++)     vec[i]=1.0*rand()/RAND_MAX;
  norm=dot(N,vec,vec);  norm=pow(norm,.5);  for(i=0;i<N;i++)     vec[i]/=norm;
  norm=dot(N,vec,vec);   norm=pow(norm,.5);
  return vec;
}


/* --------------------------------------------------------------------------------------------------------------
   Davidson algorithm. Variable names according Saad's book "Numerical Methods For Large Eigenvalue Problems", p203
   -------------------------------------------------------------------------------------------------------------- */
int main(int nargument,char **argument){

  /* INPUT
     - dimension N of the function/vector to compute and given by a the partial differential equation d^2Phi/dx^2 + gamma(x) dPhi/dx + V(x).Phi(x) = lambda.Phi(x)
    -  vectors V(x) (potential for example) and gamma(x) (damping for example)
    - initial guess vector
  */

  int N=200;
  double L=M_PI;
  double dx,dx2;
  dx=L/(N+1);
  dx2=dx*dx;
  double a,b;a=1.0/dx2;b=-.5/dx2;
  int nev;nev=4;
  int first_ev;first_ev=0;
  int nvecini;nvecini=first_ev+nev;

  int i,j,k,l;
  double **v;v=malloc(nvecini*sizeof(double*)); for(i=0;i<nvecini;i++) v[i]=malloc(N*sizeof(double));
  double norm;
  double *alpha;





  int nvec0,nvec1;nvec0=1;nvec1=1;
  double **set0,*set1;
  set0=malloc(nvec0*sizeof(double*));
  double **vec2;
  for(j=0;j<nvec0;j++)      set0[j]=set_random_vector(N,j+1);


  int loop;
  for(loop=1;loop<10;loop++){
    set1=set_random_vector(N,1000*(j+3*loop));
    vec2=GramSchmidt(N,set0,nvec0,set1);
    free(set1);
    for(j=0;j<nvec0;j++)      free(set0[j]);free(set0); nvec0++;   set0=malloc(nvec0*sizeof(double*));
    for(j=0;j<nvec0;j++){
      set0[j]=malloc(N*sizeof(double));
      for(l=0;l<N;l++) set0[j][l]=vec2[j][l];
      free(vec2[j]);
    }
    free(vec2);
    
    for(i=0;i<nvec0;i++){
      for(j=0;j<nvec0;j++){
	printf("%12.2e ",dot(N,set0[i],set0[j]));
      }
      printf("\n");
    }
  }
  
  
  exit(0);
  
  

  /* --------------------------------------------------------
     first, one initializes the first vectors of the subspace
     In the present version, They are set randomly.
     They are also normalized.
  ----------------------------------------------------------*/
  for(j=0;j<nvecini;j++){
    norm=0.0;
    for(i=0;i<N;i++) {
      v[j][i]=1.0*rand()/RAND_MAX;
      norm+=v[j][i]*v[j][i];
    }
    norm=pow(norm,.5);
    for(i=0;i<N;i++)     v[j][i]/=norm;
    norm=0.0;  for(i=0;i<N;i++)     norm+=v[j][i]*v[j][i];
    norm=pow(norm,.5);
  }
  if(nvecini>1){		/* Gram-Schmidt process if the iitial size of the subspace is larger than 1*/
    alpha=malloc((nvecini-1)*sizeof(double));
    for(j=1;j<nvecini;j++){
      for(k=0;k<j;k++) {
	alpha[k]=0.0;
	for(l=0;l<N;l++) alpha[k]+=v[k][l]*v[j][l];
      }
      for(k=0;k<j;k++) {
	for(l=0;l<N;l++) v[j][l]-=alpha[k]*v[k][l];
      }
      norm=0.0;
      for(l=0;l<N;l++) norm+=v[j][l]*v[j][l];
      norm=pow(norm,.5);
      for(l=0;l<N;l++) v[j][l]/=norm;
    }
  }
  free(alpha);


  davidson(N,v,a,b,nev,first_ev);

  
  printf("JOB DONE !\n");

  return 0;
}
