#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define TRUE 1
#define FALSE 0

/* -------------------------------------------------------------------------------------------------------------- */
void diagonalization(double **A,int n,double *lambda,double **y);
/* -------------------------------------------------------------------------------------------------------------- */
void davidson(int N,double **v,double a,double b,int nev,int first_ev){
  int i,j,k,l;

  double tol=1.0e-4;
  double etol=1.0e-6;
  int nloopmax=5000;
  int cvg; 
  double norm;  
  double **T;			/* reduced matrix */
  double **y;
  double *lambda;
  double **Ritz;
  double **r;          r=malloc(nev*sizeof(double*)); for(i=0;i<nev;i++) r[i]=malloc(N*sizeof(double));
  double *normr;   normr=malloc(nev*sizeof(double));
  double *ev;        ev=malloc(nev*sizeof(double));
  double *dev;      dev=malloc(nev*sizeof(double));  for(i=0;i<nev;i++) dev[i]=2*etol;

  double **t;         t=malloc(nev*sizeof(double*)); for(i=0;i<nev;i++) t[i]=malloc(N*sizeof(double));
  double **q;        q=malloc(nev*sizeof(double*)); for(i=0;i<nev;i++) q[i]=malloc(N*sizeof(double));
  int *q_insert;q_insert=malloc(nev*sizeof(int));
  double *beta;  if(nev>1) beta=malloc((nev-1)*sizeof(double));
  double *alpha;

  int loop;loop=1;
  FILE *log;log=fopen("debug.log","w+");fclose(log);
  int nq,iq,first_q;
  int nvecprev;

  
  int endloop=FALSE;
  int nvecini=  nvecini=first_ev+nev;
  int nvecmax; nvecmax=40;
  int nvec; nvec=nvecini;
  double **V;V=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) V[i]=malloc(N*sizeof(double));
  for(j=0;j<nvecini;j++) for(i=0;i<N;i++) V[j][i]=v[j][i];

  
  /* ----------------------
     Davidson loop 
  ------------------------ */
  while(endloop==FALSE){
    printf("----------------------------------------\n");
    /* ----------------------------------------------------------------------------------
       Computation of the reduced matrix T(nvec,nvec)=VT.A.V to be diagonalized in the 
       reduced subspace
       ------------------------------------------------------------------------------------ */
    T=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) T[i]=malloc(nvec*sizeof(double));
    for(i=0;i<nvec;i++){
      for(j=0;j<nvec;j++){
	T[i][j]=V[i][0]*(a*V[j][0]+b*V[j][1])
	  +V[i][N-1]*(a*V[j][N-1]+b*V[j][N-2]);
	for(l=1;l<N-1;l++) T[i][j]+=V[i][l]*(a*V[j][l]+b*(V[j][l-1]+V[j][l+1]));
      }
    }
    /* -------------------------------
       diagonalization of T(nvec,nvec)
       ------------------------------ */
    y=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) y[i]=malloc(nvec*sizeof(double));
    lambda=malloc(nvec*sizeof(double));
    diagonalization(T,nvec,lambda,y);
    for(i=0;i<nvec;i++) free(T[i]);free(T);
    for(i=0;i<nvec;i++) printf("%e ",lambda[i]);printf("\n");
    /* -----------------------------------------------------------------------
       expression of the |Ritz vectors> = A |y>
       nev is the number of eigenvalues to compute : nev=last_ev-first_ev+1
       There are nvec+nev Ritz vectors : 
           - nvec coming from the diagonalization of T; they are computed here
           - nev 
       ----------------------------------------------------------------------- */
    Ritz=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) Ritz[i]=malloc(N*sizeof(double));
    for(k=0;k<nvec;k++){
      for(l=0;l<N;l++) Ritz[k][l]=0.0;
      for(i=0;i<nvec;i++) 	for(l=0;l<N;l++) Ritz[k][l]+=V[i][l]*y[k][i];
    }
    for(i=0;i<nvec;i++) free(y[i]);free(y);
    /* ------------------------------------------------------------------
       computation of the residuals |r_k>=A |Ritz_k> - lambda_k |Ritz_k>
       only yhr norm of the eigenvalues we search are computed
       ------------------------------------------------------------------ */
    for(k=0;k<nev;k++){
      for(l=0;l<N;l++) {
	r[k][l]=(a-lambda[k+first_ev])*Ritz[k+first_ev][l];
	if(l>0) 	  r[k][l]+=b*Ritz[k+first_ev][l-1];
	if(l<N-1)   r[k][l]+=b*Ritz[k+first_ev][l+1];
      }
      normr[k]=0.0;
      for(l=0;l<N;l++)     normr[k]+=r[k][l]*r[k][l];
      normr[k]=pow(normr[k],.5);
    }

    for(k=0;k<nev;k++) {
      if(loop>1) dev[k]=ev[k]-lambda[k+first_ev];
      ev[k]=lambda[k+first_ev];
    }
    
    /* -------------------------------------------------------------------------
     print the logfile
    --------------------------------------------------------------------------- */
    printf("loop %d nvecmax= %d\n",loop,nvecmax);
    if(loop>1){
      printf("dev= ");for(k=0;k<nev;k++)  printf("%e ",dev[k]);printf("\n");
    }
    printf("normr= ");for(k=0;k<nev;k++)  printf("%e ",normr[k]);printf("\n");
    printf("lambda= ");for(k=0;k<nev;k++) printf("%e ",lambda[k+first_ev]);printf("\n");
    log=fopen("debug.log","a+");
    fprintf(log,"loop= %d lambda= ",loop);
    for(k=0;k<nev;k++) fprintf(log,"%e ",lambda[k+first_ev]);//fprintf(log,"\n");
    fprintf(log,"normr= ");for(k=0;k<nev;k++)  fprintf(log,"%e ",normr[k]);fprintf(log,"\n");
    fclose(log);    
    /* -------------------------------------------------------------------------
       check the convergence
    --------------------------------------------------------------------------- */
   cvg=TRUE;

    if(loop>1){
      i=0;
      while(i<nev && cvg==TRUE){
	if(normr[i]>tol) cvg=FALSE;
	//if(fabs(dev[i])>etol) cvg=FALSE;
	i++;
      }
    } else {
      cvg=FALSE;
    }

    if(cvg==FALSE){
      if(nvec+nev<=nvecmax){	/* one computes nev additional vectors to expand the Krylov's subspace */
	/* --------------------------------------------------------------------------------------
	   computation of |t_k>=M|r_k>
	   M is the preconditioning matrix. It is an approximation of (A-lambdaI)^{-1}
	   Here M=(D-lambda I)^{-1}  where D is the main diagonale of A -> Jacobi Preconditionner
	   ------------------------------------------------------------------------------------ */
	for(k=0;k<nev;k++){
	  for(i=0;i<N;i++) 	t[k][i]=r[k][i]/(a-lambda[k+first_ev]);

	  vec2=GramSchmidt(N,Ritz,nvec,t[k]);

	}



	/* --------------------------------------------------------------------------- 
	   
	                       GRAM-SCHMIDT

	   Gram-Schmidt process to insert |t> in the reduced subspace 
	   From Saad's book "Numerical Methods For Large Eigenvalue Problems", p12
	   --------------------------------------------------------------------------- */
	/* first, one removes the Ritz vectors from t[i] --> q[i] */
	alpha=malloc(nvec*sizeof(double));
	for(i=0;i<nev;i++){		/* loop over the new vectors to add */
	  for(k=0;k<nvec;k++){
	    alpha[k]=0.0;
	    for(l=0;l<N;l++) alpha[k]+=Ritz[k][l]*t[i][l];
	  }
	  for(l=0;l<N;l++){
	    q[i][l]=t[i][l];
	    for(k=0;k<nvec;k++) q[i][l]-=alpha[k]*Ritz[k][l];
	  }
	}
	free(alpha);
	/* one computes the norm of the first additional vector */
	nq=0;
	iq=0;
	while(nq==0 && iq<nev){
	  norm=0.0;
	  for(l=0;l<N;l++) norm+=q[iq][l]*q[iq][l];
	  norm=pow(norm,.5);
	  log=fopen("debug.log","a+");fprintf(log,"norm q[%d]=%g\n",iq,norm);fclose(log);
	  if(norm<1.0e-6) {
	    printf("1 - WARNING in Gram-Schmidt process : q[%d]= %g\n",iq,norm);
	    q_insert[iq]=FALSE;
	  } else {
	    first_q=iq;
	    q_insert[iq]=TRUE;
	    nq++;
	    for(l=0;l<N;l++) q[iq][l]/=norm;
	  }
	  iq++;
	}

	if(nq==0) {
	  printf("ERROR in GRAM-SCHMIDT: all <q|q>=0.0\n");

	  exit(0);
	}
	log=fopen("debug.log","a+");fprintf(log,"first q %d ",first_q);for(i=0;i<nev;i++) fprintf(log,"iq[%d]=%d ",i,q_insert[i]);printf("\n");
	for(i=first_q;i<nev;i++) fprintf(log,"ff %d ",i); fprintf(log,"\n");
	for(i=first_q+1;i<nev;i++) fprintf(log,"gg %d ",i); fprintf(log,"\n");fclose(log);
	
	for(i=first_q+1;i<nev;i++){ 
	  for(k=0;k<i;k++){
	    if(q_insert[k]==TRUE){
	      beta[k]=0.0;
	      for(l=0;l<N;l++) beta[k]+=q[k][l]*t[i][l];
	      for(l=0;l<N;l++){
		q[i][l]=t[i][l];
		for(j=0;j<i;j++) q[i][l]-=beta[j]*q[j][l];
	      }
	    }
	  }
	  norm=0.0;
	  for(l=0;l<N;l++) norm+=q[i][l]*q[i][l];
	  norm=pow(norm,.5);
	  log=fopen("debug.log","a+");fprintf(log,"norm q[%d]=%g\n",i,norm);fclose(log);
	  if(norm<1.0e-6) {
	    printf("2 - Error in Gram-Schmidt process : q[%d]= %g nq=%d\n",i,norm,nq);
	    q_insert[i]=FALSE;
	  } else {
	    nq++;
	    q_insert[i]=TRUE;
	    for(l=0;l<N;l++) q[i][l]/=norm;
	  }
	}
      

	printf("%d q to insert\n",nq);
	nvecprev=nvec;
	nvec=nvec+nq;
	for(i=0;i<nvecprev;i++) free(V[i]);free(V);
	V=malloc(nvec*sizeof(double*)); for(i=0;i<nvec;i++) V[i]=malloc(N*sizeof(double));
	for(k=0;k<nvecprev;k++) {
	  for(l=0;l<N;l++) V[k][l]=Ritz[k][l];
	}
	iq=0;
	for(k=0;k<nev;k++) {
	  if(q_insert[k]==TRUE) {
	    for(l=0;l<N;l++)	      V[nvecprev+iq][l]=q[k][l];
	    iq++;
	  }
	}



	for(i=0;i<nvecprev;i++) free(Ritz[i]);free(Ritz);
	/* /\* PRINT BLOCK *\/ */
	log=fopen("debug.log","a+");
	fprintf(log,"===================================\n");
	/* for(l=0;l<N;l++) { */
	/*   for(i=0;i<nvec;i++) printf("%12.6f ",V[i][l]); */
	/*   printf("\n"); */
	/* } */
	/* printf("\n"); */
	for(i=0;i<nvec;i++){
	  for(j=0;j<nvec;j++){
	    norm=0.0;
	    for(l=0;l<N;l++) norm+=V[i][l]*V[j][l];
	    fprintf(log,"%12.6f ",norm);
	  }
	  fprintf(log,"\n");
	}
	fprintf(log,"===================================\n");
	fclose(log);
	/* exit(0); */
	/* /\* END OF PRINT BLOCK *\/ */



      } else { // nvec+nev>nvecmax
	nvecprev=nvec;
	nvec=nvecini;
	for(i=0;i<nvecprev;i++) free(V[i]);free(V);
	V=malloc(nvec*sizeof(double*)); for(i=0;i<nvec;i++) V[i]=malloc(N*sizeof(double));
	for(k=0;k<nvec;k++) {
	  for(l=0;l<N;l++) V[k][l]=Ritz[k][l];
	}
	for(i=0;i<nvecprev;i++) free(Ritz[i]);free(Ritz);
      }
    } /* if(cvg==FALSE) */
    loop++;
    free(lambda);
    endloop=FALSE;
    if(loop>=nloopmax || cvg==TRUE) endloop=TRUE;
    /* end of Davidson loop */
  } 
  
}
