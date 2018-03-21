double pscal(int N,double *vec1,double *vec2){
  int i;
  double res=0.0;
  for(i=0;i<N;i++){
    res+=vec1[i]*vec2[i];
  }
  return res;
}
