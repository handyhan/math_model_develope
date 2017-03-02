void gradient_descent(int monitor,double (*omega)(),int dimension,double *value,double *Qvector,double precision)
{
int i,isymm,j,k,l,time;
void symmetrize_normalize();
double dQ,change,oldomega,newomega,RANDOM(),*newQ,*oldQ,*dvector(),**Qmatrix,**dmatrix();
    newQ=dvector(1,dimension); for(i=1;i<=dimension;i++) newQ[i]=Qvector[i];
    dQ=0.01;
    Qmatrix=dmatrix(0,kmax,0,kmax);
    oldQ=dvector(1,dimension);
    oldomega=omega(newQ);
    while(dQ>precision) {
        for(time=1;time<=(dimension*dimension);time++) {
                for(i=0;i<=kmax;i++) for(j=0;j<=kmax;j++) Qmatrix[i][j]=newQ[1+j+i*(kmax+1)];
                symmetrize_normalize(Qmatrix,pk,kmax);
                for(i=0;i<=kmax;i++) for(j=0;j<=kmax;j++) newQ[1+j+i*(kmax+1)]=Qmatrix[i][j];
            oldomega=omega(newQ);
            for(i=1;i<=dimension;i++) oldQ[i]=newQ[i];
            i=1+(int)(RANDOM()*dimension);
            /* printf("i=%d: ",i); */
            l=(i-1)/(kmax+1); k=i-1-l*(kmax+1);
            isymm=1+l+k*(kmax+1); /* make moves that preserve symmetry */
            newQ[i]=newQ[isymm]=newQ[i]+dQ;
                for(i=0;i<=kmax;i++) for(j=0;j<=kmax;j++) Qmatrix[i][j]=newQ[1+j+i*(kmax+1)];
                symmetrize_normalize(Qmatrix,pk,kmax);
                for(i=0;i<=kmax;i++) for(j=0;j<=kmax;j++) newQ[1+j+i*(kmax+1)]=Qmatrix[i][j];
            newomega=omega(newQ);
            newQ[i]=newQ[isymm]=oldQ[i];
            change=RANDOM()*dQ*(newomega-oldomega);
            newQ[i]=newQ[isymm]=newQ[i]-change;
                for(i=0;i<=kmax;i++) for(j=0;j<=kmax;j++) Qmatrix[i][j]=newQ[1+j+i*(kmax+1)];
                symmetrize_normalize(Qmatrix,pk,kmax);
                for(i=0;i<=kmax;i++) for(j=0;j<=kmax;j++) newQ[1+j+i*(kmax+1)]=Qmatrix[i][j];
            newomega=omega(newQ);
            if(newomega>=oldomega) {                          // if change is bad
                newQ[i]=newQ[isymm]=oldQ[i];              // reject change
                dQ=0.2*dQ;                                             // and make step smaller
                }
            else dQ=1.05*dQ;                                             // if change is good increases the step
            if(monitor) printf(">> value=%.12lf\n",newomega);
            }
        if(time==(dimension*dimension)) dQ=0.5*dQ;
        }
    for(i=1;i<=dimension;i++) Qvector[i]=newQ[i];
    free_dvector(newQ,1,dimension);
    free_dmatrix(Qmatrix,0,kmax,0,kmax);
    *value=oldomega;
}
