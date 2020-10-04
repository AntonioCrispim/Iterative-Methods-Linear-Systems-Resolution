#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#include<string.h>

#define n 63

int preencheMatriz(double **A, int k){
     
     int i,j;
     
     for (i=0; i<k; i++){
               for(j=0; j<k;j++){
                                 if (i==j){
                                           A[i][j]= 2.0;
                                 }else{
                                       A[i][j]=0.0;
                                 }
               }
                           
     }
                           
      for (i=0; i<k; i++){
                          A[i][i+1]=-1.0;
                                 if(i!=0){
                                          A[i][i-1]=-1.0;
                                 }
      }   
      return 0;          
}

int preencheVetor(double *v, int l){
    
    int j;
    double pi;
                         
    pi = 2.0*asin(1);
                         
    for(j=0; j<=l+1; j++){
                          v[j] = sin((16*j*pi)/(l+1))+sin((40*j*pi)/(l+1));
                          }                  
    return 0;                         
}

int GaussSeidel(double **A ,double *f, double *V, int k, int l ){ 
      
      int i, j, cont;
      double *Aux;
      cont = 0;
      
      Aux = malloc(n * sizeof(double));
      if(Aux == NULL){
             printf("Erro de alocacao de memoria\n");
             exit(1);
      }
      
      do{                 
          for(i=0;i<l;i++){                                             
                           Aux[i]=0.0;
                           for(j=0;j<l;j++){
                                            if(i!=j){
                                                     Aux[i] = Aux[i]+ A[i][j]*V[j];
                                            }
                           }
                           V[i]=(f[i]-Aux[i])/A[i][i];                                 
          }
          cont = cont + 1;
      }while(cont<k); 
      return 0;
}

int Residuo(double **A, double *V, double *f, double *r, int l){
    
    int i, j;
    double *Aux;
    
    Aux = malloc(n * sizeof(double));
      if(Aux == NULL){
             printf("Erro de alocacao de memoria\n");
             exit(1);
      }
    
    for(i=0; i<l; i++){
             Aux[i]=0.0;
             for(j=0; j<l; j++){
                      Aux[i]= Aux[i] + A[i][j]*V[j];
             }
             r[i]= f[i] - Aux[i];    
    }
}

int InterpolacaoRestricao(double **I21, double **I12, double *r, int k, int IR){
    
   int i, j;
   double *Aux;
   
   Aux = malloc(n * sizeof(double));
   if(Aux == NULL){
             printf("Erro de alocacao de memoria\n");
             exit(1);
   }
   
   for(j=0; j<(k-1)/2 ; j++){ 
            for (i=0; i<k-1; i++){
                                 if (i < (2*(j+1)+1) && i>=2*j && i != (2*(j+1)-1)){
                                                                                    I21[i][j]= 1.0;
                                 }else if (i < (2*(j+1)+1) && i>=2*j &&  i == (2*(j+1)-1)){
                                                                                          I21[i][j]=2.0;
                                       }else{
                                             I21[i][j]=0.0;
                                       }   
            }
                           
   }
   
   for(i=0; i<10; i++){                           
           printf("\n");
           for(j=0;j<(10-1)/2;j++){
                                   printf("%lf  ",I21[i][j]);
     }  
   }
     
   for(i=0; i<(k-1) ; i++){ 
            for (j=0; j<(k-1)/2; j++){
                                      I12[j][i] = I21[i][j];
            }
   } 
   
   if(IR==0){
             for(i=0; i<(k-1) ; i++){
                      Aux[i]=0.0; 
                      for (j=0; j<(k-1)/2; j++){
                                Aux[i]= Aux[i] + I12[i][j]*r[j];
                      }
                      r[i] = Aux[i];
                      printf("%lf", r[i]);    
             }
   }else if (IR==1){
                    for(i=0; i<(k-1)/2 ; i++){
                             Aux[i]=0.0; 
                             for (j=0; j<(k-1); j++){
                                       Aux[i]= Aux[i] + I21[i][j]*r[j];
                             }
                             r[i] = Aux[i];
                             //printf("%lf", r[i]);    
                    }
   }
   
   return 0;   
}

int Correcao(double *V, double *E, int l){
    
    int i;
    
    for(i=0; i<l; i++){
             V[i] = V[i] + E[i];
             //printf("E=%lf, V=%lf ", E[i], V[i]);
    }
}

void plotGraphic(const char *gnucommand){
    char syscommand[1024];

    sprintf(syscommand, "echo \%s\ | gnuplot -persist", gnucommand);
    system(syscommand);
}

int main(void){
     double **A, *f, *V, *r, **I12, **I21, *E, somaE, Er, max, h;
     int i, j, k, l, cont;
     FILE  *GrafG1, *GrafG2;
     
         
  
  GrafG1 = fopen("GrafG1.txt", "w");
  if (GrafG1 == NULL){
    printf("erro na apertura do arquivo GrafG1.txt\n");
    exit(EXIT_FAILURE);
  }
  
  GrafG2 = fopen("GrafG2.txt", "w");
  if (GrafG2 == NULL){
    printf("erro na apertura do arquivo GrafG2.txt\n");
    exit(EXIT_FAILURE);
  }
  
  A = malloc(n * sizeof(double *));              
  if(A == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<n; i++){
    A[i] = malloc(n * sizeof(double));
    if(A[i] == NULL){
      printf("Erro de alocacao de memoria\n");
      exit(1);
    }
  }
  
  f = malloc(n * sizeof(double));
  if(f == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  V = malloc(n * sizeof(double));
  if(V == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
  
  r = malloc(n * sizeof(double));
  if(r == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
  
  I12 = malloc((n-1) * sizeof(double *));              
  if(I12 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<(n-1)/2; i++){
    I12[i] = malloc(n * sizeof(double));
    if(I12[i] == NULL){
      printf("Erro de alocacao de memoria\n");
      exit(1);
    }
  }
  
  I21 = malloc((n-1)/2 * sizeof(double *));              
  if(I21 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<(n-1); i++){
      I21[i] = malloc(n * sizeof(double));
      if(I21[i] == NULL){
         printf("Erro de alocacao de memoria\n");
         exit(1);
    }
  }

  E = malloc(n * sizeof(double));
  if(E == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
                   
  
  preencheMatriz(A, n);
  preencheVetor(V, n); 
  h=0.0;
  for(i=0; i<63; i++){
           
            h = i*1.0/(63-1);               
            fprintf(GrafG1,"%lf %lf\n", h, V[i]);
  }
  

  
  for(i=0; i<n; i++){
           
           f[i]=0.0;
           E[i]=0.0;
  }

            GaussSeidel(A, f, V, 3, 10);
            Residuo(A, V, f, r, 10);
            InterpolacaoRestricao(I21, I12, r, 10, 1);
            GaussSeidel(A, r, E, 3, 5);
            InterpolacaoRestricao(I21, I12, E, 10, 0);
            Correcao(V, E, n);
  
   for(i=0; i<63; i++){
           
            h = i*1.0/(63-1);               
            fprintf(GrafG2,"%lf %lf\n", h, V[i]);
   }
  
  fclose(GrafG1);
  fclose(GrafG2);
  
 /* free(A);
  free(f);
  free(V);
  free(r);
  free(E);
  free(I21);
  free(I12);
  */
  A=NULL;
  f=NULL;
  V=NULL;
  r=NULL;
  E=NULL;
  I21=NULL;
  I12=NULL;
  
  plotGraphic("set title 'Bigrid'; set xlabel 'h'; set ylabel 'u'; set xrange [0:1]; set yrange [-2:2]; plot 'GrafG1.txt' title 'Sem ajustamento' lt 1 w lp; pause 100");
  
  
return 0;
}  


