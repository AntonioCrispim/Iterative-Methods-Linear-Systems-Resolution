#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#include<string.h>

#define n 257


int preencheMatriz(double **A, int k){                        /* Gera a matriz tridiagonal representante do sistema */
     
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

int preencheVetor(double *v, int l){                            /* Gera o vetor com chute inicial da solução do sistema com dois componentes, um de baixa e um de alta frequência */
    
    int j;
    double pi;
                         
    pi = 2.0*asin(1);
                         
    for(j=0; j<l; j++){
                          v[j] = 0.5*(sin((90*(j+1)*pi)/(l+1)) + sin((210*(j+1)*pi)/(l+1)));
                          //printf("%lf  \n", v[j]);
                          }                  
    return 0;                         
}

int GaussSeidel(double **A ,double *f, double *V, int k, int l ){  /* Método de GaussSeidel utilizado para relaxação do sistema*/
      
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

int Residuo(double **A, double *V, double *f, double *r, int l){ /* Calcula o resíduo  r = f - Av */
    
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

int GeraMatrizInterpolacaoRestricao(double **I21, double **I12, int k){ /* Gera as matrizes de interpolação e restrição, e faz a interpolação e restrição*/ 
    
   int i, j;
  
   
   for(j=0; j<(k-1)/2 ; j++){                                                      /* Gera a matriz de interpolação */
            for (i=0; i<k; i++){
                                 if (i < (2*(j+1)+1) && i>=2*j && i != (2*(j+1)-1)){        /*  j=0  ->  0<= i <3  e para i=0 -> I21[0][0]= 1.0, para i=1 -> I21[1][0]= 2.0 ...*/    
                                                                                    I21[i][j]= 1.0;
                                 }else if (i < (2*(j+1)+1) && i>=2*j &&  i == (2*(j+1)-1)){
                                                                                          I21[i][j]=2.0;
                                       }else{
                                             I21[i][j]=0.0;
                                       }   
            }
                           
   }
   
  /*for(i=0; i<k; i++){                           
           printf("\n");
           for(j=0;j<(k-1)/2;j++){
                                   printf("%lf  ",I21[i][j]);
     }  
   }*/
    
   for(i=0; i< k ; i++){                            /* Gera a matriz de restrição que é a transposta da de interpolação */ 
            for (j=0; j<(k-1)/2; j++){
                                      I12[j][i] = I21[i][j];
            }
   } 
   
   /* for(i=0; i<(k-1)/2; i++){                           
           printf("\n");
           for(j=0;j<k;j++){
                                   printf("%lf  ",I12[i][j]);
     }  
   } */
   
   return 0; 
   
     
}



int Correcao(double *V, double *E, int l){  /* Faz  a correção do chute inicial para aproximar da solução */
    
    int i;
    
    for(i=0; i<l; i++){
             V[i] = V[i] + E[i];
             //printf("E=%lf, V=%lf ", E[i], V[i]);
    }
}

int Redimensionamento(double **A1, double **A2, double **I12, double **I21, double *r1, double *r2, int k, int IR){        /* Redimensiona a matriz do sistema, fazendo R*A1*P=A2 ou I12*A1*I21=A2  e faz interpolação I21*r1 ou restrição I12*r2 */
    
   int i, j, m;
   double **Aux, *Aux1, *Aux2;
   
   Aux1 = malloc(k * sizeof(double));
   if(Aux1 == NULL){
             printf("Erro de alocacao de memoria\n");
             exit(1);
   }
   
   Aux2 = malloc((k-1)/2 * sizeof(double));
   if(Aux2 == NULL){
             printf("Erro de alocacao de memoria\n");
             exit(1);
   }
   
   Aux = malloc(k * sizeof(double *));              
         if(Aux == NULL){
         printf("Erro de alocacao de memoria\n");
         exit(1);
        }

   for(i=0; i<k; i++){
      Aux[i] = malloc(k * sizeof(double));
         if(Aux[i] == NULL){
          printf("Erro de alocacao de memoria\n");
          exit(1);
        }
   }
   
    if(IR==0){                                              /* Caso IR=0 na chamada da subrotina faz a interpolação */
             for(i=0; i<k ; i++){
                      Aux1[i]=0.0; 
                      for (j=0; j<(k-1)/2; j++){
                                Aux1[i]= Aux1[i] + 0.5*I21[i][j]*r2[j];
                                
                      }
                      r1[i] = Aux1[i];
                      //printf("%lf \n", r1[i]); 
                      
       
                      
                         
             }
   }else if (IR==1){
                    for(i=0; i<(k-1)/2 ; i++){              /* Caso IR=1 na chamada da subrotina faz a restrição */
                             Aux2[i]=0.0; 
                             for (j=0; j<k; j++){
                                       Aux2[i]= Aux2[i] + 0.25*I12[i][j]*r1[j];
                             }
                             r2[i] = Aux2[i];
                             //printf("%lf", r2[i]);    
                    }
   }
   
   for(i=0; i<(k-1)/2 ; i++){ 
             for (j=0; j<k; j++){
                       Aux[i][j]=0.0;
                       for (m=0; m<k; m++){          
                          Aux[i][j]= Aux[i][j] + I12[i][m]*A1[m][j];
                       }
             }     
    }
   
    /*for(i=0; i<(k-1)/2; i++){                           
           printf("\n");
           for(j=0;j<k;j++){
                                   printf("%lf  ",Aux[i][j]);
     }  
   }*/
   
    for(i=0; i<(k-1)/2 ; i++){ 
             for (j=0; j<(k-1)/2; j++){
                       A2[i][j]=0.0;
                       for (m=0; m<k; m++){          
                          A2[i][j]= A2[i][j] + Aux[i][m]*I21[m][j];
                       }
             }     
    }
    
   /* for(i=0; i<(k-1)/2; i++){                           
           printf("\n");
           for(j=0;j<(k-1)/2;j++){
                                   printf("%lf  ",A2[i][j]);
     }  
   }*/
}    

/*Rotina que faz técnica multigrid recursivamente*/
void ZeraV(double *V, int v){
     
     int i;
     
     for(i=0; i<v; i++){
           V[i]=0.0;
     }
}
void MG(double **A1, double **A2, double **I12, double **I21, double *f, double *V, double *r1, double *r2, double *E1, double *E2, int k, int i){
    
       while(i<1){
                  if(k<2){
                         
                          GaussSeidel(A1, f, V, 3, n);
                          Residuo(A1, V, f, r1, n);
                          GeraMatrizInterpolacaoRestricao(I21, I12, n);
                          Redimensionamento(A1, A2, I12, I21, r1, r2, n, 1);
                          ZeraV(E2, (n-1)/2);
                          k = k + 1;         
                          MG(A1, A2, I12, I21, f, V, r1, r2, E1, E2, k, i);
                         
                          Redimensionamento(A1, A2, I12, I21, E1, E2, n, 0);
                          Correcao(V, E1, n);
                          }else{
                              
                              GaussSeidel(A2, r2, E2, 3, (n-1)/2);
                              i = i + 1;
                 }
       }
}

void plotGraphic(const char *gnucommand){
    char syscommand[1024];

    sprintf(syscommand, "echo \%s\ | gnuplot -persist", gnucommand);
    system(syscommand);
}

int main(void){
     double **A1, **A2, *f, *V, *r1, *r2, *E1, *E2, **I12, **I21, somaE, Er, max, h;
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
  
 /* alocacao2(A1, n, n);
  alocacao2(A2, n, n);
  alocacao2(I21, n-1, (n-1)/2);
  alocacao2(I12, (n-1)/2, n-1);
  alocacao1(f, n);
  alocacao1(V1, n);
  alocacao1(r1, n);
  alocacao1(E1, n);
  */
  
  A1 = malloc(n * sizeof(double *));              
  if(A1 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<n; i++){
    A1[i] = malloc(n * sizeof(double));
    if(A1[i] == NULL){
      printf("Erro de alocacao de memoria\n");
      exit(1);
    }
  }
  
   A2 = malloc((n-1)/2 * sizeof(double *));              
  if(A2 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<(n-1)/2; i++){
    A2[i] = malloc(n * sizeof(double));
    if(A2[i] == NULL){
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
  
  r1 = malloc(n * sizeof(double));
  if(r1 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
  
   r2 = malloc((n-1)/2 * sizeof(double));
  if(r2 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
  
  I12 = malloc((n-1)/2 * sizeof(double *));              
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
  
  I21 = malloc(n * sizeof(double *));              
  if(I21 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i< n; i++){
      I21[i] = malloc((n-1)/2 * sizeof(double));
      if(I21[i] == NULL){
         printf("Erro de alocacao de memoria\n");
         exit(1);
    }
  }

  E1 = malloc(n * sizeof(double));
  if(E1 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
  
   E2 = malloc((n-1)/2 * sizeof(double));
  if(E2 == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
                   
  
  preencheMatriz(A1, n);
  preencheVetor(V, n); 
  h=0.0;
  for(i=0; i<n; i++){
           
            h = i*1.0/(n-1);               
            fprintf(GrafG1,"%lf %lf\n", h, V[i]);
  }
  fclose(GrafG1);


  
  for(i=0; i<(n-1)/2; i++){
           E2[i]=0.0;
  }
  for(i=0; i<n; i++){
           f[i]=0.0;
  }
  
  MG(A1, A2, I12, I21, f, V, r1, r2, E1, E2, 0, 0);
   
   for(i=0; i<n; i++){
           
            h = i*1.0/(n-1);               
            fprintf(GrafG2,"%lf %lf\n", h, V[i]);
   }
  
  
  fclose(GrafG2);
  
  /*free(A);
  free(f);
  free(V);
  free(r);
  free(E);
  free(I21);
  free(I12);
   
  A=NULL;
  f=NULL;
  V=NULL;
  r=NULL;
  E=NULL;
  I21=NULL;
  I12=NULL;
  */
  
  
   plotGraphic("set title 'Bigrid'; set xlabel ''; set ylabel ''; set xrange [0:0.1]; set yrange [-0.1:0.1]; plot 'GrafG1.txt' title 'Sem ajustamento' lt 1 w lp,\ 'GrafG2.txt' title 'Com ajustamento' lt 2 w lp; pause 100");
  
  
 scanf("%d",i); 
return 0;
}  


