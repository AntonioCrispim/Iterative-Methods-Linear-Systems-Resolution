#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>

#define n 3

//double A[n][n]={{6.0, 2.0, -3.0}, {-1.0, 8.0, 3.0}, {1.0, 4.0, 12.0}}; double B[n]={5.0, -10.0, 12.0}; double X[n]={1.0, 1.0, 1.0};
double A[n][n]={{2.0, 0.8, -0.6}, {0.8, 1.0, 0.1}, {-0.6, 0.1, 4.0}}; double B[n]={1.0, -1.0, 1.0}; double X[n]={1.0, 1.0, 1.0};

void diagonalInversa(int N, double **D){                    /* A partir da matriz A obtém a inversa da diagonal de A */
     int i, j;

     for (i=0; i<N; i++){
               for (j=0; j<N; j++){
                         if(i==j){
                                  D[i][j] = A[i][j];
                                  if(D[i][j]!=0.0){
                                                   D[i][j]= 1.0/D[i][j];
                                  }
                         }else{
                               D[i][j] = 0.0;
                         }
               }
     }
}
void matrizLU(int N, double **LU){   /* A partir da matriz A obtém a diagonal superior somada a diagonal inferior de A */
     int i, j;

      for (i=0; i<N; i++){
                for (j=0; j<N; j++){
                          if (i==j){
                                    LU[i][j]= 0.0;
                          }
                          if (i>j){
                                   LU[i][j]= -A[i][j];
                          }
                          if (i<j){
                                   LU[i][j]= -A[i][j];
                          }
                }
      }
}

int main(void){
     double **D, **LU, **M, *N, *Aux, *El, min,p, e;
     int i, j, k, cont;

  D = malloc(n * sizeof(double *));              /* Aloca na memoria as ordens da matrizes */
  if(D == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<n; i++){
    D[i] = malloc(n * sizeof(double));
    if(D[i] == NULL){
      printf("Erro de alocacao de memoria\n");
      exit(1);
    }
  }

  LU = malloc(n * sizeof(double *));
  if(LU == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<n; i++){
    LU[i] = malloc(n * sizeof(double));
    if(LU[i] == NULL){
      printf("Erro de alocacao de memoria\n");
      exit(1);
    }
  }

  M = malloc(n * sizeof(double *));
  if(M == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  for(i=0; i<n; i++){
    M[i] = malloc(n * sizeof(double));
    if(M[i] == NULL){
      printf("Erro de alocacao de memoria\n");
      exit(1);
    }
  }

  N = malloc(n * sizeof(double));
  if(N == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  Aux = malloc(n * sizeof(double));
  if(Aux == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  El = malloc(n * sizeof(double));
  if(El == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

     diagonalInversa(n, D);              /*Chama o procedimento que calcula a inversa*/
     matrizLU( n, LU);                   /*Chama o procedimento que calcula as diogonais*/

       /*for(i=0; i<n; i++){
                 for(j=0; j<n; j++){
                          printf("D[%d][%d] = %lf \n", i, j, D[i][j]);
                          }
              } */

           /*   for(i=0; i<n; i++){
                 for(j=0; j<n; j++){
                          printf("LU[%d][%d] = %lf \n", i, j, LU[i][j]);
                          }
              }*/

     for(i=0; i<n; i++){                        /*Calcula a matriz M=D*LU.  Onde D é a inversa da diagonal e LU é a diagonal superior somada a inferior negativa*/
              for(j=0; j<n; j++){
                       M[i][j]=0.0;
                       for(k=0; k<n; k++){
                                M[i][j]=M[i][j] + D[i][k]*LU[k][j];
                       }
              }
     }

     for(i=0; i<n; i++){
                 for(j=0; j<n; j++){
                          printf("M[%d][%d] = %lf \n", i, j, M[i][j]);
                          }
              }

     printf("\n");


     for(i=0; i<n; i++){                    /* Cálcula N= D*B */
              N[i]=0.0;
              for(k=0; k<n; k++){
                       N[i]= N[i] + D[i][k]*B[k];
              }
     }

      for(i=0; i<n; i++){
               printf("N[%d] = %lf \n", i, N[i]);

               }
               printf("\n");

     printf("Qual precisao?");
     scanf("%lf",&p);
     k=0;

    do{                /*Faz as iteracoes X = M*X + N*/

    for(i=0; i<n; i++){
             Aux[i]=0.0;
             for(j=0; j<n; j++){
                      Aux[i]= Aux[i] + M[i][j]*X[j];
             }
             Aux[i]= Aux[i]+ N[i];
    }
    
    for(i=0; i<n; i++){
             X[i] = Aux[i];
                     }

    for(i=0; i<n; i++){
             Aux[i] = 0.0;
             for(j=0; j<n; j++){
                      Aux[i] = Aux[i] + A[i][j]*X[j];
                             }
                      El[i] = abs(Aux[i]-B[i]);
                            //printf("%lf \n", El[i]);
                      }

    
    for(i=0; i<n; i++){
             if(El[i]<=p){
                          k=k+1;
                         }
                      }
    }while(k<3);

     for(i=0; i<n; i++){                  /*Mostra os valores de xi*/

              printf("X%d = %lf \n", i, X[i]);

              }//*/
      scanf("%d", cont);

return 0;
}

