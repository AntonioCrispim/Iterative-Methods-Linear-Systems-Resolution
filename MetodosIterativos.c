#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#include<string.h>

#define n 20


//double A[n][n]={{2.0, 0.8, -0.6}, {0.8, 1.0, 0.1}, {-0.6, 0.1, 4.0}}; double B[n]={1.0, -1.0, 1.0}; double X[n]={0.0, 0.0, 0.0};

void plotGraphic(const char *gnucommand){
    char syscommand[1024];

    sprintf(syscommand, "echo \%s\ | gnuplot -persist", gnucommand);
    system(syscommand);
}

int preencheMatriz(double **A){

     int i,j;

     for (i=0; i<n; i++){
               for(j=0; j<n;j++){
                        if (i==j){
                                  A[i][j]= 2.0;
                                  }else{
                                        A[i][j]=0.0;
                                        }
                                  }

                        }
      for (i=0; i<n; i++){
                  A[i][i+1]=-1.0;
                           if(i!=0){
                                    A[i][i-1]=-1.0;
                                    }
                           }
      return 0;
}

int preencheVetorB(double *B){
    int i;

    for(i=0; i<n; i++){
            if(i==0||i==n-1){
                             B[i]=1.0;
                             }else{
                                    B[i]= 0.0;
            }
    }
}



void diagonalInversa(int N, double **D, double **A){                    /* A partir da matriz A obtém a inversa da diagonal de A */
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
void matrizLU(int N, double **LU, double **A){   /* A partir da matriz A obtém a diagonal superior somada a diagonal inferior de A */
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
     double **A, *B, *X, **D, **LU, **M, *N, *Aux, *El, p, S, SomaE, E, max, w;
     int i, j, k, cont;
     FILE *GrafJ, *GrafG, *GrafS;


   A = malloc(n * sizeof(double *));              /* Aloca na memoria as ordens da matrizes */
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

    B = malloc(n * sizeof(double));
  if(B == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

      X = malloc(n * sizeof(double));
  if(X == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  GrafJ = fopen("GrafJ.txt", "w");
  if (GrafJ == NULL){
    printf("erro na apertura do arquivo GrafJ.txt\n");
    exit(EXIT_FAILURE);
  }

  GrafG = fopen("GrafG.txt", "w");
  if (GrafG == NULL){
    printf("erro na apertura do arquivo GrafG.txt\n");
    exit(EXIT_FAILURE);
  }

  GrafS = fopen("GrafS.txt", "w");
  if (GrafS == NULL){
    printf("erro na apertura do arquivo GrafS.txt\n");
    exit(EXIT_FAILURE);
  }

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

/*-------------------------------------------------Método de Jacobi----------------------------------------------------------------- */
     preencheMatriz(A);
     preencheVetorB(B);
     diagonalInversa(n, D, A);              /*Chama o procedimento que calcula a inversa*/
     matrizLU( n, LU, A);                   /*Chama o procedimento que calcula as diogonais*/

       /*for(i=0; i<n; i++){
                 for(j=0; j<n; j++){
                          printf("D[%d][%d] = %lf \n", i, j, D[i][j]);
                          }
              }

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

    /* for(i=0; i<n; i++){
                 for(j=0; j<n; j++){
                          printf("M[%d][%d] = %lf \n", i, j, M[i][j]);
                          }
              }

     printf("\n");
*/

     for(i=0; i<n; i++){                    /* Cálcula N= D*B */
              N[i]=0.0;
                       for(k=0; k<n; k++){
                                N[i]= N[i] + D[i][k]*B[k];
                       }
     }

      /*for(i=0; i<n; i++){
               printf("N[%d] = %lf \n", i, N[i]);

               }

               printf("\n");

               */

     printf("Qual precisao?");
     scanf("%lf",&p);
     printf("Qual fator de relaxamento?");
     scanf("%lf",&w);
     k=0;
     cont=0;

     for(i=0;i<n;i++){
                      X[i]=0.0;                                    /*Define vetor solução inicial*/
     }

     do{                /*Faz as iteracoes X = M*X + N*/
        cont = cont+1;
        for(i=0; i<n; i++){
                 Aux[i]=0.0;
                       for(j=0; j<n; j++){
                                Aux[i]= Aux[i] + M[i][j]*X[j];
                       }
                       Aux[i]= Aux[i]+ N[i];
        }

        for(i=0; i<n; i++){
                 X[i] = Aux[i];
                 //printf("%lf \n", X[i]);
                    }

        for(i=0; i<n; i++){
                 Aux[i] = 0.0;
                        for(j=0; j<n; j++){
                                Aux[i] = Aux[i] + A[i][j]*X[j];
                            }
                            //printf("%lf \n", Aux[i]);
                                El[i] = fabs(Aux[i]-B[i]);


        }
        j=0;
        for(i=0; i<n; i++){
                 max=El[j];
                 if(El[i]>max){
                              max = El[i];
                              j=i;
                              }
        }                                 /*grava dados no arquivo para fazer o grafico*/
          fprintf(GrafJ,"%d %lf\n", cont, max);
          //printf("%d", cont);
           //printf("%lf\n", E);

        }while(p<max);

        printf("Jacobi \n");

        for(i=0; i<n; i++){                  /*Mostra os valores de xi*/

              printf("X%d = %lf \n", i, X[i]);

              }//*/

/*-----------------------------------------------------Método de Gauss-Seidel--------------------------------------------------------*/
          for(i=0;i<n;i++){
                           X[i]=0.0;                                    /*Define vetor solução inicial*/
          }
       k=0;
       cont=0;
       E=0.0;

       do{
          cont= cont+1;
          for(i=0;i<n;i++){                                             /*Metodo de Gauss-Seidel*/
                           Aux[i]=0.0;
                           for(j=0;j<n;j++){
                                            if(i!=j){
                                            Aux[i] = Aux[i]+ A[i][j]*X[j];
                                            }
                           }
                           X[i]=(B[i]-Aux[i])/A[i][i];
                            //  printf("%lf \n", X[i]);
          }
           for(i=0; i<n; i++){                                             /*Obtem os erros*/
                 Aux[i] = 0.0;
                        for(j=0; j<n; j++){
                                Aux[i] = Aux[i] + A[i][j]*X[j];
                            }
                            //printf("%lf \n", Aux[i]);
                                El[i] = fabs(Aux[i]-B[i]);
                            //printf("%lf \n", El[i]);
        }
        j=0;
        for(i=0; i<n; i++){                                                /*Obtem o erro maximo*/
                 max=El[j];
                 if(El[i]>max){
                              max = El[i];
                              j=i;
                              }
        }
                                        /*grava dados no arquivo para fazer o grafico*/
          fprintf(GrafG,"%d %lf\n", cont, max);
              //printf("%d", cont);
               //printf("%lf\n", E);



        }while(p<max);

        printf("Gauss-Seidel \n");

        for(i=0; i<n; i++){                  /*Mostra os valores de xi*/

              printf("X%d = %lf \n", i, X[i]);

              }

/*----------------------------------------------------------Metodo de SOR--------------------------------------------------------------------*/
       for(i=0;i<n;i++){
          X[i]=0.0;                                    /*Define vetor solução inicial*/
          }

       cont=0;

       do{
          cont= cont+1;
          for(i=0;i<n;i++){                                             /*Método de SOR*/
                           Aux[i]=0.0;
                           for(j=0;j<n;j++){
                                            if(j!=i){
                                            Aux[i] = Aux[i] + A[i][j]*X[j];
                                            }
                           }
                           X[i] = X[i] + w*((B[i]-Aux[i])/A[i][i]-X[i]);
                              //printf("%lf \n", X[i]);
          }
           for(i=0; i<n; i++){                                             /*Obtem os erros*/
                 Aux[i] = 0.0;
                        for(j=0; j<n; j++){
                                Aux[i] = Aux[i] + A[i][j]*X[j];
                            }
                            //printf("%lf \n", Aux[i]);
                                El[i] = fabs(Aux[i]-B[i]);
                            //printf("%lf \n", El[i]);
        }
        j=0;
        for(i=0; i<n; i++){                                                /*Obtem o erro maximo*/
                 max=El[j];
                 if(El[i]>max){
                              max = El[i];
                              j=i;
                              }
        }

         fprintf(GrafS,"%d %lf\n", cont, max);

}while(p<max);


        printf("SOR \n");

        for(i=0; i<n; i++){                  /*Mostra os valores de xi*/

              printf("X%d = %lf \n", i, X[i]);

              }


          fclose(GrafJ);
          fclose(GrafG);
          fclose(GrafS);

           /*Chama o gnuplot com a configuração contida na chamada*/
          plotGraphic("set title 'Metodos iterativos'; set xlabel 'Iteracoes'; set ylabel 'Erro'; set xrange [1:*]; plot 'GrafJ.txt' title 'Metodo de Jacobi' lt 1 w lp,\  'GrafG.txt' title 'Metodo de Gauss-Seidel' lt 2 w lp,\  'GrafS.txt' title 'Metodo SOR' lt 3 w lp; pause 100");



return 0;
}

