#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#include<string.h>

#define n 63



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

int preencheVetor(double *v, int k){
    
    int j;
    double pi;
                         
    pi = 2.0*asin(1);
                         
    for(j=0; j<=n+1; j++){
                          v[j] = sin((k*j*pi)/(n+1));
                          }                  
    return 0;                         
}     

void plotGraphic(const char *gnucommand){
    char syscommand[1024];

    sprintf(syscommand, "echo \%s\ | gnuplot -persist", gnucommand);
    system(syscommand);
}




int main(void){
     double **A, *V, *B, *Aux, *E, somaE, Er, max;
     int i, j, k, K, cont;
     FILE  *GrafG1, *GrafG2, *GrafG3;
     
         
  
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
  
  GrafG3 = fopen("GrafG3.txt", "w");
  if (GrafG3 == NULL){
    printf("erro na apertura do arquivo GrafG3.txt\n");
    exit(EXIT_FAILURE);
  }
  

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

  V = malloc(n * sizeof(double));
  if(V == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
  
  B = malloc(n * sizeof(double));
  if(B == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  Aux = malloc(n * sizeof(double));
  if(Aux == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }

  E = malloc(n * sizeof(double));
  if(E == NULL){
    printf("Erro de alocacao de memoria\n");
    exit(1);
  }
preencheMatriz(A);  

for(K=0; K<3; K++){
           
printf("Escolha a frequencia: ");
scanf("%d", &k);

preencheVetor(V, k);

/*-----------------------------------------------------Método de Gauss-Seidel--------------------------------------------------------*/
      for(i=0;i<n;i++){
                       B[i]=0.0;                                    /*Define vetor solução inicial*/
       }
       cont=0;
       Er=0.0;
       
       do{
          cont= cont+1;         
          for(i=0;i<n;i++){                                             /*Metodo de Gauss-Seidel*/
                           Aux[i]=0.0;
                           for(j=0;j<n;j++){
                                            if(i!=j){
                                                     Aux[i] = Aux[i]+ A[i][j]*V[j];
                                                     }
                                            }
                           V[i]=(B[i]-Aux[i])/A[i][i];
                           E[i] = fabs(-V[i]);                                  /* Toma o erro igual a V, mas no livro deveria ser E=-V */
                          
          }
        
          j=0;
          for(i=0; i<n; i++){
                             max=E[j];
                             if(E[i]>max){
                                          max = E[i];
                                          j=i;
                             }
          }
      
                                                  /*grava dados no arquivo para fazer o grafico*/
        
        if(K==0){
        fprintf(GrafG1,"%d %lf\n", cont, max);
        } 
        if(K==1){
        fprintf(GrafG2,"%d %lf\n", cont, max);
        }
        if(K==2){
        fprintf(GrafG3,"%d %lf\n", cont, max);
        }                    
             
           
        }while(cont<101);
}     
       

              
         
          fclose(GrafG1);    
          fclose(GrafG2);
          fclose(GrafG3);
          
          
           /*Chama o gnuplot com a configuração contida na chamada*/    
          plotGraphic("set title 'Gauss-Seidel em Grid'; set xlabel 'Iteracoes'; set ylabel 'Erro'; set xrange [1:100]; set yrange [0:1]; plot 'GrafG1.txt' title 'k=1' lt 1 w lp,\ 'GrafG2.txt' title 'k=3' lt 2 w lp,\ 'GrafG3.txt' title 'k=6' lt 3 w lp; pause 100"); 

return 0;              
}             

