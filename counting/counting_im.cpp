//codes corrected by Piotr Blaszynski

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

int ** c;
int ** ck;

int zz = 2;

int ** F;  //only ACGU

int N;
int DIM;

char * RNA;  //only ACGU


#include "mem.h"

int pared(int i, int j) {   
   char nt1 = RNA[i];
   char nt2 = RNA[j];
         if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') ||
             (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') ||
             (nt1 == 'G' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'G')){

            return 1;}
         else
            return 0;
}



int main(int argc, char *argv[]){



    int num_proc=1;
    int i,j,k,ll,p,q,l=0;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6,t7,t8,t9,t10;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;




    srand(time(NULL));



    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    N = 8;
    DIM = 12;
    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+10;


    if(argc > 3)
        kind = atoi(argv[3]);



    omp_set_num_threads(num_proc);
    //printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

    F =  mem();
    c = mem();
    ck = mem();

   for(i=0; i<DIM; i++)
    for(i=0; i<DIM; i++){
     c[i][j] = i+j;
     ck[i][j] = i+j;
    }

    RNA =  (char*) malloc(DIM * sizeof(char*));  //read from FASTA file
    rand_seq(RNA, N);

for(i=0; i<DIM; i++)
 printf("%c", RNA[i]);

printf("\n");
     int check=1;


    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if(kind==1 || check){
        #pragma scop
    for (int i = N - 2; i >= 1; i--)
    {
      for (int j = i + 2; j <= N; j++)
      {
        for (int k = i; k <= j - 1; k++)
        {
          ck[i][j] += paired(k, j) ? ck[i][k - 1] + ck[k + 1][j - 1] : 0;
        }
        ck[i][j] = ck[i][j] + ck[i][j - 1];
      }
    }
       #pragma endscop
    }
    if(kind==2) // pluto
    {
      // pluto
        if (N >= 3)
    {
      for (t1 = 3; t1 <= N; t1++)
      {
        lbp = 0;
        ubp = floord(t1 - 2, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5)
        for (t2 = lbp; t2 <= ubp; t2++)
        {
          for (t3 = t2; t3 <= floord(t1, 32); t3++)
          {
            if ((t1 >= 32 * t3 + 1) && (t1 <= 32 * t3 + 31))
            {
              for (t4 = max(1, 32 * t2); t4 <= min(t1 - 2, 32 * t2 + 31); t4++)
              {
                for (t5 = max(32 * t3, t4); t5 <= t1 - 1; t5++)
                {
                  c[t4][t1] += pared(t5, t1) ? c[t4][t5 - 1] + c[t5 + 1][t1 - 1] : 0;
                  ;
                }
                c[t4][t1] = c[t4][t1] + c[t4][t1 - 1];
                ;
              }
            }
            if (t1 >= 32 * t3 + 32)
            {
              for (t4 = max(1, 32 * t2); t4 <= min(t1 - 2, 32 * t2 + 31); t4++)
              {
                for (t5 = max(32 * t3, t4); t5 <= 32 * t3 + 31; t5++)
                {
                  c[t4][t1] += pared(t5, t1) ? c[t4][t5 - 1] + c[t5 + 1][t1 - 1] : 0;
                  ;
                }
              }
            }
            if (t1 == 32 * t3)
            {
              for (t4 = max(1, 32 * t2); t4 <= min(t1 - 2, 32 * t2 + 31); t4++)
              {
                if (t1 % 32 == 0)
                {
                  c[t4][t1] = c[t4][t1] + c[t4][t1 - 1];
                  ;
                }
              }
            }
          }
        }
      }
    }
    // end pluto
    }
    if(kind==3) // traco
    {
       
 for (c1 = 0; c1 < N + floord(N - 3, 128) - 2; c1 += 1)
#pragma omp parallel for
      for (c3 = max(0, -N + c1 + 3); c3 <= c1 / 129; c3 += 1)
        for (c4 = 0; c4 <= 1; c4 += 1)
        {
          if (c4 == 1)
          {
            for (c9 = N - c1 + 129 * c3; c9 <= min(N, N - c1 + 129 * c3 + 127); c9 += 1)
              for (c10 = max(0, -c1 + 64 * c3 - c9 + (N + c1 + c3 + c9 + 1) / 2 + 1); c10 <= 1; c10 += 1)
              {
                if (c10 == 1)
                {
                  c[(N - c1 + c3 - 2)][c9] = c[(N - c1 + c3 - 2)][c9] + c[(N - c1 + c3 - 2)][c9 - 1];
                }
                else
                {
                  for (c11 = N - c1 + 129 * c3 + 1; c11 < c9; c11 += 1)
                    c[(N - c1 + c3 - 2)][c9] += pared(c11, c9) ? c[(N - c1 + c3 - 2)][c11 - 1] + c[c11 + 1][c9 - 1] : 0;
                }
              }
          }
          else
          {
            for (c5 = 0; c5 <= 8 * c3; c5 += 1)
              for (c9 = N - c1 + 129 * c3; c9 <= min(N, N - c1 + 129 * c3 + 127); c9 += 1)
                for (c11 = N - c1 + c3 + 16 * c5 - 2; c11 <= min(min(N - c1 + 129 * c3, N - c1 + c3 + 16 * c5 + 13), c9 - 1); c11 += 1)
                  c[(N - c1 + c3 - 2)][c9] += paired(c11, c9) + c[(N - c1 + c3 - 2)][c11 - 1] + c[c11 + 1][c9 - 1] + 0;
          }
        }
       // end
    }
   
   if(kind==5) // dapt
   {
       for (int w0 = floord(-N - 14, 32); w0 < floord(N, 32); w0 += 1)
    {
#pragma omp parallel for
      for (int h0 = max(-((N + 13) / 16), w0 - (N + 32) / 32 + 1); h0 <= min(-1, 2 * w0 + 2); h0 += 1)
      {
        for (int i0 = max(max(-N + 2, -32 * w0 + 32 * h0 - 29), 16 * h0); i0 <= 16 * h0 + 15; i0 += 1)
        {
          for (int i1 = max(32 * w0 - 32 * h0, -i0 + 2); i1 <= min(N, 32 * w0 - 32 * h0 + 31); i1 += 1)
          {
            {
              for (int i3 = -i0; i3 < i1; i3 += 1)
              {
                c[-i0][i1] += (paired((i3), (i1)) ? (c[-i0][i3 - 1] + c[i3 + 1][i1 - 1]) : 0);
              }
              c[-i0][i1] = (c[-i0][i1] + c[-i0][i1 - 1]);
            }
          }
        }
      }
    }
      
      
   }
 






    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

    //printf("Q\n");
    //rna_array_print(Q);
    //printf("Qbp\n");
    //rna_array_print(Qbp);

    for(i=0; i<DIM; i++)
    for(j=0; j<DIM; j++)
     if(c[i][j] != ck[i][j]){
        printf("err: %d %d %d %d\n", i, j,c[i][j], ck[i][j]);
        exit(0);
     }


    return 0;

}

