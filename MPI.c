#include <stdio.h>
#include <time.h>                   /* clock_t, clock, CLOCKS_PER_SEC */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


//Declaration of the global variables
const int Total_Len= 2000000;
int Num_Proc;
int Proc_Rank;
double Elapsed_Parallel;
double Elapsed_Serial;
clock_t Serial_Start, Serial_Finish, Parallel_Start, Parallel_End;

int main(void) {
    
    //Declaration of Methods
    void ParallelPearsonCC ();
    void SerialPearsonCC();
    
    //Initialization of MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &Num_Proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &Proc_Rank);
    
    if (Proc_Rank == 0) {
        SerialPearsonCC();
    }
    
    if (Num_Proc > 1) {
        ParallelPearsonCC();
    }
    
    MPI_Finalize();
    
}


//Implementation of the Serial Algorithm
void SerialPearsonCC () {
    
    int j; //declares the "j'th" elements of the arrays
    double *a= malloc(Total_Len * sizeof(a));
    double *b= malloc(Total_Len * sizeof(b));
    double Total_Array_a= 0; 
    double Total_Array_b= 0; 
    double Mean_Array_a= 0; 
    double Mean_Array_b= 0; 
    double Stan_Dev_a_sum= 0; 
    double Stan_Dev_a= 0; 
    double Stan_Dev_b_sum= 0; 
    double Stan_Dev_b= 0; 
    double Pear_Sum= 0; 
    double Pear= 0; 
    
    clock_t Serial_Start = clock();
    
    //Initializing Arrays
    for(j= 0; j < Total_Len; j++) {
        a[j]= sin(j);
        b[j]= sin(j +5); 
    }
    
    //Computing sum of array a and b
    for(j= 0; j < Total_Len; j++) {
        Total_Array_a= Total_Array_a + a[j];
        Total_Array_b= Total_Array_b + b[j];
    }
    
    //Computing mean
    Mean_Array_a= Total_Array_a/ Total_Len;
    Mean_Array_b= (Total_Array_b/ Total_Len);
    
    
    //Computing standard deviation of array a and b
    for(j= 0; j < Total_Len; j++) {
        Stan_Dev_a_sum= Stan_Dev_a_sum + (a[j] - Mean_Array_a) * (a[j] - Mean_Array_a);
        Stan_Dev_b_sum= Stan_Dev_b_sum + (b[j] - Mean_Array_b) * (b[j] - Mean_Array_b);
    }
    
    Stan_Dev_a= sqrt(Stan_Dev_a_sum/Total_Len);
    Stan_Dev_b= sqrt(Stan_Dev_b_sum/Total_Len);
    
    //Computing Pearson Correlation between the two arrays
    
    for(j= 0; j < Total_Len; j++) {
        Pear_Sum= Pear_Sum + (a[j] - Mean_Array_a) * (b[j] - Mean_Array_b);
    }
    
    Pear= (Pear_Sum / Total_Len) / (Stan_Dev_a * Stan_Dev_b);
    
    clock_t Serial_Finish = clock();
    Elapsed_Serial = (double)(Serial_Finish - Serial_Start) / CLOCKS_PER_SEC;
    
    printf("\n\n~~~~~~~~~~~     Results of the Serial Execution     ~~~~~~~~~~\n\n");
    printf("\n     mean_a= %lf    mean_b= %lf\n", Mean_Array_a, Mean_Array_b);
    printf("\n     Standard Deviation a= %lf    Standard Deviation b= %lf\n", Stan_Dev_a, Stan_Dev_b);
    printf("\n     Pearson Correlation Coefficient = %lf\n", Pear);
    printf("\n     Execution time is %lf s\n", Elapsed_Serial);
}

// Methods for Parellelization
int Mean(double a[], double b[],int L_N, int First_a, int Last_b,
         double *Temp_Mean_a, double *Temp_Mean_b) {
    
    int i;
    
    double Total_Array_a= 0; //local sum_a
    double Total_Array_b= 0; //local sum_a
    
    
    for(i= First_a; i <= Last_b; i++) {
        a[i]=sin(i);
        b[i]=sin(i+5); 
        Total_Array_a= Total_Array_a + a[i];
        Total_Array_b= Total_Array_b + b[i];
    }
    
    *Temp_Mean_a= Total_Array_a/L_N;
    *Temp_Mean_b = Total_Array_b/L_N;
    
    return 0;
    
}


int Standard_Deviation (double a[], double b[],int local_n, int First_a, int Last_b, double Final_Mean_a, double Final_Mean_b, double *Temp_Stan_Dev_a, double *Temp_Stan_Dev_b) {
    
    int i;
    for(i= First_a; i <= Last_b; i++) {
        *Temp_Stan_Dev_a += (a[i] - Final_Mean_a) * (a[i] - Final_Mean_a);
        *Temp_Stan_Dev_b += (b[i] - Final_Mean_b) * (b[i] - Final_Mean_b);
    }
    
    return 0;
    
}


int Pearson_Correlation_Coefficient (double a[], double b[], int L_N, int First_a, int Last_b, double Final_Mean_a, double Final_Mean_b, double *Temp_Pear) {
    
    int i;
    for(i= First_a; i <= Last_b; i++) {
        *Temp_Pear += (a[i] - Final_Mean_a) * (b[i] - Final_Mean_b);
    }
    
    return 0;
    
}

//Parallel Method

void ParallelPearsonCC () {
    
    double *a= malloc(Total_Len * sizeof(a));
    double *b= malloc(Total_Len * sizeof(b));
    int First_a= 0; //1st element of the subset
    int Last_b= 0; //last element of the subset 
    int L_N= 0; //number of elements in the subset
    double Temp_Mean_a= 0; 
    double Temp_Mean_b= 0; 
    double Final_Mean_a= 0; 
    double Final_Mean_b= 0; 
    double Temp_Stan_Dev_a= 0; 
    double Temp_Stan_Dev_b= 0; 
    double Final_Stan_Dev_a= 0; 
    double Final_Stan_Dev_b= 0; 
    double Temp_Pear= 0; 
    double Pear= 0; 
    
    //Methods used
    int Mean(double a[], double b[],int L_N, int First_a, int Last_b,
             double *Temp_Mean_a, double *Temp_Mean_b);
    
    int Standard_Deviation (double a[], double b[],int L_N, int Temp_Mean_a, int Temp_Mean_b, double Final_Mean_a, double Final_Mean_b, double *Temp_Stan_Dev_a, double *Temp_Stan_Dev_b);
    
    int Pearson_Correlation_Coefficient (double a[], double b[], int L_N, int Temp_Mean_a, int Temp_Mean_b, double Final_Mean_a, double Final_Mean_b, double *Temp_Pear);
    
    clock_t begin_parallel= clock();
    
    //Partioning of Blocks (AUC)
    L_N= Total_Len/Num_Proc;
    
    
    //How to calculate the local subsets of the arrays
    First_a= Proc_Rank * L_N;
    Last_b= First_a+ L_N-1;
    
    
    Mean(a, b, L_N, First_a, Last_b, &Temp_Mean_a, &Temp_Mean_b);
    
    if (Proc_Rank != 0) {
        
        MPI_Send(&Temp_Mean_a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&Temp_Mean_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        
    }
    
    else {
        
        Final_Mean_a += Temp_Mean_a;
        Final_Mean_b += Temp_Mean_b;
        
        
        for (int q = 1; q < Num_Proc; q++) {
            
            
            MPI_Recv(&Temp_Mean_a, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&Temp_Mean_b, 1, MPI_DOUBLE, q, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            Final_Mean_a += Temp_Mean_a;
            Final_Mean_b += Temp_Mean_b;
            
        }
        
        Final_Mean_a= Final_Mean_a/Num_Proc;
        Final_Mean_b= Final_Mean_b/Num_Proc;
        
    }
    
    //Broadcasting the means to all the processes
    
    MPI_Bcast(&Final_Mean_a, 1 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Final_Mean_b, 1 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    Standard_Deviation(a, b, L_N, First_a, Last_b, Final_Mean_a, Final_Mean_b, &Temp_Stan_Dev_a, &Temp_Stan_Dev_b);
   
    
    if (Proc_Rank != 0) {
        
        MPI_Send(&Temp_Stan_Dev_a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&Temp_Stan_Dev_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        
    }
    
    else {
        
        Final_Stan_Dev_a += Temp_Stan_Dev_a;
        Final_Stan_Dev_b += Temp_Stan_Dev_b;
        
        
        for (int q = 1; q < Num_Proc; q++) {
            
            MPI_Recv(&Temp_Stan_Dev_a, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&Temp_Stan_Dev_b, 1, MPI_DOUBLE, q, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Final_Stan_Dev_a += Temp_Stan_Dev_a;
            Final_Stan_Dev_b += Temp_Stan_Dev_b;
            
        }
        
        Final_Stan_Dev_a=sqrt(Final_Stan_Dev_a/Total_Len);
        Final_Stan_Dev_b= sqrt(Final_Stan_Dev_b/Total_Len);
        
    }
    
    Pearson_Correlation_Coefficient(a, b, L_N, First_a, Last_b, Final_Mean_a, Final_Mean_b, &Temp_Pear);
    
    if (Proc_Rank != 0) {
        MPI_Send(&Temp_Pear, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);   
    }
    
    else {
        
        Pear += Temp_Pear;
        for (int q = 1; q < Num_Proc; q++) {
            MPI_Recv(&Temp_Pear, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Pear += Temp_Pear;
        }
        
        clock_t end_parallel = clock();
        Elapsed_Parallel = (double)(end_parallel - begin_parallel) / CLOCKS_PER_SEC;
        
        printf("\n\n~~~~~~~~~~~     Results of the Parallel Execution     ~~~~~~~~~~\n\n");
        printf("\n     Mean a= %lf      Mean b= %lf\n", Final_Mean_a, Final_Mean_b);
        printf("\n     Standard Deviation a= %lf        Standard Deviation b= %lf\n", Final_Stan_Dev_a, Final_Stan_Dev_b);
        Pear= (Pear / Total_Len) / (Final_Stan_Dev_a * Final_Stan_Dev_b);
        printf("\n     Pearson Correltion Coefficient = %lf\n", Pear);
        printf("\n     Execution time: %lf s \n", Elapsed_Parallel);
        printf("\n     Speedup Achieved =  %lf\n", Elapsed_Serial /Elapsed_Parallel);
        printf("\n     Percentage Speedup Achieved = %0.2lf%%\n\n", (1-(Elapsed_Parallel/Elapsed_Serial))*100);
    }
    
}