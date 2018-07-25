#include <stdio.h>
#include <time.h>                   /* clock_t, clock, CLOCKS_PER_SEC */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

//Declaration of the global variables
const int Total_Len= 5000000;
double Elapsed_Parallel;
double Elapsed_Serial;
int j;

//Declaration of Methods
void ParallelPearsonCC ();
void SerialPearsonCC();

// Main method 
int main(void) {
    ParallelPearsonCC ();
    SerialPearsonCC();
    }
//-------------------------------------------------------------------------------------------------------------  
//--------------------------------------------- Serial Algorithm ----------------------------------------------
//-------------------------------------------------------------------------------------------------------------

//Implementation of the Serial Algorithm
void SerialPearsonCC () {
    
    
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
    
    double Begin= 0; //Begin timing
    double Finish= 0; //Finish Timing
    
    Begin= omp_get_wtime();
    
    //Initializing Arrays
    for(j= 0; j < Total_Len; j++) {
        a[j]= sin(j);
        b[j]= sin(j +5); 
        
         //Computing sum of array a and b
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
        
        //Computing Pearson Correlation between the two arrays
        Pear_Sum= Pear_Sum + (a[j] - Mean_Array_a) * (b[j] - Mean_Array_b);
    }
    
    Stan_Dev_a= sqrt(Stan_Dev_a_sum/Total_Len);
    Stan_Dev_b= sqrt(Stan_Dev_b_sum/Total_Len);
    
    
    Pear= (Pear_Sum / Total_Len) / (Stan_Dev_a * Stan_Dev_b);
    
    Finish= omp_get_wtime();
    Elapsed_Serial= (double)(Finish-Begin);
    
    printf("\n\n~~~~~~~~~~~     Results of the Serial Execution     ~~~~~~~~~~\n\n");
    printf("\n     Mean a = %lf    Mean b = %lf\n", Mean_Array_a, Mean_Array_b);
    printf("\n     Standard Deviation a = %lf    Standard Deviation b = %lf\n", Stan_Dev_a, Stan_Dev_b);
    printf("\n     Pearson Correlation Coefficient = %lf\n", Pear);
    printf("\n     Execution time (which includes array initialization) is %lf s\n", Elapsed_Serial);
}
//-------------------------------------------------------------------------------------------------------------  
//--------------------------------------------------- End -----------------------------------------------------
//-------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------  
//-------------------------------------------- Parallel Algorithm ---------------------------------------------
//-------------------------------------------------------------------------------------------------------------

//Implementation of the Parallel Method

void ParallelPearsonCC () {
 	
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
    
    double Begino= 0; //Begin timing
    double Finisho= 0; //Finish Timing
    
    Begino= omp_get_wtime();
    
    #pragma omp parallel for reduction(+:Total_Array_a) reduction(+:Total_Array_b) num_threads(8)
    //Initializing Arrays
    for(j= 0; j < Total_Len; j++) {
        a[j]= sin(j);
        b[j]= sin(j +5); 
        
        //Computing sum of array a and b
        Total_Array_a= Total_Array_a + a[j];
        Total_Array_b= Total_Array_b + b[j];
    }
    
   
    //Computing mean
    #pragma omp parallel sections num_threads(2)
    {
    #pragma omp section {
    Mean_Array_a= Total_Array_a/ Total_Len;
    }
    #pragma omp section {
    Mean_Array_b= (Total_Array_b/ Total_Len);
    }
    
    #pragma omp parallel for reduction(+:Stan_Dev_a_sum) reduction(+:Stan_Dev_b_sum) reduction(+:Pear_Sum) num_threads(8)
    //Computing standard deviation of array a and b
    for(j= 0; j < Total_Len; j++) {
        Stan_Dev_a_sum+= (a[j] - Mean_Array_a) * (a[j] - Mean_Array_a);
        Stan_Dev_b_sum+= (b[j] - Mean_Array_b) * (b[j] - Mean_Array_b);
        
        //Computing Pearson Correlation between the two arrays
        Pear_Sum+= (a[j] - Mean_Array_a) * (b[j] - Mean_Array_b);
    }
    #pragma omp parallel sections num_threads(2){

	#pragma omp section{
    Stan_Dev_a= sqrt(Stan_Dev_a_sum/Total_Len);
    }
    #pragma omp section{
    Stan_Dev_b= sqrt(Stan_Dev_b_sum/Total_Len);
    }
    
    }
    
    Pear= (Pear_Sum / Total_Len) / (Stan_Dev_a * Stan_Dev_b);
    
    Finisho= omp_get_wtime();
    Elapsed_Serial= (double)(Finisho-Begino);
   
    printf("\n\n~~~~~~~~~~~     Results of the Parallel Execution     ~~~~~~~~~~\n\n");
    printf("\n     Mean a = %lf    Mean b = %lf\n", Mean_Array_a, Mean_Array_b);
    printf("\n     Standard Deviation a = %lf    Standard Deviation b = %lf\n", Stan_Dev_a, Stan_Dev_b);
    printf("\n     Pearson Correlation Coefficient = %lf\n", Pear);
    printf("\n     Execution time (which includes array initialization) is %lf s\n", Elapsed_Serial);
}
//-------------------------------------------------------------------------------------------------------------  
//--------------------------------------------------- End -----------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
    
    