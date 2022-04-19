#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<string.h>

#include "common.h"
#include "jacobi.h"


#define jacobi_eps 1.0e-15
#define jacobi_max_iter 100

double sign(double x);
void calc_c_and_t(double **A, int i, int j, double *c, double *s);
void find_Pivot(int n,double **A, int *i, int *j);
void reset_mat(int n,double **mat);
double calc_off2(int n,double **mat);
void copy_mat(int n, double** src,double **dest);
void Build_P_and_A_tag(int n,double **A, double **A_tag, double **P,double *c, double *s);
double **duplicate_mat(int n, double** mat);
void mult_mats(int n,double **X, double **Y, double **dest);
void swap_mats(double ***X,double ***Y);
void init_Imat(int n,double **X);
void free_memory_of_matrix(double** matrix,int number_of_rows);
double** Jacobi_algorithm(int n,double **A, int *res_rows, int *res_cols);


double sign(double x){
    if (x >= 0)
        return 1;
    else
        return -1;
}


void calc_c_and_t(double **A, int i, int j, double *c, double *s){

    double teta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    double t = sign(teta)/(fabs(teta) + sqrt(pow(teta,2) + 1));
    *c = (1)/(sqrt(pow(t,2) + 1));
    *s = (*c)*t;

}

void find_Pivot(int n,double **A, int *i, int *j){
    int row,col;
    double max_val;
    *i = 0, *j = 1;
    /*Assert n > 1*/
    max_val = fabs(A[0][1]);
    for ( row = 0; row < n; ++row) {
        for ( col = row + 1; col < n; ++col) {
            if(fabs(A[row][col]) > max_val){
                max_val = fabs(A[row][col]);
                *i = row, *j = col;
            }
        }
    }
}

void reset_mat(int n,double **mat){
    int i,j;
    for (i = 0;  i < n ; i++) {
        for (j = 0; j < n; ++j) {
            mat[i][j] = 0;
        }
    }
}

double calc_off2(int n,double **mat){
    int i,j;
    double sum = 0;
    for (i = 0; i < n; ++i) {
        for (j = i+1; j < n; ++j) {
            sum += pow(mat[i][j],2);
        }
    }

    sum *=2;

    return sum;
}

void copy_mat(int n, double** src,double **dest){
    int row,col;
    for (row = 0; row < n; ++row) {
        for (col = 0; col < n; ++col) {
            dest[row][col] = src[row][col];
        }
    }
}


void Build_P_and_A_tag(int n,double **A, double **A_tag, double **P,double *c, double *s){
    /*Assumption - p is all zero matrix*/
    int i,j,k;
    /*Build P*/
    find_Pivot(n,A, &i, &j);
    calc_c_and_t(A , i, j, c, s);

    init_Imat(n,P);

    P[i][i] = *c;
    P[i][j] = *s;
    P[j][i] = -*s;
    P[j][j] = *c;


    copy_mat(n,A,A_tag);
    for (k = 0; k < n; ++k) {
            if (k != i && k != j){
                A_tag[k][i] = (*c) * A[k][i] - (*s) * A[k][j];
                A_tag[i][k] = A_tag[k][i];
                A_tag[k][j] = (*c) * A[k][j] + (*s) * A[k][i];
                A_tag[j][k] = A_tag[k][j];
            }
        }

    A_tag[i][i] = ((*c)*(*c)) * A[i][i] + ((*s) * (*s)) * A[j][j] - 2 * (*s) * (*c) * A[i][j];
    A_tag[j][j] = ((*s)*(*s)) * A[i][i] + ((*c) * (*c)) * A[j][j] + 2 * (*s) * (*c) * A[i][j];
    A_tag[i][j] = 0;
    A_tag[j][i] = 0;

}

double ** duplicate_mat(int n, double** mat){
    int row,col;
    double **new_mat = (double**)calloc(n,sizeof(double*));
    /* add assert*/

    for (row = 0; row < n; row++) {
        new_mat[row] = (double *) calloc(n, sizeof(double));
        /*Add assert*/
    }
    for (row = 0; row < n; ++row) {
        for (col = 0; col < n; ++col) {
            new_mat[row][col] = mat[row][col];
        }
    }
    return new_mat;
}


void mult_mats(int n,double **X, double **Y, double **dest){
    double sum;
    int i,j,k;
    for ( i = 0; i < n; ++i) {
        for ( j = 0; j < n; ++j) {
            sum = 0;
            for ( k = 0; k < n; ++k) {
                sum += (X[i][k])*(Y[k][j]);
            }
            dest[i][j] = sum;
        }
    }
}

void swap_mats(double ***X,double ***Y){
    double **tmp;
    tmp = *X;
    *X = *Y;
    *Y = tmp;
}

void init_Imat(int n,double **X){
    int i;
    reset_mat(n,X);
    for (i = 0; i < n; ++i) {
        X[i][i] = 1;
    }
}

void free_memory_of_matrix(double** matrix,int number_of_rows){
    int i;
    for(i=0; i<number_of_rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}

double** Jacobi_algorithm(int n,double **A, int *res_rows, int *res_cols){
    int row, col,curr_iter = 0;
    double c,s;
    /* just as first value to enter main loop*/
    double A_off = (double)(2*jacobi_eps),A_tag_off = 0;

    double **A_tag = (double**)calloc(n,sizeof(double*));
    double **P = (double**)calloc(n,sizeof(double*));
    double **V = (double**)calloc(n,sizeof(double*));
    double **tmp_V = (double**)calloc(n,sizeof(double*));
    double **jacobi_res = (double**)calloc(n,sizeof(double*));
    double **A_clone = duplicate_mat(n,A);

    /* Add assert*/
    for (row = 0; row < n; row++) {
        P[row] = (double*) calloc(n,sizeof(double));
        A_tag[row] = (double*) calloc(n,sizeof(double));
        V[row] = (double*) calloc(n,sizeof(double));
        tmp_V[row] = (double*) calloc(n,sizeof(double));
        jacobi_res[row] = (double*) calloc(n+1,sizeof(double));
        /* Add assert*/
    }

    init_Imat(n,V);

    while(curr_iter < jacobi_max_iter && A_off-A_tag_off > jacobi_eps){

        Build_P_and_A_tag(n,A_clone,A_tag, P, &c, &s);
        A_off = calc_off2(n,A_clone);
        A_tag_off = calc_off2(n,A_tag);
        curr_iter++;
        mult_mats(n,V,P,tmp_V);
        swap_mats(&V,&tmp_V);
        swap_mats(&A_clone,&A_tag);


    }

    /*Build jacobi_res*/
    
    for (row = 0; row < n; ++row) {
        jacobi_res[row][0] = A_clone[row][row];
        for ( col = 0; col < n; ++col) {
            jacobi_res[row][col+1] = V[col][row];
        }
    
    }
    *res_rows =n;
    *res_cols = n+1;

    free_memory_of_matrix(P,n);
    free_memory_of_matrix(V,n);
    free_memory_of_matrix(tmp_V,n);
    free_memory_of_matrix(A_clone,n);
    free_memory_of_matrix(A_tag,n);

    return jacobi_res;


}


int main() {

    int resrow,rescol;
    double **mat,**res;
    int number_of_cords,number_of_lines;
    FILE *fp;

     
    fp = fopen("lnorm_input_to_jacobi_ramis_input.txt","r");
    number_of_cords = compute_number_of_cord(fp);
    number_of_lines = compute_number_of_x(fp);    
    mat = read_data_from_file(fp,number_of_cords,number_of_lines);
    fclose(fp);
    printf("input matrix: \n");
    print_matrix(mat,number_of_lines,number_of_cords);

    res = Jacobi_algorithm(number_of_lines,mat,&resrow,&rescol);
    res = transpose(res,resrow,rescol);
    printf("jacobi matrix: \n");
    print_matrix(res,rescol,resrow);

    return 0;
}