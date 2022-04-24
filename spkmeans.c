#include "spkmeansmodule.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "spkmeans.h"
#include "common.h"
#include "Kmeans.h"

#include <stdio.h>



#define jacobi_eps pow(10,-5)
#define jacobi_max_iter 100

typedef enum {spk,wam,ddg,lnorm,jacobi} goals;

/* function declaration*/
int initial_rows(double** W, int number_of_lines,int number_of_columns);
double compute_wij(double* xi, double* xj, int number_of_cords);
double** weighted_adjacency_matrix(double** X, int number_of_lines,int number_of_cord);
double degree_of_line_i(double* line, int number_of_cords);
void raise_D_minus_half(double* D, int number_of_lines);
double* diagonal_degree_matrix(double** W, int number_of_lines);
double** subtract_M_from_I(double** M, int number_of_lines);
void array_matrix_representation(double** mult, double**W, double* D, int row, int column);
void matrix_array_representation(double** mult, double**M, double* D, int row, int column);
double** mult_matrices(double** M, double* D, int number_of_lines, void representation(double**, double**, double*, int, int));
int find_k(double** V, int number_of_lines);
double** transpose(double** matrix, int number_of_lines, int number_of_columns);
double** create_eigenvectors_matrix(double** V, int number_of_cords, int number_of_lines);
double sum_of_row(double* row, int number_of_column);
double** normalize_U(double** U, int number_of_lines, int number_of_columns);
void print_matrix(double** matrix, int number_of_rows,int number_of_columns,char* name_of_matrix);
void printmatrix(double** matrix, int number_of_rows,int number_of_columns);
double** calculate_Lnorm(double* D, double**W, int number_of_lines);
int submit_args_spkmeans(int argc, char **argv,char **input_file, int *goal);
double** Spkmeans(int* rows,int*cols, char* input_filename, int purpose, int k);
int get_enum_val(char* str);
void fix_zeros(double** V, int number_of_columns);



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ jacobi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*perform matrix multiplication  x*y=dest*/
void mult_mats_for_jacobi(int n,double **X, double **Y, double **dest){
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

/* returns 1 if x>=0 else -1 */
double sign(double x){
    if (x >= 0)
        return 1;
    else
        return -1;
}

/* calculate c and t values for jacobi algorithm */
void calc_c_and_t(double **A, int i, int j, double *c, double *s){

    double theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    double t = sign(theta)/(fabs(theta) + sqrt(pow(theta,2) + 1));
    *c = (1)/(sqrt(pow(t,2) + 1));
    *s = (*c)*t;
}

/* finds pivot for jacobi algorithm, asserts n>1 */
int find_Pivot(int n,double **A, int *i, int *j){
    int row,col;
    double max_val = 0;
    *i = 0, *j = 1;

    for ( row = 0; row < n; ++row) {
        for ( col = row + 1; col < n; ++col) {

            if(fabs(A[row][col]) > max_val){
                max_val = fabs(A[row][col]);
                *i = row, *j = col;
            }
        }

    }

    if(max_val ==0)
        return -1;
    else
        return 0;

}


/* reset mat to 0 matrix */
void reset_mat(int n,double **mat){
    int i,j;
    for (i = 0;  i < n ; i++) {
        for (j = 0; j < n; ++j) {
            mat[i][j] = 0;
        }
    }
}

/* calculate offset squared for jacobi algorithm */
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

/* copy src -> dest */
void copy_mat(int n, double** src,double **dest){
    int row,col;
    for (row = 0; row < n; ++row) {
        for (col = 0; col < n; ++col) {
            dest[row][col] = src[row][col];
        }
    }
}

/* init X as I mat */
void init_Imat(int n,double **X){
    int i;
    reset_mat(n,X);
    for ( i = 0; i < n; ++i) {
        X[i][i] = 1;
    }
}

/* Build P and A_tag matrices for jacobi algorithm
 * return 0 on success else -1
 * Assert n > 1
 * */
int Build_P_and_A_tag(int n,double **A, double **A_tag, double **P){

    int i,j,k;
    double c,s;
    /*Build P*/
    if(find_Pivot(n,A, &i, &j) == -1)
        return -1;
    calc_c_and_t(A , i, j, &c, &s);
    init_Imat(n,P);
    P[i][i] = c;
    P[i][j] = s;
    P[j][i] = -s;
    P[j][j] = c;

    /*Build A'*/
    copy_mat(n,A,A_tag);
    for (k = 0; k < n; ++k) {
                A_tag[k][i] = (c * A[k][i]) - (s * A[k][j]);
                A_tag[i][k] = A_tag[k][i];
                A_tag[k][j] = (c * A[k][j]) + (s * A[k][i]);
                A_tag[j][k] = A_tag[k][j];
        }
    A_tag[i][i] = (pow(c,2) * A[i][i]) + (pow(s,2) * A[j][j]) - (2 * s * c * A[i][j]);
    A_tag[j][j] = (pow(s,2) * A[i][i]) + (pow(c,2) * A[j][j]) + (2 * s * c * A[i][j]);
    A_tag[i][j] = 0;
    A_tag[j][i] = 0;

    return 0;
}


/* Allocate a new matrix and copy mat to it
 * returns NULL on failing*/
double ** duplicate_mat(int n, double** mat){
    int row;
    double **new_mat = (double**)calloc(n,sizeof(double*));

    if(new_mat == NULL)
        return NULL;

    for (row = 0; row < n; row++) {
        new_mat[row] = (double *) calloc(n, sizeof(double));

        if(new_mat[row] == NULL)
            return NULL;
    }

    copy_mat(n,mat,new_mat);
    return new_mat;
}

/* preform A swap, expect pointers to matrices */
void swap_mats(double ***X,double ***Y){
    double **tmp;
    tmp = *X;
    *X = *Y;
    *Y = tmp;
}

/* free memory of matrix */
void free_memory_of_matrix(double** matrix,int number_of_rows){
    int i;
    for(i=0; i<number_of_rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}

/* Calculate jacobi eigenvalues and eigenvectors
 * returns matrix (n+1 X n) like:
 *
 *  e.value(1) | e.vector(1,1) ... e.vector(1,n)
 *                  .
 *                  .
 *                  .
 *  e.value(n) | e.vector(n,1) ... e.vector(n,n)
 */
double** Jacobi_algorithm(int* res_row, int*res_col, int n,double **A){

    int row,col, curr_iter = 0;

    /*st as first value to enter main loop*/
    double A_off = (double)(2*jacobi_eps),A_tag_off = 0;

    double **A_tag = (double**)calloc(n,sizeof(double*));
    double **P = (double**)calloc(n,sizeof(double*));
    double **V = (double**)calloc(n,sizeof(double*));
    double **tmp_V = (double**)calloc(n,sizeof(double*));
    double **jacobi_res = (double**)calloc(n,sizeof(double*));
    double **A_clone = duplicate_mat(n,A);

    if(A_tag==NULL || P==NULL || V == NULL || tmp_V == NULL || jacobi_res == NULL || A_clone == NULL)
        return NULL;

    for (row = 0; row < n; row++) {
        P[row] = (double*) calloc(n,sizeof(double));
        A_tag[row] = (double*) calloc(n,sizeof(double));
        V[row] = (double*) calloc(n,sizeof(double));
        tmp_V[row] = (double*) calloc(n,sizeof(double));
        jacobi_res[row] = (double*) calloc(n+1,sizeof(double));

        if(A_tag[row] == NULL || P[row] == NULL || V[row] == NULL || tmp_V[row] == NULL || jacobi_res[row] == NULL)
            return NULL;

    }

    init_Imat(n,V);

    while(curr_iter < jacobi_max_iter && (A_off-A_tag_off) > jacobi_eps){

        A_off = calc_off2(n,A_clone);
        Build_P_and_A_tag(n,A_clone,A_tag, P);

        A_tag_off = calc_off2(n,A_tag);
        curr_iter++;
        mult_mats_for_jacobi(n,V,P,tmp_V);

        /* Swap matrices as a tool to save matrix allocations */
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

    *res_row = n;
    *res_col = n+1;

    free_memory_of_matrix(P,n);
    free_memory_of_matrix(V,n);
    free_memory_of_matrix(tmp_V,n);
    free_memory_of_matrix(A_clone,n);
    free_memory_of_matrix(A_tag,n);
    free_memory_of_matrix(A,n);

    return jacobi_res;


}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end jacobi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ spkmeans ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*assign memory to each line in matrix W
if succeeded returns 1, else 0*/
int initial_rows(double** W, int number_of_lines, int number_of_column){
    int row;
    for(row=0; row<number_of_lines;row++){
        W[row] = (double*)calloc(number_of_column,sizeof(double));
        if(W[row] == NULL) return 0;
    }
    return 1;
}

/*calculate element in the weighted adjacency matrix*/
double compute_wij(double* xi, double* xj, int number_of_cords){
    double power= -euqlide_norm(xi,xj,number_of_cords)/2;
    return exp(power);
}

/*creates the weighted adjacency matrix
 if didn't succeed returns NULL*/
double** weighted_adjacency_matrix(double** X, int number_of_lines,int number_of_cord){
    int row,column,s;
    double** W = (double**)calloc(number_of_lines,sizeof(double*));
    if(W == NULL) return NULL;
    s = initial_rows(W,number_of_lines, number_of_lines);
    if(s == 0) return NULL;
    for(row =0;row<number_of_lines;row++){
        for(column=row;column<number_of_lines;column++){
            if(row == column){
                W[row][column] = 0;
            }else{
                W[row][column] = compute_wij(X[row],X[column], number_of_cord);
                W[column][row] = W[row][column]; /*W is symmetric*/
            }
        }
    }
    free_memory_of_matrix(X,number_of_lines);
    return W;
}

/*compute the degree of a vertex xi -> sum of the line i*/
double degree_of_line_i(double* line, int number_of_cords){
    int col;
    double sum=0;
    for(col =0; col< number_of_cords; col++){
        sum += line[col];
    }
    return sum;
}

/*compute D^(-0.5) -> for every element d in the main diagonal of D calculate d^(-0.5)*/
void raise_D_minus_half(double* D, int number_of_lines){
    int row;
    for(row =0; row<number_of_lines;row++){
        /*D[row] = pow(D[row], -0.5);*/
        D[row] = 1/sqrt(D[row]);
    }
}

/*creates the diagonal degree matrix
if didn't succeed returns NULL*/
double* diagonal_degree_matrix(double** W, int number_of_lines){
    int row;
    double* D = (double*)calloc(number_of_lines, sizeof(double));
    if(D == NULL) return NULL;
    for(row =0; row<number_of_lines; row++){
        D[row] = degree_of_line_i(W[row], number_of_lines);
    }
    return D;
}

/*subtract a square matrix M from the identity matrix
if didn't succeed returns NULL*/
double** subtract_M_from_I(double** M, int number_of_lines){
    int row, column,s;
    double** sub = (double**)calloc(number_of_lines,sizeof(double*));
    if(sub == NULL) return NULL;
    s = initial_rows(sub,number_of_lines,number_of_lines);
    if(s == 0) return NULL;
    for(row =0; row<number_of_lines; row++){
        for(column =0; column<number_of_lines;column++){
            if(column==row){
                sub[row][column] = 1 - M[row][column];
            }else{
                sub[row][column] = -M[row][column];
            }
        }
    }
    free_memory_of_matrix(M,number_of_lines);
    return sub;
}

/*calculate element in the case of array*matrix */
void array_matrix_representation(double** mult, double**W, double* D, int row, int column){
    mult[row][column] = D[row] * W[row][column];
}

/*calculate element in the case of matrix*array */
void matrix_array_representation(double** mult, double**M, double* D, int row, int column){
    mult[row][column]= M[row][column]*D[column];
}

/*calculate M*D or D*M
if didn't succeed returns NULL*/
double** mult_matrices(double** M, double* D, int number_of_lines, void representation(double**, double**, double*, int, int)){
    int row, column,s;
    double** mult = (double**)calloc(number_of_lines,sizeof(double*));
    if(mult == NULL) return NULL;
    s= initial_rows(mult,number_of_lines,number_of_lines);
    if(s == 0) return NULL;
    for(row = 0; row<number_of_lines;row++){
        for(column =0; column<number_of_lines; column++){
            representation(mult,M,D,row,column);
        }  
    }
    return mult;
}


int comparator(const void* row1, const void* row2){
    double* line1 = *(double* const*)row1;
    double* line2 = *(double* const*)row2;
    if(line1[0]>line2[0]) return 1;
    if(line1[0]<line2[0]) return -1;
    return 0;
}

/*finding K in the heuristic way - find the max eigen gap and returns number of eigenvalues before that gap*/
int find_k(double** V, int number_of_lines){
    int curr,k =0;
    double max_delta =0,delta;
    for(curr=0; curr<(int)(number_of_lines/2); curr++){
        delta = fabs(V[curr][0]- V[curr+1][0]);
        if(max_delta < delta){
            max_delta = delta;
            k = curr;
        }
    }
    return k+1;
}

/*transpose a matrix
if didn't succeed returns NULL*/
double** transpose(double** matrix, int number_of_lines, int number_of_columns){
    int row, column,s;
    double** T = (double**)calloc(number_of_columns,sizeof(double*));
    if(T == NULL) return NULL;
    s = initial_rows(T,number_of_columns,number_of_lines);
    if(s == 0) return NULL;
    for(row =0; row<number_of_lines;row++){
        for(column=0; column<number_of_columns;column++){
            T[column][row] = matrix[row][column];
        }
    }
    free_memory_of_matrix(matrix, number_of_lines);
    return T;
}

/*create eigenvectors matrix from the matrix V. 
in each row of V - first element is eigenvalue, next is its eigenvector
and also does the transpose
if didn't succeed returns NULL*/
double** create_eigenvectors_matrix(double** V, int number_of_cords, int number_of_lines){
    int row,columns,s;
    double** U = (double**)calloc(number_of_lines,sizeof(double*));
    if(U == NULL) return NULL;
    s = initial_rows(U,number_of_lines,number_of_cords-1);
    if(s == 0) return NULL;

    for(row=0;row<number_of_lines;row++){
        for(columns = 1; columns<number_of_cords; columns++){
            U[row][columns-1] = V[row][columns];
        }
    }
    free_memory_of_matrix(V,number_of_lines);
    U = transpose(U,number_of_lines,number_of_cords-1);
    return U;
}

/*calculate the sum og the element squares in a row*/
double sum_of_row(double* row, int number_of_column){
    double sum=0;
    int cord;
    for(cord=0; cord<number_of_column; cord++){
        sum+=pow(row[cord],2);
    }
    return sum;
}

/*normalizes a matrix - sum of each row = 1
if didn't succeed returns NULL*/
double** normalize_U(double** U, int number_of_lines, int number_of_columns){
    int row,column,s;
    double denominator;
    double** T = (double**)calloc(number_of_lines,sizeof(double*));
    if(T == NULL) return NULL;
    s= initial_rows(T,number_of_lines,number_of_columns);
    if(s == 0) return NULL;
    for(row=0;row<number_of_lines;row++){
        denominator = sqrt(sum_of_row(U[row],number_of_columns));
        for(column=0;column<number_of_columns;column++){
            if(denominator ==0){
                T[row][column] = denominator;
            }else{
                T[row][column] = U[row][column]/denominator;
            }
        }
    }
    free_memory_of_matrix(U,number_of_lines);
    return T;
}


/*Auxiliary function that prints matrices*/
void printmatrix(double** matrix, int number_of_rows,int number_of_columns){
    int row,column;
    for(row=0;row<number_of_rows;row++){
        for(column=0;column<number_of_columns;column++){
            if(column!=number_of_columns-1){
                printf("%.4f,",matrix[row][column]);
            }else{
                printf("%.4f\n",matrix[row][column]);
            }
        }
    }
}

/*calculate Lnorm = I - D^(-0.5)WD^(-0.5)
if didn't succeed returns NULL */
double** calculate_Lnorm(double* D, double**W, int number_of_lines){
   double** temp1, **temp2;
    raise_D_minus_half(D,number_of_lines);
    temp1 = mult_matrices(W,D,number_of_lines,array_matrix_representation);
    if(temp1 == NULL) return NULL;
    temp2 = mult_matrices(temp1,D,number_of_lines, matrix_array_representation);
    if(temp2 == NULL) return NULL;
    free_memory_of_matrix(temp1,number_of_lines);
    return subtract_M_from_I(temp2,number_of_lines);
}

/*Turns an array that represents a diagonal of a matrix into a matrix
if didn't succeed returns NULL*/
double** D_asMatrix(double* D, int rows){
    int row,col,s;
    double **D_matrix = (double**)calloc(rows,sizeof(double*));
    if(D_matrix == NULL) return NULL;
    s= initial_rows(D_matrix,rows,rows);
    if(s == 0) return NULL;
    for(row =0;row<rows;row++){
        for(col=0; col<rows;col++){
            if(row == col)
                D_matrix[row][col] = D[row];
            else
                D_matrix[row][col] = 0;
        }
    }
    free(D);
    return D_matrix;
}


int check_input_symmetric_matrix(double** input,int number_of_lines, int number_of_cords){
    int row,col;
    if(number_of_cords != number_of_lines){
        return -1;
    }

    for(row =0; row<number_of_lines;row++){
        for(col =0; col <number_of_cords;col++){
            if(input[row][col] != input[col][row])
                return -1;
        }
    }
    return 0;
}
/*reset to zero each eigen value has a value bigger than -0.00005*/
void fix_zeros(double** V, int number_of_columns){
    int i;
      for (i = 0; i < number_of_columns; ++i) {
        if(V[0][i] < 0 && V[0][i] > -0.00005) V[0][i] = 0;    
    } 
}



/*gets data from file, goal and k -returns the appropriate matrix to the given data
if didn't succeed returns NULL*/
double** Spkmeans (int* row, int*col, char* input_filename, int purpose, int k_from_py){
    int number_of_cords, number_of_lines,k,res_row,res_col,symmetric;
    double **X, **W, **Lnorm, **eigen_values_vectors, **k_eigenvectors,**T,**Jacobi;
    double *D ;
    FILE *data_file;

    data_file = fopen(input_filename,"r");

    if(data_file == NULL)
        return NULL;

    number_of_cords = compute_number_of_cord(data_file);
    number_of_lines = compute_number_of_x(data_file);
    X = read_data_from_file(data_file,number_of_cords,number_of_lines);
    if(X == NULL) return NULL;
    fclose(data_file);

    *row = number_of_lines;
    *col = number_of_lines;
    if(purpose != jacobi){
        W = weighted_adjacency_matrix(X,number_of_lines,number_of_cords);
        if(W == NULL) return NULL;
        if(purpose == wam){
            return W;
        }

        D = diagonal_degree_matrix(W,number_of_lines);
        if(D == NULL) return NULL;
        if(purpose== ddg){
            free_memory_of_matrix(W,number_of_lines);
            return D_asMatrix(D,number_of_lines);
        }

        Lnorm = calculate_Lnorm(D,W,number_of_lines);
        if(Lnorm == NULL) return NULL;
        if(purpose ==lnorm){
            free(D);
            free_memory_of_matrix(W, number_of_lines);
            return Lnorm;
        }
        
        /*goal = spk*/
        eigen_values_vectors = Jacobi_algorithm(&res_row,&res_col,number_of_lines,Lnorm);
        if(eigen_values_vectors == NULL) return NULL;

        *row = res_row;
        *col = res_col;

        qsort(eigen_values_vectors,res_row,sizeof(double*),comparator);
        
        if(k_from_py==0){
            k = find_k(eigen_values_vectors,res_row);
            if(k ==1 || k == 0){
                return NULL;
            }
        }
        else
            k = k_from_py;
        
        
        /*U- containing the eigenvectors u1, . . . , uk of Lnorm as columns */
        k_eigenvectors = create_eigenvectors_matrix(eigen_values_vectors,res_col,k);
        if(k_eigenvectors == NULL)return NULL;
                
        T = normalize_U(k_eigenvectors,res_col-1,k);
        if(T == NULL) return NULL;
               
        *row = res_col-1;
        *col = k;
        return T;
    }

    /*purpose == jacobi*/ 
    symmetric = check_input_symmetric_matrix(X,number_of_lines,number_of_cords);
    if(symmetric == -1) return NULL;
    
    eigen_values_vectors = Jacobi_algorithm(&res_row,&res_col,number_of_lines,X);
    if(eigen_values_vectors == NULL) return NULL;

    Jacobi= transpose(eigen_values_vectors,res_row,res_col); 
    if(Jacobi == NULL) return NULL;

    *row = res_col;
    *col = res_row;

    fix_zeros(Jacobi,*col);

    return Jacobi;
     
}


int get_enum_val(char* str){
    if(!strcmp( "wam", str)) return wam;
    if(!strcmp( "ddg", str)) return ddg;
    if(!strcmp( "lnorm", str)) return lnorm;
    if(!strcmp( "jacobi", str)) return jacobi;
    return 5;
}

/* submit args to vars, returns 0 if succeed else 1 */
int submit_args_spkmeans(int argc, char **argv, char **input_file, int *goal){

    if (argc != 3){
        printf("Invalid Input!");
        return 1;
    }

    *goal = get_enum_val(argv[1]);
    if(*goal > 4){
        printf("Invalid Input!");
        return 1; 
    }
    
    *input_file = argv[2];
    
    return 0;
}

int main(int argc, char** argv){
    int goal,row,col,s;
    char* filename;
    double** mat;
    s = submit_args_spkmeans(argc,argv,&filename,&goal);
    if(s==1)
        return 1;

    mat = Spkmeans(&row,&col,filename,goal,0);
    if(mat == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    printmatrix(mat,row,col);
    free_memory_of_matrix(mat,row);
    return 0 ;
}


