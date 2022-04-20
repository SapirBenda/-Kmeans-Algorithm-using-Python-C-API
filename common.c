#include "common.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> 


/* create an array of pointers to xi calld X
 output: read file and update  all_x_array -  an array of vertex , initial all mus to null*/
double** read_data_from_file(FILE* fp, int number_of_cord, int number_of_lines){
    int line_number;
    char *curr_number_start, *line;
    double* xi;
    double** X = (double**)calloc(number_of_lines , sizeof(double*));
    rewind(fp);
    line = (char*)calloc(2000000,sizeof(char));
    if (line == NULL || X == NULL){
        return NULL;
    }
    for(line_number =0 ; line_number < number_of_lines; line_number++){
        curr_number_start = line;
        if(fscanf(fp, "%s",line)== EOF){
            return NULL;
        }
        xi = (double*) calloc(number_of_cord,sizeof(double));
        if (xi == NULL){
            return NULL;
        }
        X[line_number]= read_line_into_xi(number_of_cord,curr_number_start,xi);
    }
    free(line);
    return X;

}

/* read line from file and save it to xi */
double* read_line_into_xi(int number_of_cord, char* curr_number_start, double* xi ){
    char* curr_comma;
    int word_num;
    for ( word_num =0; word_num<number_of_cord -1; word_num++){
        curr_comma = strchr((char*)curr_number_start,','); /* point to the first , in line */
        *curr_comma = '\0'; 
        xi[word_num] = atof(curr_number_start); /* insert word to xi */
        curr_number_start = ++curr_comma; /* update line to tne next char after , */
    }
    xi[number_of_cord-1] = atof(curr_number_start);
    return xi;
}

/*euqlide_norm
output: for every cord (sum of (xi[cord]**2 -xj[cord])**2 )**0.5*/
double euqlide_norm(double* xi, double* xj, int number_of_cords){
    double sum =0, power;
    int i;
    for(i =0; i< number_of_cords ; i++){
        power = pow(fabs(xi[i] - xj[i]),2);
        sum += power;
    }
    /*sum = pow(sum,0.5);*/
    sum = sqrt(sum);
    return sum;
}


/*output: returns number of lines in file fp*/
int compute_number_of_x(FILE *fp ){
    char ch;
    int xi_counter =0;
    rewind(fp); /* goes to the beginning of th file*/
     do{
        ch = fgetc(fp);
        if(ch =='\n') xi_counter ++;
    } while (ch != EOF);
    return xi_counter ;
}

/* output: return the number of the coordinates in every line in fp*/
int compute_number_of_cord(FILE *fp){
    char ch = '0';
    int cords_counter = 1;
    rewind(fp); /* goes to the beginning of th file*/
    while (ch!='\n'){
        ch = fgetc(fp);
        if(ch == ','){
           cords_counter++;
        } 

    }
  return cords_counter;
}
