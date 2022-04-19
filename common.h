#ifndef E2C1FEE5_A883_4A38_839E_445DBDFBB306
#define E2C1FEE5_A883_4A38_839E_445DBDFBB306


#include <stdio.h>

double** read_data_from_file(FILE* fp, int number_of_cord, int number_of_lines);
double* read_line_into_xi(int number_of_cord, char* curr_number_start, double* xi );
double euqlide_norm(double* old_mu, double* new_mu, int number_of_cords);
int compute_number_of_x(FILE *fp );
int compute_number_of_cord(FILE *fp);

#endif /* E2C1FEE5_A883_4A38_839E_445DBDFBB306 */
