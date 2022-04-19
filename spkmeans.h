#ifndef F9F1FA83_BEE0_41EF_928F_7D521D8EC307
#define F9F1FA83_BEE0_41EF_928F_7D521D8EC307

double** Spkmeans (int* rows, int* cols, char* input_filename, int purpose, int k_from_py);
int K_mean(int K, int max_iter, double epsilon, char* data_filename, char* mus_filename);
int initial_rows(double** W, int number_of_lines,int number_of_columns);


#endif /* F9F1FA83_BEE0_41EF_928F_7D521D8EC307 */
