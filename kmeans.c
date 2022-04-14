#include "common.h"

#include <math.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct node
{
    double* xi;
    struct node* next;
    struct node* prev;
} node;
typedef struct mu{
    double* mui;
    node* xi_list;
} mu;

typedef struct change{
    node* xi_node;
    int old_mu_index;
    int new_mu_index;
} change;

#define DEF_MAX_ITER 300
#define EPSILON 0.000

/* function decleration*/
/*double** read_all_lines_to_X(int number_of_lines, FILE* fp, int number_of_cord );*/
mu* initialze_mus_array(double** X, int K, int number_of_cord);
int argmin(double* xi, mu* mus, int K, int number_of_cords);
void add_x_to_mu(node* xi_node, mu* mui);
void delete_x_from_mu(node* xi_node, mu* mui);
void swap_nodes(node* xi_node, int old_mu, int new_mu, mu* mus);
double update_mus(mu* mus, int K, int number_of_cords);
int initial_xi_liked_list( int number_of_lines, int K, int number_of_cords, mu* mus, double** X);
change* update_changes_array(mu* mus, change* change_array, int number_of_cords, int K );
void implementing_changes(mu* mus, int number_of_lines, change* change_array);
int K_mean(int K, int max_iter, double epsilon, char* data_filename, char* mus_filename);
void free_memory(double** X,double** Y, mu* mus, int num_of_X, int k);

/* free memory at the end of k_means */
void free_memory(double** X,double** Y, mu* mus, int num_of_X, int k){
    node* curr;
    node* step;
    int i;
    for (i = 0; i < num_of_X; i++){
        free(X[i]);
    }
    free(X);

    for (i = 0; i < k; i++){
        free(Y[i]);
    }
    free(Y);

    for (i = 0; i < k; i++){
        free(mus[i].mui);
        curr = mus[i].xi_list;
        while(curr){
            step = curr->next;
            free(curr);
            curr = step;
        }
    }
    free((void*)mus);
}

/* create mu_array from first K xi from X
 output: return mu_array - the firt K xi from X */ 
mu* initialze_mus_array(double** X, int K, int number_of_cord){
    mu* mus = (mu*)calloc(K,sizeof(mu));
    int i, word;
    for(i =0; i< K; i++){
        mus[i].mui= (double*)calloc(number_of_cord,sizeof(double));
        if(mus[i].mui == NULL){
            printf("An Error Has Occurred");
            return NULL;
        }
        mus[i].xi_list = NULL;
        for(word =0; word < number_of_cord; word++){
            mus[i].mui[word] = X[i][word];
        }
    }
    return mus;
}


/* outpot: the index of its clostest mu of xi*/
int argmin(double* xi, mu* mus, int K, int number_of_cords){
   double min_sum, sum;
   int min_index=0, mu_index, cord;

    for(mu_index =0; mu_index < K; mu_index++ ){
       sum =0;
       for(cord = 0; cord< number_of_cords; cord++){
           sum += pow((xi[cord] - mus[mu_index].mui[cord]),2);
       }
       if (mu_index ==0) min_sum = sum;
        else if(sum<min_sum) {
            min_sum = sum;
            min_index = mu_index;
        } 
    }
    return min_index;
}

/* output: add node of xi to xi_list of mui*/
void add_x_to_mu(node* xi_node, mu* mui){
    xi_node->prev = NULL;
    if(mui->xi_list == NULL){ /* xi_list empty*/
        mui->xi_list = xi_node; 
        xi_node->next = NULL;
    }else{
    mui->xi_list->prev = xi_node;
    xi_node->next = mui->xi_list;
    mui->xi_list = xi_node;
    }
}

/* output: remove node ox xi frox xi_list of mui*/
void delete_x_from_mu(node* xi_node, mu* mui){
    if(xi_node->prev == NULL){ /* first node in list*/
        if(xi_node->next == NULL){ /* xi_node is the only node in list*/
            mui->xi_list = NULL;
        }else{
            xi_node->next->prev = NULL;
            mui->xi_list= xi_node->next;
        }
    } 
    else if (xi_node->next ==NULL){ /* last node in list*/
        xi_node->prev->next=  NULL;
    }
    else{ 
        xi_node->prev->next = xi_node->next;
        xi_node->next->prev = xi_node->prev;
    }
}

/* output: remote xi_node from old_mu list and add it to the start of new_mu list*/
void swap_nodes(node* xi_node, int old_mu, int new_mu, mu* mus){
    delete_x_from_mu(xi_node, &mus[old_mu]);
    add_x_to_mu(xi_node,&mus[new_mu]);
}


/* output: compute and update new_mus_array and returns delta max */
double update_mus(mu* mus, int K, int number_of_cords){
    double deltamax =0, delta;
    int mu_index, xi_list_len, cord;
    double* new_mu;
    double* old_mu;
    node* curr_xi;

    for(mu_index =0; mu_index< K; mu_index++){
        xi_list_len =0;
        new_mu = (double*)calloc(number_of_cords, sizeof(double));
        if(new_mu == NULL){
            printf("An Error Has Occurred");
            return -1;
        }
        curr_xi = mus[mu_index].xi_list;
        while (curr_xi != NULL)
        {
            for(cord =0; cord<number_of_cords; cord++){
                new_mu[cord] += curr_xi->xi[cord];
            }
            xi_list_len++;
            curr_xi = curr_xi->next;
        }
        /* we assume xi_list_len != 0 becuse mui dosent have an empty list as noted in the forum */
        for(cord =0; cord< number_of_cords; cord++){
            new_mu[cord] = new_mu[cord]/xi_list_len;
        }
        old_mu =  mus[mu_index].mui;
        delta = euqlide_norm(old_mu, new_mu, number_of_cords);
        if (delta > deltamax) deltamax= delta;
        free(old_mu);
        mus[mu_index].mui = new_mu;
    }
    return deltamax;
}

/* output: create output.txt and write mus_array in it*/
int write_to_outputfile(mu* mus, int K, FILE* fp_out, int number_of_cords ){
    int mu, cord;
    for(mu =0; mu <K; mu++){
        for(cord =0; cord< number_of_cords; cord++){
             fprintf(fp_out,"%.4f",mus[mu].mui[cord]);
            if (cord == number_of_cords -1 && mu != K-1 ){
                fprintf(fp_out,"%s","\n");
            }else if (cord != number_of_cords-1){
                fprintf(fp_out,"%s",",");
            }     
        }
    }
    fprintf(fp_out,"%s","\n");
    return 0; 
}


/*create node for each xi and insert it to it's closest mu linked list */
int initial_xi_liked_list( int number_of_lines, int K, int number_of_cords, mu* mus, double** X){
    int xi, new_mu;
    node* new_xi_node;
    for(xi = 0; xi< number_of_lines; xi++){
        new_mu = argmin(X[xi],mus,K,number_of_cords);
        new_xi_node = (node*)malloc(sizeof(node));
        if(new_xi_node == NULL ) {
            printf("An Error Has Occurred");
            return 0;
        }
        new_xi_node->xi = X[xi];
        add_x_to_mu(new_xi_node, &mus[new_mu]);
    }
    return 1;
}

/* for evry xi save its old mu & new mu in changes array */ 
change* update_changes_array(mu* mus, change* change_array, int number_of_cords, int K ){
    int new_mu, mu_index, xi_counter =0;
    node* curr_xi;
    for(mu_index =0; mu_index< K; mu_index++){
        curr_xi = mus[mu_index].xi_list;
        while (curr_xi != NULL){
            new_mu = argmin(curr_xi->xi,mus,K,number_of_cords);
            change_array[xi_counter].new_mu_index = new_mu;
            change_array[xi_counter].old_mu_index=mu_index;
            change_array[xi_counter].xi_node = curr_xi;
            xi_counter++;
            curr_xi = curr_xi->next;
        }
    }
    return change_array;
}

/* update mu's linked list according changes array */ 
void implementing_changes(mu* mus, int number_of_lines, change* change_array){
    int change_num,new_mu,old_mu;
    for(change_num = 0; change_num< number_of_lines; change_num++){
        new_mu = change_array[change_num].new_mu_index;
        old_mu = change_array[change_num].old_mu_index;
        if(new_mu != old_mu ){
            swap_nodes(change_array[change_num].xi_node, old_mu,new_mu, mus);
        }
    }
}


int K_mean(int K, int max_iter, double epsilon, char* data_filename, char* mus_filename){

    double maxdelta = epsilon;
    int iter =1, number_of_cords, number_of_lines, memory_alocate;
    change* change_array;
    double** X, **Y; 
    mu* mus;
    FILE *data_file, *mus_file, *final_mus;

    data_file = fopen(data_filename,"r");
    mus_file = fopen(mus_filename,"r");

    number_of_cords = compute_number_of_cord(data_file);
    number_of_lines = compute_number_of_x(data_file);

    if(K > number_of_lines){
        printf("Invalid Input!");
        return 1;
    }
    X = read_data_from_file(data_file,number_of_cords,number_of_lines);
    if(X == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    Y = read_data_from_file(mus_file,number_of_cords,K);
    mus = initialze_mus_array(Y,K,number_of_cords);
    if(mus == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    memory_alocate= initial_xi_liked_list(number_of_lines, K, number_of_cords, mus,X);
    if(memory_alocate == 0) {
        printf("An Error Has Occurred");
        return 1;
    }
    fclose(mus_file);
    fclose(data_file);
    
    while (iter <= max_iter && maxdelta >= epsilon){
        change_array =(change*) calloc(number_of_lines,sizeof(change));
        if(change_array == NULL ) {
            printf("An Error Has Occurred");
            return 1;
        }
        change_array = update_changes_array(mus,change_array,number_of_cords,K); /* update changes_array according new mus */
        implementing_changes(mus, number_of_lines,change_array); /* link evrey xi to its new mu according to changes_array */
        maxdelta = update_mus(mus,K,number_of_cords); /* compute new max delts */
        if(maxdelta == -1) {
            printf("An Error Has Occurred");
            return 1;
        }
        iter++;
        /*free change array*/
        free(change_array);   
    }

    final_mus = fopen("data_from_kmeans_c.txt","w");
    write_to_outputfile(mus,K,final_mus,number_of_cords);
    fclose(final_mus);
    free_memory(X,Y,mus,number_of_lines,K);
    return 0;
}