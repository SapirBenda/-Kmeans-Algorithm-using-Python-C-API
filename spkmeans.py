import sys
import numpy as np
import pandas as pd
import mykmeanssp
import enum

MAX_ITER = 300
EPSILON =0

class goals(enum.Enum):
    spk =0
    wam=1
    ddg =2
    lnorm=3
    jacobi =4


def print_file(filename):
    f = open(filename)
    for line in f:
        print(line[:-1])
    f.close()
    
def print_jacobi_matrix(matrix):
    for row in range(len(matrix)):
        if(row ==0):
            line = [(str("{:.4f}".format(cord))).replace("-0.0000","0.0000") for cord in matrix[row]]
        else: 
            line = [(str("{:.4f}".format(cord))) for cord in matrix[row]]
        print(",".join(line))
        
def print_matrix(matrix):
    for row in matrix:
        row = [(str("{:.4f}".format(cord))) for cord in row]
        print(",".join(row))



def calc_min_dist(rows, mus, curr_num_of_mus, row_by_index_d):
    d_arr = np.zeros(rows)
    for l in range(rows):
        dl = min([np.sum((row_by_index_d[l] - mus[i]) ** 2) for i in range(curr_num_of_mus)])
        d_arr[l] = dl
    return d_arr


def calc_probs(d_arr):
    d_sum = np.sum(d_arr)
    return d_arr / d_sum


def create_row_by_index_d(np_data):
    d = {}
    for row in np_data:
        d[int(row[0])] = row[1:]
    return d


def create_data_file(filename, data):
    file = open(filename, "w")
    for row in data:
        mu = [(str(cord)) for cord in row]
        file.write(','.join(mu) + "\n")
    file.close()
    
def create_mus_file(filename, data):
    file = open(filename, "w")
    for row in data:
        lst_row = row.tolist()
        mu = [(str(cord)) for cord in lst_row]
        file.write(','.join(mu) + "\n")
    file.close()


def find_mus(indexes, cols, np_data,k, rows):
    mus_indexes = []
    mus = np.zeros([k, cols - 1])
    chosen_xi = np.random.choice(indexes)
    mus_indexes.append(chosen_xi)
    row_by_index_d = create_row_by_index_d(np_data)

    mus[0] = row_by_index_d[chosen_xi]
    for i in range(1, k):
        dist_arr = calc_min_dist(rows, mus, i, row_by_index_d)
        probs_arr = calc_probs(dist_arr)
        chosen_xi = np.random.choice(indexes, p=probs_arr)
        mus[i] = row_by_index_d[chosen_xi]
        mus_indexes.append(chosen_xi)

    mus_indexes_str = [(str(cord)) for cord in mus_indexes]
    return mus_indexes_str, mus

def create_indexes(np_data):
    indexes = np_data[:, 0].astype(int)
    indexes.sort()
    return indexes


def k_means_pp(matrix):
    np.random.seed(0)
    matrix = np.array(matrix)
    rows_number = len(matrix)
    k = len(matrix[0])
    indexes = np.arange(rows_number)
    matrix = np.insert(matrix,0,indexes,axis =1)
    cols_number = len(matrix[0])
    mus_indexes_str, mus = find_mus(indexes, cols_number, matrix,k, rows_number)
    create_mus_file("mus_file.txt", mus)
    print(','.join(mus_indexes_str))
    return mus


def submit_args():
    if len(sys.argv) != 4:
        print("Invalid Input!",end="")
        return 1
    try:
        k = int(sys.argv[1])
        goal = sys.argv[2]
        intgoal = int(goals[goal].value)
        input_filename = sys.argv[3]
        f_input_1 = open(input_filename,"r")
        length_file = len(f_input_1.readlines())
        if k == 1 or k < 0 or k>=length_file:
            print("Invalid Input!",end="")
            return 1
        f_input_1.close()

    except (ValueError,KeyError ,OSError):
        print("Invalid Input!",end="")
        return 1
    return k, intgoal, input_filename


def main():
    args = submit_args()
    mus_filename = 'mus_file.txt'
    data_filename = "matrix_from_C.txt"
    if args == 1:
        return 1
    k, goal, input_filename = args
    try:
        matrix = mykmeanssp.Spkmeans(k,goal,input_filename)
        if matrix!= None:
            if goal == 0: # goal = spk
                create_data_file(data_filename,matrix)
                matrix = k_means_pp(matrix)
                k = len(matrix[0])
                success = mykmeanssp.k_means(k,MAX_ITER,EPSILON,data_filename,mus_filename)
                if(success != 0): return 1
                else: print_file("data_from_kmeans_c.txt")
            else:
                matrix = np.array(matrix)
                if goal == 4 : # goal = jacobi
                    print_jacobi_matrix(matrix)
                else: # gaol = wam / ddg / lnorm
                    print_matrix(matrix)    
            return 0   
        else: # matrix from Spkmeans = NULL
            print('An Error Has Occurred',end="")
            return 1
    except Exception as e:
        print('An Error Has Occurred',end="")
        return 1
    

if __name__ == '__main__':
    main()
