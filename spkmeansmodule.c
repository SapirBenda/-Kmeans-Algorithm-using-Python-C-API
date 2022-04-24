#include <Python.h>
#include "spkmeansmodule.h"
#include "spkmeans.h"

/* free memory of matrix */
void free_memory_of_matrix1(double** matrix,int number_of_rows){
    int i;
    for(i=0; i<number_of_rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}

PyObject *c_matrix_to_python_matrix(int rows, int cols, double **c_mat) {
    PyObject *py_mat;
    PyObject *py_list;
    int row, col;

    py_mat = PyList_New(rows);

    for (row = 0; row < rows; row++) {
        py_list = PyList_New(cols);

        for (col = 0; col < cols; col++) {
            PyList_SetItem(py_list, col, Py_BuildValue("f", c_mat[row][col]));
        }

        PyList_SetItem(py_mat, row, py_list);
    }
    free_memory_of_matrix1(c_mat,rows);

    return py_mat;
}

PyObject *c_square_matrix_to_python_square_matrix(int order, double **c_mat) {
    return c_matrix_to_python_matrix(order, order, c_mat);
}

PyObject* fit(PyObject *self, PyObject *args){
    int k, goal;
    char *input_filename;
    double** matrix;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "iis", &k,&goal, &input_filename)) {
        printf("An Error Has Occurred2\n");
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    int rows;
    int cols;
    matrix = Spkmeans(&rows, &cols, input_filename,goal,k);
    if(matrix == NULL) return NULL;
    return c_matrix_to_python_matrix(rows, cols, matrix); /*  Py_BuildValue(...) returns a PyObject*  */
}

PyObject* fit_kmeans(PyObject *self, PyObject *args){
    int k,max_iter,epsilon;
    char *data_filename, *mus_filename;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "iidss", &k,&max_iter,&epsilon,&data_filename, &mus_filename)) {
        printf("An Error Has Occurred\n");
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    return Py_BuildValue("i", K_mean(k,max_iter,epsilon,data_filename,mus_filename)); /*  Py_BuildValue(...) returns a PyObject*  */
}

static PyMethodDef capiMethods[] = {
        {"Spkmeans",                   /* the Python method name that will be used */
                (PyCFunction) fit, /* the C-function that implements the Python function and returns static PyObject*  */
                METH_VARARGS,           /* flags indicating parametersaccepted for this function */
                        PyDoc_STR("spkmeans")}, /*  The docstring for the function */
        {"k_means",                   /* the Python method name that will be used */
                (PyCFunction) fit_kmeans, /* the C-function that implements the Python function and returns static PyObject*  */
                METH_VARARGS,           /* flags indicating parametersaccepted for this function */
                        PyDoc_STR("kmeans ++")}, /*  The docstring for the function */
        {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

/* This initiates the module using the above definitions. */
static PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp", /* name of module */
        NULL, /* module documentation, may be NULL */
        -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};

/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
