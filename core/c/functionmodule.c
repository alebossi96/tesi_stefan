//problema i_bound
#define PY_SSIZE_T_CLEAN
#include "python3.8/Python.h"
#include<stdio.h>
#include "function.h"
#include<stdlib.h>
#include "functionmodule.h"
#include <numpy/arrayobject.h>
#define PRINT 1
static PyObject* py_f_roots(PyObject* self, PyObject* args) {
  PyObject *di_py;
  PyObject *d_o_py;
  struct DataInput di;
  struct DataOutput_raw d_o;
  if (!PyArg_ParseTuple(args, "OO", &di_py, &d_o_py))
    return NULL;
  parse_di_py(di_py, &di);
  //inizializzare d_o arrays
  d_o.kappa_j= (double*)calloc(di.n_kj,sizeof(double));
  d_o.i_bound = (int*)calloc(di.n_kj*3,sizeof(int));
  d_o.kappa_z0 = (double*)calloc(di.n_kj*3*di.ni,sizeof(double));
  set_double_array_from_np1D(d_o_py, "kappa_j",d_o.kappa_j);
  set_double_array_from_np1D(d_o_py, "kappa_z0",d_o.kappa_z0);
  set_int_array_from_np1D(d_o_py, "i_bound",d_o.i_bound);
  f_roots(di, &d_o);
  //TODO return & fill do_py
  fill_do_py(d_o_py, d_o);
  //TODO fare i free dei malloc
  return Py_None;
}


static PyObject* py_f_plot(PyObject* self, PyObject* args) {
  PyObject *di_py;
  PyObject *d_o_py;
  PyObject *dp_py;
  int i_r;
  struct DataInput di;
  struct DataOutput_raw d_o;
  struct Data_plot dp_inner;
  if (!PyArg_ParseTuple(args, "OOOi", &di_py, &d_o_py, &dp_py, &i_r))
    return NULL;
  parse_di_py(di_py, &di);
  //inizializzare d_o arrays con 
  d_o.kappa_j= (double*)calloc(di.n_kj,sizeof(double));
  d_o.i_bound = (int*)calloc(di.n_kj*3,sizeof(int));
  d_o.kappa_z0 = (double*)calloc(di.n_kj*3*di.ni,sizeof(double));
  set_double_array_from_np1D(d_o_py, "kappa_z0",d_o.kappa_z0);
  set_int_array_from_np1D(d_o_py, "i_bound",d_o.i_bound);
  set_double_array_from_np1D(d_o_py, "kappa_j",d_o.kappa_j);
  for (int j=0; j<N_REC_MAX; j++)
	{
	for (int i=0; i<N_PUNTI_TPSF_MAX; i++)
		{
		dp_inner.r_tpsf[i][j]=0.;
		dp_inner.t_tpsf[i][j]=0.;
		}	
	}

  f_plot(di,&d_o,&dp_inner,i_r);
  //TODO return & fill do_py
  fill_do_py(d_o_py, d_o);
  fill_dp_py(dp_py,dp_inner);
  //TODO fare i free dei malloc
  return Py_None;
}

void fill_do_py(PyObject *d_o_py,const struct DataOutput_raw d_o){

  set_int_np1D_from_array(d_o_py,"i_bound", d_o.i_bound);
  set_double_np1D_from_array(d_o_py,"kappa_z0", d_o.kappa_z0);
  set_double_np1D_from_array(d_o_py,"kappa_j", d_o.kappa_j);
  free(d_o.kappa_j);
  free(d_o.i_bound);
  free(d_o.kappa_z0);

}

void fill_dp_py(PyObject *dp_py, struct Data_plot dp){
    PyArrayObject *data_py_r_tpsf = (PyArrayObject *) PyObject_GetAttrString(dp_py, "r_tpsf");
    double *data_r_tpsf = (double *) PyArray_DATA(data_py_r_tpsf); //sono tutti in filap i
    int cont = 0;
    for(int i = 0; i<PyArray_DIMS(data_py_r_tpsf)[0]; i++){
        for(int j = 0; j<PyArray_DIMS(data_py_r_tpsf)[1]; j++){
            data_r_tpsf[cont] = dp.r_tpsf[i][j];
            cont++;
        }
    }
    PyArrayObject *data_py_t_tpsf = (PyArrayObject *) PyObject_GetAttrString(dp_py, "t_tpsf");
    double *data_t_tpsf = (double *) PyArray_DATA(data_py_t_tpsf); //sono tutti in fila
    cont = 0;
    for(int i = 0; i<PyArray_DIMS(data_py_t_tpsf)[0]; i++){
        for(int j = 0; j<PyArray_DIMS(data_py_t_tpsf)[1]; j++){
            data_t_tpsf[cont] = dp.t_tpsf[i][j];
            cont++;
        }
    }
    //set_double_np2d_from_matrix(dp_py, "r_tpsf", dp.r_tpsf);
    //set_double_np2d_from_matrix(dp_py, "t_tpsf", dp.t_tpsf);
}


void parse_di_py(PyObject *di_py,struct DataInput *di){
    set_int_from_py(di_py, "i_type", &di->i_type);
    set_int_from_py(di_py, "n_rec", &di->n_rec);
    set_double_array_from_np1D(di_py, "rec",di->rec);//TODO allocare spazio!
    set_double_from_py(di_py, "R", &di->R);
    set_double_from_py(di_py, "n0", &di->n0);
    set_double_from_py(di_py, "n1", &di->n1);
    set_double_from_py(di_py, "ne", &di->ne); 
    set_double_from_py(di_py, "c", &di->c);
    set_double_from_py(di_py, "tmin", &di->tmin); 
    set_double_from_py(di_py, "dr", &di->dr);
    set_double_from_py(di_py, "rmin", &di->rmin); 
    set_double_from_py(di_py, "xacc", &di->xacc);
    set_double_from_py(di_py, "precisione", &di->precisione); 
    set_int_from_py(di_py, "ni", &di->ni);
    set_int_from_py(di_py, "n_kj", &di->n_kj); 
    set_int_from_py(di_py, "nh", &di->nh);
    set_int_from_py(di_py, "nhi", &di->nhi); 
    set_double_from_py(di_py, "n0e", &di->n0e);
    set_double_from_py(di_py, "n1e", &di->n1e); 
    set_double_from_py(di_py, "A0e", &di->A0e);
    set_double_from_py(di_py, "A1e", &di->A1e); 
    set_double_from_py(di_py, "pi", &di->pi);
    set_double_from_py(di_py, "v0", &di->v0);     
    set_double_from_py(di_py, "v1", &di->v1);
    set_double_from_py(di_py, "dt", &di->dt); 
    set_int_from_py(di_py, "n_tpsf", &di->n_tpsf);
    set_double_from_py(di_py, "z0", &di->z0); 
    set_double_from_py(di_py, "z1", &di->z1);
    set_double_from_py(di_py, "ua0", &di->ua0);
    set_double_from_py(di_py, "ua1", &di->ua1); 
    set_double_from_py(di_py, "ud0", &di->ud0);
    set_double_from_py(di_py, "ud1", &di->ud1);
    set_double_from_py(di_py, "usR0", &di->usR0); 
    set_double_from_py(di_py, "usR1", &di->usR1);
}


void set_int_from_py(PyObject *obj, char *name, int *to_set){
    PyObject *attr = PyObject_GetAttrString(obj, name);
    PyArg_Parse(attr, "i", to_set);
    
}
void set_double_from_py(PyObject *obj, char *name, double *to_set){
    PyObject *attr = PyObject_GetAttrString(obj, name);
    PyArg_Parse(attr, "d", to_set);    
}
void set_double_np2d_from_matrix(PyObject *obj, char *name, double **to_set){
    PyArrayObject *data_py = (PyArrayObject *) PyObject_GetAttrString(obj, name);
    double *data = (double *) PyArray_DATA(data_py); //sono tutti in fila
    int cont = 0;
    for(int i = 0; i<PyArray_DIMS(data_py)[0]; i++){
        for(int j = 0; j<PyArray_DIMS(data_py)[1]; j++){
            data[cont] = to_set[i][j];
            cont++;
        }
    }
}
void set_double_array_from_np1D(PyObject *obj, char *name, double *to_set){
    PyArrayObject *data_py = (PyArrayObject *) PyObject_GetAttrString(obj, name);
    //TODO vediamo se va bene!
    double *data = (double *) PyArray_DATA(data_py);
    for (size_t i = 0; i<PyArray_DIMS(data_py)[0]; i++){ 
        to_set[i] = data[i];
    }
}
void set_int_array_from_np1D(PyObject *obj, char *name, int *to_set){
    PyArrayObject *data_py = (PyArrayObject *) PyObject_GetAttrString(obj, name);
    //TODO vediamo se va bene!
    long *data = (long *) PyArray_DATA(data_py);
    for (size_t i = 0; i<PyArray_DIMS(data_py)[0]; i++){ 
        to_set[i] = (int) data[i];
    }
}
void set_double_array_from_list(PyObject *obj, char *name, double **to_set, int len){
    PyObject *attr = PyObject_GetAttrString(obj, name);
    for(size_t i = 0; i<len; i++){
        PyObject * item = PyList_GET_ITEM(obj,i);//TODO non devo incrementare?
        PyArg_Parse(item, "d", *to_set[i]);
    } 
}
void set_double_np1D_from_array(PyObject *obj,char *name, const double *data){
    PyArrayObject *to_set = (PyArrayObject *) PyObject_GetAttrString(obj, name);
    double *data_out = PyArray_DATA(to_set);
    for (size_t i = 0; i<PyArray_DIMS(to_set)[0]; i++)
        data_out[i] = data[i];
}

void set_int_np1D_from_array(PyObject *obj,char *name, const int *data){
    PyArrayObject *to_set = (PyArrayObject *) PyObject_GetAttrString(obj, name);
    long *data_out = PyArray_DATA(to_set);
    for (size_t i = 0; i<PyArray_DIMS(to_set)[0]; i++)
        data_out[i] = (long) data[i];
}

void set_list_double_from_array(PyObject *list, double *array, int len){
    for (size_t i = 0; i<len; i++)
        PyList_Append(list, Py_BuildValue("d",  array[i]));
}
// Define a collection of methods callable from our module
static PyMethodDef PyFunctionMethods[] = {
  {"f_roots", py_f_roots, METH_VARARGS, "Funzione che calcola le radici dell'equazione trascendente"},
  {"f_plot", py_f_plot, METH_VARARGS, "questa funzione costruisce il vettore del profilo temporale uttilizzando le radici calcolate dalla funzione f_roots"},
  {NULL,NULL}
};

// Module definition
static struct PyModuleDef functionmodule = {
  PyModuleDef_HEAD_INIT,
  "functionmodule",
  "This module calculates the roots of transcendent function",
  -1,
  PyFunctionMethods
};


PyMODINIT_FUNC PyInit_functionmodule(void) 
{ 
  return PyModule_Create(&functionmodule); 
} 
