#ifndef FUNCTIONMODULE_H
#define FUNCTIONMODULE_H

static PyObject* py_f_roots(PyObject* self, PyObject* args);
static PyObject* py_f_plot(PyObject* self, PyObject* args);
void fill_do_py(PyObject *d_o_py,const struct DataOutput_raw d_o);
void parse_di_py(PyObject *di_py,struct DataInput *di);
static PyObject* py_f_plot(PyObject* self, PyObject* args);
void set_int_from_py(PyObject *obj, char *name, int *to_set);
void set_double_from_py(PyObject *obj, char *name, double *to_set);
void set_double_array_from_np1D(PyObject *obj, char *name, double *to_set);
void set_double_array_from_list(PyObject *obj, char *name, double **to_set, int len);
void set_double_np1D_from_array(PyObject *obj,char *name, const double *data);
void set_int_array_from_np1D(PyObject *obj, char *name, int *to_set);
void set_int_np1D_from_array(PyObject *obj,char *name, const int *data);
void set_list_double_from_array(PyObject *list, double *array, int len);
void set_double_np2d_from_matrix(PyObject *obj, char *name, double **to_set);
void fill_dp_py(PyObject *dp_py, struct Data_plot dp);
#endif
