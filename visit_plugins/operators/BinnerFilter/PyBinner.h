#ifndef PY_BINNER_H
#define PY_BINNER_H
#include <Python.h>
#include <BinnerFilter.h>

//
// Functions exposed to the VisIt module.
//
void            PyBinner_StartUp(Binner *subj, void *data);
void            PyBinner_CloseDown();
PyMethodDef    *PyBinner_GetMethodTable(int *nMethods);
bool            PyBinner_Check(PyObject *obj);
Binner *PyBinner_FromPyObject(PyObject *obj);
PyObject       *PyBinner_NewPyObject();
PyObject       *PyBinner_WrapPyObject(const Binner *attr);
void            PyBinner_SetDefaults(const Binner *atts);
std::string     PyBinner_GetLogString();
std::string     PyBinner_ToString(const Binner *, const char *);

#endif

