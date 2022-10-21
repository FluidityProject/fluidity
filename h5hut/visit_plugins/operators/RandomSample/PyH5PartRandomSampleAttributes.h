#ifndef PY_H5PARTRANDOMSAMPLEATTRIBUTES_H
#define PY_H5PARTRANDOMSAMPLEATTRIBUTES_H
#include <Python.h>
#include <H5PartRandomSampleAttributes.h>

//
// Functions exposed to the VisIt module.
//
void            PyH5PartRandomSampleAttributes_StartUp(H5PartRandomSampleAttributes *subj, void *data);
void            PyH5PartRandomSampleAttributes_CloseDown();
PyMethodDef    *PyH5PartRandomSampleAttributes_GetMethodTable(int *nMethods);
bool            PyH5PartRandomSampleAttributes_Check(PyObject *obj);
H5PartRandomSampleAttributes *PyH5PartRandomSampleAttributes_FromPyObject(PyObject *obj);
PyObject       *PyH5PartRandomSampleAttributes_NewPyObject();
PyObject       *PyH5PartRandomSampleAttributes_WrapPyObject(const H5PartRandomSampleAttributes *attr);
void            PyH5PartRandomSampleAttributes_SetDefaults(const H5PartRandomSampleAttributes *atts);
std::string     PyH5PartRandomSampleAttributes_GetLogString();
std::string     PyH5PartRandomSampleAttributes_ToString(const H5PartRandomSampleAttributes *, const char *);

#endif

