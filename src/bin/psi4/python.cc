// Include Python's header file.
// Use "python-config --includes" to determine this location.
#include <cstdio>
#include <boost/python.hpp>
//#include <Python.h>
#include <psiconfig.h>
#include <sstream>
#include "script.h"

using namespace psi;
using namespace boost;

namespace psi {
    extern void psiclean(void);
    extern FILE *outfile;
}

bool py_psi_input()
{
    // Need to modify input to not take argc and argv
    // And Options object to be global.
    // input::input3();
    return true;
}

char const* py_psi_version()
{
    return PSI_VERSION;
}

void py_psi_clean()
{
    psiclean();
}

BOOST_PYTHON_MODULE(psi)
{
    using namespace boost::python;
    def("version", py_psi_version);
    def("clean", py_psi_clean);

    // modules
    def("input",py_psi_input);
}

//static PyObject *
//py_psi_version(PyObject *self, PyObject *args)
//{
//    return Py_BuildValue("s", PSI_VERSION);
//}
//
//static PyObject *
//py_psi_input(PyObject *self, PyObject *args)
//{
////    int result = input::input(options, 0, NULL);
//    return Py_BuildValue("i", 0);
//}
//
//static PyObject *
//py_psi_clean(PyObject *self, PyObject *args)
//{
//    psiclean();
//    return Py_BuildValue("i", 0);
//}
//
//static PyMethodDef PsiMethods[] = {
//    { "version", py_psi_version, METH_NOARGS, "Obtain version information from PSI."},
//    { "input",   py_psi_input,   METH_NOARGS, "Run the input module."},
//    { "clean",   py_psi_clean,   METH_NOARGS, "Run the psiclean module."},
//    { NULL, NULL, 0, NULL } /* Sentinel */
//};
//
//static struct PyModuleDef PsiModule = {
//    PyModuleDef_HEAD_INIT,
//    "psi",   /* name of the module */
//    NULL,    /* something with documentation */
//    -1,      /* size of the per-interpreter state of the module,
//                or -1 if the module keeps state in global variables. */
//    PsiMethods
//};
//
//PyMODINIT_FUNC
//PyInit_psi()
//{
//    return PyModule_Create(&PsiModule);
//}

Python::Python() : Script()
{
    
}

Python::~Python()
{
    
}

void Python::initialize()
{
}

void Python::finalize()
{
}

void Python::run(FILE *input)
{
    using namespace boost::python;
    char *s;
    if (input == NULL)
        return;
    if (!Py_IsInitialized()) {
        s = strdup("psi");
//        PyImport_AppendInittab(s, initpsi);
        Py_Initialize();
        #if PY_VERSION_HEX >= 0x03000000
        Py_SetProgramName(L"psi");
        #else
        Py_SetProgramName(s);
        #endif
    }
    if (Py_IsInitialized()) {
        char line[256];
        std::stringstream file;
        while(fgets(line, sizeof(line), input)) {
            file << line;
        }
        printf("Input file:\n%s", file.str().c_str());
        str strStartScript(file.str().c_str());

        try {
            PyImport_AppendInittab(s, initpsi);
            object objectMain(handle<>(borrowed(PyImport_AddModule("__main__"))));
            object objectDict = objectMain.attr("__dict__");
            s = strdup("import psi; from psi import *;");
            PyRun_SimpleString(s);
            
            object objectScriptInit = exec( strStartScript, objectDict, objectDict );
        }
        catch (error_already_set const& e)
        {
            PyErr_Print();
        }
    } else {
        fprintf(stderr, "Unable to run Python input file.\n");
        return;
    }
}