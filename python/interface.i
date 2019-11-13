%module msr
 %{
 /* Includes the header in the wrapper code */
 #include "../src/python_api.h"
 %}
 
 /* Parse the header file to generate wrappers */
 %include "../src/python_api.h"
