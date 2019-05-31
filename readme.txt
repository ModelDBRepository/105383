This is the readme for the model associated with

Cortassa S, Aon MA, O'Rourke B, Jacques R, Tseng HJ, Marban E, Winslow
RL (2006) An integrated model of cardiac mitochondrial energy
metabolism and calcium dynamics. Biophys J 84:2734-55


Our ECME model in C++ contains the source files and the CVODE version
that we use as integrator to run it.

In order to run it, you need to establish the model parameters in file
"parameter", the running parameter in file "control" and it requires
also the file "initial_conditions".

It will save the output in the files "states", "currents",
"derivatives" and "ss_conditions". All those files are ASCII (text)
files.

The model code is actually in the file "model.cpp" there are some
files like "exp.cpp" that correspond to older versions of the code and
I believe we do not use them any longer.

There are three modes of running: the numeric mode, the average mode
and the MCA mode. They differ in the output. I guess what you need to
use is the numeric mode that gives you the states, currents, etc as a
function of time.

In order to compile it we used Microsoft visual studio 2003. And as
accessory library we need the "Boost library" and the version of it
that we are using is: 1.31.0.  After you install the boost library you
may need to include the paths in the project if it is installed in a
different location than the default.

I am not including these here for reasons of room.

Please, do not hesitate to contact me when you need further
assistance.

Sonia Cortassa
