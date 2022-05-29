Compiling:
All of the task codes can be compile by using the Makefile.
For MPI programs, please load mpi version 4.0.1 before compiling
The command to load that version is 'module load openmpi/4.0.1'
Use the command 'Make' to compile everything and 'Make clean' to delete all the compiled programs

Naming:
Each task code is named as TaskB-Y and Y is the variant.
Below is the variant description:
TaskB-1 : Static scheduling
TaslB-2 : Dynamic scheduling
** Difference between variants is detailed in the report.

Running the program
For MPI program
mpirun -n <num_threads> TaskB-Y <input_testcase_filename> <output_filename>
Eg: mpirun -n 16 ./TaskB-2 ../testcases/input_testcase_8 ../testcases/output

Testcases:
The testcase are all store in the testcases directory
The provided testcase has the following naming convention:
input_testcase_D , where D is the dimension
Provided dimensions are : 8,16,24,32,64,96,384,800,1000,1536
Unfornately all programs except for task1 only accepts squared matrix input