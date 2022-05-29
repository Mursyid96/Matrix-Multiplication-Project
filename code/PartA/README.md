Compiling:
All of the task codes can be compile by using the Makefile.
For MPI programs, please load mpi version 4.0.1 before compiling
The command to load that version is 'module load openmpi/4.0.1'
Use the command 'Make' to compile everything and 'Make clean' to delete all the compiled programs

Naming:
Each task code is named as TaskX-Y where X is the task number and Y is the variant.
Below is the variant description:
Task1
1: Static scheduling
2: Dynamic scheduling
Task2
1: Blocking send and recv
2: Non-blocking Isend and Irecv
Task3
1: Blocking send and recv
2: Non-blocking Isend and Irecv
**Since task 4 has two sub-tasks, the naming convention is a bit different.
**Task4A is SUMMA and Task4B is cannon's algorithm
Task4A
1: Buffered
2: Direct
Task4B
1: Buffered
2: Direct
** Difference between variants is detailed in the report.

Running the program
For normal C program:
TaskX-Y <num_threads> <input_testcase_filename> <output_filename>
Eg: ./Task1-2 16 ../testcases/input_testcase_8 ../testcases/output

For MPI program (Task2 & Task3)
mpirun -n <num_threads> TaskX-Y <input_testcase_filename> <output_filename>
Eg: mpirun -n 16 ./Task3-2 ../testcases/input_testcase_8 ../testcases/output

Testcases:
The testcase are all store in the testcases directory
The provided testcase has the following naming convention:
input_testcase_D , where D is the dimension
Provided dimensions are : 8,16,24,32,64,96,384,800,1000,1536
Unfornately all programs except for task1 only accepts squared matrix input
