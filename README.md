# 2-D Wave Equation Solver and Simulator
![image](https://user-images.githubusercontent.com/88569855/155702140-11ecd0dc-85b5-4365-a39e-243908158e42.png)

This piece of software is designed to implement a 2-D Wave Solver using parallelising concepts from MPI. The software aims to solve a 2-D domain through the discretisation of its PDE in both time and space. The main input parameters include the duration of the run, the grid size and the boundary conditions(among others). There are three boundary conditions implemented for this software - Dirichlet, Neumann and Periodic and the user may pick from amongst these. The initial disturbance is by default to a single point type and can be further expanded to include multiple points or a sinusoidal wave. Finally, we have the post-processing and visualisations of these results in the form of animated gifs. Read the report linked [here](report/report.pdf) for further details on the analysis.

## Repository Hierarchies
All the source file including the main are contained in the `src` folder and the associated headers are in the `include` folder. The `serial` sub-folder in output folder contains grid outputs from serial runs as well as images of timing tests for various grid inputs, subfolder `loadbalancing` conatins work load balancing outputs to ensure all processes were equally engaged and finally the `parallel` subfolder coantains an example each of the the domain outputs for each of the 3 boundary conditions and most importantly the `hpc` folder conatins all the .o files from the multiple HPC runs. \
The `postprocessing` folder contains the jupter notebook used for post processing. The `animations` folder in the postprocessing directory is the one in which the postprocessing outputs are stores


## Installation and Dependencies
The user is expected to have the MPI libraries and paths loaded- either locally or on the HPC systems

## Set Parameters
The user is given then ability to set the required input parameters using a Parameters.txt file. The user must strictly follow the instructions on the first line of the file to load the input conditions correctly. The parameteres to be input(in this order) are : \
 - `Maxinumn Domain sizes` : imax jmax \
 - `Domain Dimensions` : x_max y_max \
 - `Start Time and End Time for the Iterations` : t t_max \
 - `Intervals to obtain output` : dt_out \
 - `Wave Speed` : c \
 - `Boundary Conditions` : boundary(0 = Diriclet, 1 = Neumann, 2 = Periodic) \
 - `Disturbance parameters` : disturbance(x location, y location, radius)\
```
600 600 10 10 0.0 30.0 0.04 1.0 2 3.0 3.0 1.0
```
## Run the Code
This code is designed to run on the supercomputer systems. \
Once the user has logged on to the HPC system, they can first compile the code by first loading the intel-suite/2019.4 and the mpi module, from the list of available modules then compiling the files with the following command : 
```
mpicxx TimingTest.cpp Domain.cpp MPIDatatype.cpp WaveSolver.cpp Main.cpp -o Main
```
The user is expected to have access to such a system and be able to load a pbs script containing the job decription(number of nodes, number of cpus and number of processes, see example [here](HPC_loading_example.pbs)). \
To test the code on the local system you can compile in a similar manner locally and run it using the following command:
```
mpiexec -n (number of nodes) ./Main
```

An example of the expected output on the HPC will look like the following : \

<img width="485" alt="hpcoutput" src="https://user-images.githubusercontent.com/88569855/155711333-ab6c7457-b682-4a4a-8a70-7b893441e6a9.png">



## Post Processing and Visualizations

In the postprocessing folder, open the jupyter notebook on your local system. Within the code navigate to the conditions you wish to obtain outputs for and use out compilation outputs to add in total number of outputs in `num_outputs`, the number of rows in `sub_row`, the number of columns in `sub_cols` and the number of processors used in `processors` . Look for commenting on the files for guidance. \

An example of the output(for periodic boundary condition) generated is as follows :

<img src="postprocessing\animations\animate_parallel_periodic.gif">

## Version History
See the linked [Version Files](version.md)
## Documentation
See the linked [documentation files](/documentation/documentation.pdf) for the pdf version of the documentation. These documentation files were generated using Doxygen and you can also render then in an HTML format.
## License
This project is licensed under the MIT License - see the [LICENSE.md](License.md) file for details

