
#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <chrono> 
#include <sstream>
#include <fstream>

#include <string>
#include "Parameters.h"
#include "WaveSolver.h"

#include "MPIDatatype.h"
#include "TimingTest.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------|
//  DECLARING SOFTWARE PARAMETERS FOR USE IN MAIN
//----------------------------------------------------------------------------------------------------------------|
// Declaring processor parameters 
int id, p;
int tag_num = 0;

// Declaring spatial grid parameters
int imax, jmax ;
double y_max, x_max, dx, dy;

// Declaring time stepping parameters
double t_max;
double t , dt;

// Declaring parameters to print output to file
double t_out, dt_out;

// Declaring wave equation parameters
double c;

// Declaring boundary condition parameters
int boundary;					
//----------------------------------------------------------------------------------------------------------------|

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 10);

	//------------------------------------------------------------------------------------------------------------|
	//  READ FROM PARAMETERS TEXT FILE
	// -----------------------------------------------------------------------------------------------------------|

	fstream file;
	file.open("Parameters.txt", std::ios::in);
	// If file is open, read the text in the file and write it out
	if (file.is_open())
	{
		std::string line;
		// Skip the first line with instructions
		getline(file, line);

		// Store variables
		file >> imax >> jmax >> x_max >> y_max >> t >> t_max >> dt_out >> c >> boundary ;
	}

	// Close file
	file.close();
	//------------------------------------------------------------------------------------------------------------|
	//============================================================================================================|
	// START TIMING FOR WHOLE PROCESS
	//============================================================================================================|
	auto time_total_start = chrono::high_resolution_clock::now();

	//------------------------------------------------------------------------------------------------------------|
	// Ceate a domain object
	auto* domain_instance = new Domain<double>();
	// Setup Timing Test object
	auto* Timer = new TimingTest;

	//------------------------------------------------------------------------------------------------------------|
	// Set global parameters to our domain
	domain_instance->global_imax = imax;
	domain_instance->global_jmax = jmax;

	//============================================================================================================|
	// START TIMING FOR DOMAIN DECOMPOSITION
	//============================================================================================================|
	 auto time_domain_start = chrono::high_resolution_clock::now();
	//============================================================================================================|
	
	// Find the dimeansions of the sub-domain
	domain_instance->createSubDomains(p);

	//------------------------------------------------------------------------------------------------------------|
	// Allocate the pointers for all the grids
	domain_instance->setSubDomainSize();

	//------------------------------------------------------------------------------------------------------------|
	// Allocate memory to the process domains
	domain_instance->allocatePointers();

	//------------------------------------------------------------------------------------------------------------|
	// Set Initial Disturbance
	setInitialDisturbance(*domain_instance);
	//fillSubDomainsWithProcessId(*domain_instance);


	//------------------------------------------------------------------------------------------------------------|
	// Find the neighbours of the process
	domain_instance->findNeighbours();

	//============================================================================================================|
	auto time_domain_end = std::chrono::high_resolution_clock::now();
	Timer->domain_setup = time_domain_end - time_domain_start;
	auto count_domain_setup = Timer->domain_setup.count();
	Timer->domain_setup_time.push_back(count_domain_setup);
	//cout << "Domain Setup Time for process " <<  id << " is " << count_domain_setup << "s " << endl;
	//Timer->timing_domain_setup(id, p);
	//============================================================================================================|

	

	//------------------------------------------------------------------------------------------------------------|
	// Calculate time stepping parameters
	dt = 0.1 * min(dx, dy) / c;

	//------------------------------------------------------------------------------------------------------------|
	// Set output count and iteration paramters
	int out_cnt = 0, it = 0;

	//------------------------------------------------------------------------------------------------------------|
	// Output Initial Disturbance to file and update the time counters
	grid_to_file(*domain_instance, out_cnt, boundary);
	out_cnt++;
	t_out += dt_out;


	//------------------------------------------------------------------------------------------------------------|
	//  MPI DATATYPES FOR COMMUNICATIONS
	// -----------------------------------------------------------------------------------------------------------|
	//============================================================================================================|
	// START TIMING FOR MPI DATATYPES
	//============================================================================================================|
	auto time_MPI_create_start = chrono::high_resolution_clock::now();
	//============================================================================================================|

	// Create instance of MPI_Data class to access MPI Datatypes and send and receive functions
	MPI_Data<double>* MData = new MPI_Data<double>();

	// Create MPI Data Types
	MData->CreateMPI_Types(*domain_instance);

	// Setup requests for non-blocking communications
	MPI_Request* request_list;
	request_list = nullptr;


	//============================================================================================================|
	auto time_MPI_create_end = std::chrono::high_resolution_clock::now();
	Timer->time_MPI_create = time_MPI_create_end - time_MPI_create_start;
	Timer->mpi_create_time.push_back(Timer->time_MPI_create.count());
	//cout << "MPI Create Time for process " <<  id << " is " << time_MPI_create.count() << "s " << endl;
	//Timer->timing_MPI_Create_setup(id, p);
	//============================================================================================================|



	while (t < t_max)
	{
		//---------------------------------------------------------------------------------------------------------|
		//		START MPI COMMUNICATIONS
		//---------------------------------------------------------------------------------------------------------|
	
		// Setup counter for requests
		int count = 0;

		// Initialize with 8 requests - 4 forr send and 4 for receive for each process
		request_list = new MPI_Request[8 * p];

		//============================================================================================================|
		// START TIMING FOR MPI COMMUNICATIONS
		//============================================================================================================|
		auto time_MPI_comms_start = chrono::high_resolution_clock::now();
		//============================================================================================================|
		//---------------------------------------------------------------------------------------------------------|
		// Do MPI Communications between the processes to send and receive boundary data
		doMPICommunications(*domain_instance, *MData, count, request_list);


		//============================================================================================================|
		auto time_MPI_comms_end = std::chrono::high_resolution_clock::now();
		Timer->time_MPI_comms = time_MPI_comms_end - time_MPI_comms_start;
		Timer->mpi_comms_time.push_back(Timer->time_MPI_comms.count());
		//cout << "MPI Ccommunication Time for process " <<  id << " is " << time_MPI_comms.count() << "s " << endl;
		//Timer->timing_MPI_comms(id, p);
		//============================================================================================================|


		//============================================================================================================|
		// START TIMING FOR  INNER ITERATIONS
		//============================================================================================================|
		auto time_inner_iterations_start = chrono::high_resolution_clock::now();

		// Do Iterations for the Inner Domain of the process grid
		doInnerIterations(*domain_instance);

		//============================================================================================================|
		auto time_inner_iterations_end = std::chrono::high_resolution_clock::now();
		Timer->time_inner_iterations = time_inner_iterations_end - time_inner_iterations_start;
		Timer->inner_iterations_time.push_back(Timer->time_inner_iterations.count());
		//cout << "Average inner iteration time for all process " <<  id << " is " << time_inner_iterations.count() << "s " << endl;
		//Timer->timing_inner_iterations(id, p);*/
		//============================================================================================================|


		// Wait for all MPI processes to finish communicating
		MPI_Waitall(count, request_list, MPI_STATUSES_IGNORE);



		//============================================================================================================|
		// START TIMING FOR  BOUNDARY SELECTION AND ITERATIONS
		//============================================================================================================|
		auto time_edge_iterations_start = chrono::high_resolution_clock::now();

		//------------------------------------------------------------------------------------------------------------|
		// If not periodic, set Boundary Conditions depening upon selected input(above)
		if(boundary != 2) // Used to avoid additional checks, however if its periodic, this is automatically negated
			selectBoundary(*domain_instance, boundary); // This is set to neumann, since boundary = 1 by default


		// Do Iterations on the Edgaes
		doBoundaryIterations(*domain_instance);
		
		//============================================================================================================|
		auto time_edge_iterations_end = std::chrono::high_resolution_clock::now();
		Timer->time_edge_iterations = time_edge_iterations_end - time_edge_iterations_start;
		Timer->edge_iterations_time.push_back(Timer->time_edge_iterations.count());
		//cout << "Edge iteration time for process " <<  id << " is " << time_edge_iterations.count() << "s " << endl;
		//Timer->timing_edge_iterations(id, p);
		//============================================================================================================|		

		//---------------------------------------------------------------------------------------------------------|
		// Swap pointers of old grid with new grid and new grid with current grid
		domain_instance->swapGridPointers(domain_instance->old_grid_2d, domain_instance->current_grid_2d, 
										  domain_instance->new_grid_2d, domain_instance->old_grid_1d,
										  domain_instance->current_grid_1d, domain_instance->new_grid_1d);
		//---------------------------------------------------------------------------------------------------------|
		// Increment time step
		t += dt;
		
		//---------------------------------------------------------------------------------------------------------|
		// Output the grid contents to files
		if (t_out <= t)
		{
			cout << "output: " << out_cnt << "\tt: " << t << "\titeration: " << it << endl;
			grid_to_file(*domain_instance, out_cnt, boundary);
			out_cnt++;
			t_out += dt_out;
		}
		
		//---------------------------------------------------------------------------------------------------------|
		// Incerement the number of iterations
		it++;
		
		//---------------------------------------------------------------------------------------------------------|
		// Delete pointer  the request list
		delete[] request_list;

	}


	//============================================================================================================|
	auto time_total_end = std::chrono::high_resolution_clock::now();
	Timer->time_total = time_total_end - time_total_start;
	Timer->total_time.push_back(Timer->time_total.count());

	//Timer->timing_total(id, p, imax, jmax);
	if (id == 0)
	{
		cout << " ---------------------------------------------------------------------------------" << endl;
		cout << " SIMULATION PARAMETERS  " << endl;
		cout << " ---------------------------------------------------------------------------------" << endl;
		cout << " Domain Size : " << imax << " X " << jmax << endl;
		cout << " Total number of processors " << p << endl;
		cout << " The domain divides into" << domain_instance->sub_rows << " and "  << domain_instance->sub_cols << " columns." << endl;
		cout << " Total Time : " << t_max << "s. Total Iterations : " << it << endl;
		cout << " Boundary Type ( 0 = Dirichlit, 1 = Neumann, 2 = Periodic) : " << boundary << endl;
		cout << "----------------------------------------------------------------------------------" << endl;
		cout << endl;
		cout << endl;
		cout << " ---------------------------------------------------------------------------------" << endl;
		cout << " TIMING ANALYSIS RESULTS " << endl;
		cout << " ---------------------------------------------------------------------------------" << endl;
		cout << " Time spech by a process to setup the domain : " << count_domain_setup << endl;
		cout << " Time spent by a process to create MPI Data Types : " << Timer->time_MPI_create.count() << endl;
		cout << " Time spent by a process to carry out MPI communications : " << Timer->time_MPI_comms.count() << endl;
		cout << " Time spent by a process to finish a set of inner iterations : " << Timer->time_inner_iterations.count() << endl;
		cout << " Time spent by a process to finish a set of edge iterations :  " << Timer->time_edge_iterations.count() << endl;
		cout << "Total time for a processor is : " << Timer->time_total.count() << "s. " << endl;
	}
	//============================================================================================================|		


	// End MPI Functions
	MPI_Finalize();
	return(0);
}
