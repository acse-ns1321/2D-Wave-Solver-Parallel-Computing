#pragma once
#include <vector>
#include <chrono>

using namespace std;

class TimingTest
{
public:
	// Vector to store time spent in setting up the domain for each process and creating MPI_DataTypes
	vector<double> domain_setup_time;

	// Vector to store total run time for all iterations
	vector<double> total_time;

	// Vector to store timings for a single iteration
	vector<double> mpi_create_time;

	// Timing to doMPICommunications
	vector<double> mpi_comms_time;

	// Timing to complete innerIterations
	vector<double> inner_iterations_time;


	// Timings to complete set boundary conditions and complete outerIterations
	vector<double> edge_iterations_time;

	chrono::duration < double> time_MPI_create;
	chrono::duration <double> time_total;
	chrono::duration < double> domain_setup;
	chrono::duration <double> time_MPI_comms;
	chrono::duration < double> time_inner_iterations;
	chrono::duration < double> time_edge_iterations;

	//--------------------------------------------------------------------------------|
	//  FUNCTION DEALING TIMING DOMAIN SETUP 
	// -------------------------------------------------------------------------------|
	void timing_domain_setup(int id, int p);


	//--------------------------------------------------------------------------------|
	//  FUNCTIONS DEALING WITH MPI CREATION AND COMMUNICATIONS 
	// -------------------------------------------------------------------------------|
	void timing_MPI_Create_setup(int id, int p);

	void timing_MPI_comms(int id, int p);

	//--------------------------------------------------------------------------------|
	//  FUNCTION DEALING TIMING TIME STEP ITERATIONS
	// -------------------------------------------------------------------------------|

	void timing_inner_iterations(int id, int p);

	void timing_edge_iterations(int id, int p);
	//--------------------------------------------------------------------------------|
	//  FUNCTION DEALING TIMING TOTAL TIME FOR ALL PROCESSES 
	// -------------------------------------------------------------------------------|

	void timing_total(int id, int p, int imax, int jmax);

};
