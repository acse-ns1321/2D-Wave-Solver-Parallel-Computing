#include <vector>
#include<iostream>

#include "TimingTest.h"

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------|

void TimingTest::timing_domain_setup(int id, int p)
{
	if (id == 0)
	{
		double avg_domain_time = 0;
		for (int i  = 0; i < domain_setup_time.size(); i++)
		{
			avg_domain_time += domain_setup_time[i];
		}
		avg_domain_time = avg_domain_time / p;

		cout << "The Average Domain Setup Time for all the processes is " << avg_domain_time << "s. " << endl;
	}
}

// -------------------------------------------------------------------------------------------------------------------------|

void TimingTest::timing_MPI_Create_setup(int id, int p)
{
	if (id == 0)
	{
		double avg_MPI_Create_time = 0;
		for (int i = 0; i < mpi_create_time.size(); i++)
		{
			avg_MPI_Create_time += mpi_create_time[i];
		}
		avg_MPI_Create_time = avg_MPI_Create_time / p;

		cout << "The Average MPI Create Time for all the processes is " << avg_MPI_Create_time << "s. " << endl;
	}

}

// -------------------------------------------------------------------------------------------------------------------------|

void TimingTest::timing_MPI_comms(int id, int p)
{
	if (id == 0)
	{
		double avg_MPI_comms_time = 0;
		for (int i = 0; i < mpi_comms_time.size(); i++)
		{
			avg_MPI_comms_time += mpi_comms_time[i];
		}
		avg_MPI_comms_time = avg_MPI_comms_time / p;

		cout << "The Average MPI Communication Time for all the processes is " << avg_MPI_comms_time << "s. " << endl;
	}

}

// -------------------------------------------------------------------------------------------------------------------------|

void TimingTest::timing_edge_iterations(int id, int p)
{
	if (id == 0)
	{
		double avg = 0;
		for (int i = 0; i < edge_iterations_time.size(); i++)
		{
			avg += edge_iterations_time[i];
		}
		avg = avg / p;

		cout << "The Average Edge Iteration Time for all the processes is " << avg << "s. " << endl;
	}

}

// -------------------------------------------------------------------------------------------------------------------------|

void TimingTest::timing_inner_iterations(int id, int p)
{
	if (id == 0)
	{
		double avg = 0;
		for (int i = 0; i < inner_iterations_time.size(); i++)
		{
			avg += inner_iterations_time[i];
		}
		avg = avg / p;

		cout << "The Average Inner Iteration for all the processes is " << avg << "s. " << endl;
	}

}

// -------------------------------------------------------------------------------------------------------------------------|

void TimingTest::timing_total(int id, int p, int imax, int jmax)
{
	if (id == 0)
	{
		double avg = 0;
		for (int i = 0; i < this->total_time.size(); i++)
		{
			avg += total_time[i];
		}
		avg = avg / p;
		cout << "The Average Total time for all the processes is " << avg << "s. " << endl;
	}

}

// -------------------------------------------------------------------------------------------------------------------------|