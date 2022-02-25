
#define _USE_MATH_DEFINES
#include<iostream>
#include<string>
#include<mpi.h>
#include <sstream>
#include <fstream>
#include <cmath>


#include "WaveSolver.h"
#include "Parameters.h"

using namespace std;

void readWaveSolverParameters(int &i , int &j , double &x_m, double &y_m, double &tee, double &t_m,double &dt_o , double &cee , int &bo,double &dist, double *disturbance_p)
{

	stringstream data;
	data << "Parameters.txt";
	ifstream file;

	// Open Parameters file
	file.open(data.str().c_str(), ios_base::in);

	// Read in value by value
	if (file.good())
	{
		file >> i >> j >> x_m >> y_m >> t >> t_m >> dt_o >> c >> bo >> dist;
	}

	// Initialize distrbance array with 3 parameters(radius, x loc, y loc)
	disturbance_p = new double[dist * 3];

	while (file.good())
	{
		// Fill te disturbance array with the 3 parameters
		for (int i = 0; i < 3 * dist; i++)
		{
			file >> disturbance_p[i];
		}
	}
	file.close();
	cout << i << j;
}

// ------------------------------------------------------------------------------------------------|
// Select Boundary
void selectBoundary(Domain<double> &Dom, int boundary)
{
	if (boundary == 0)
	{
		dirichlitBoundary(Dom);
	}
	if (boundary == 1)
	{
		neumannBouandry(Dom);
	}

}

// Set Dirichlit Boundry Conditions
void dirichlitBoundary(Domain<double>& Dom)
{
	// If the left neighbours dont exist
	if (Dom.neighbours[0] == -1)
	{
		// Set left boundaries
		for (int i = 1; i < Dom.sub_imax - 1; i++)
		{
			Dom.current_grid_2d[i][0] = 0;
		}
	}

	// If the right neighbours dont exist
	if (Dom.neighbours[2] == -1)
	{
		// Set right boundaries 
		for (int i = 1; i < Dom.sub_imax - 1; i++)
		{
			Dom.current_grid_2d[i][Dom.sub_jmax - 1] = 0;
		}
	}
	
	//If the top neighbours dont exist
	if (Dom.neighbours[1] == -1)
	{
		// Set top boundaries 
		for (int i = 1; i < Dom.sub_jmax - 1; i++)
		{
			Dom.current_grid_2d[1][i] = 0;
		}
	}

	// If the bottom neighbours dont exist
	if (Dom.neighbours[3] == -1)
	{
		// Set bottom boundaries
		for (int i = 1; i < Dom.sub_jmax - 1; i++)
		{
			Dom.current_grid_2d[Dom.sub_imax - 1][i] = 0;
		}
	}
}

// Set Neumann Boundary Conditions
void neumannBouandry(Domain<double>& Dom)
{
	
	// If the left neighbours dont exist
	if (Dom.neighbours[0] == -1)
	{
		// Set left boundaries
		for (int i = 1; i < Dom.sub_imax - 1; i++)
		{
			Dom.current_grid_2d[i][1] = Dom.current_grid_2d[i][2];
		}
	}

	// If the right neighbours dont exist
	if (Dom.neighbours[2] == -1)
	{
		// Set right boundaries
		for (int i = 1; i < Dom.sub_imax - 1; i++)
		{
			Dom.current_grid_2d[i][Dom.sub_jmax - 1] = Dom.current_grid_2d[i][Dom.sub_jmax - 2];
		}
	}

	//If the top neighbours dont exist
	if (Dom.neighbours[1] == -1)
	{
		// Set top boundaries 
		for (int i = 1; i < Dom.sub_jmax - 1; i++)
		{
			Dom.current_grid_2d[1][i] = Dom.current_grid_2d[2][i];
		}
	}

	// If the bottom neighbours dont exist
	if (Dom.neighbours[3] == -1)
	{
		// Set bottom boundaries
		for (int i = 1; i < Dom.sub_jmax - 1; i++)
		{
			Dom.current_grid_2d[Dom.sub_imax-1][i] = Dom.current_grid_2d[Dom.sub_imax - 2][i];
		}
	}
}


// ------------------------------------------------------------------------------------------------|
// Create Disturbance
void setInitialDisturbance(Domain<double>& Dom)
{

		// Create an array to store x,y and radius for each disturbance
		// Calculate the grid spacing
		dx = x_max / (Dom.global_imax - 1);
		dy = y_max / (Dom.global_jmax - 1);
		
		double x_splash = 3;
		double y_splash = 3;
		double r_splash = 1;

		// Start looping on the non-padded region
		for (int i = 1; i < Dom.sub_imax - 1; i++)
		{
			for (int j = 1; j < Dom.sub_jmax - 1; j++)
			{

				// map global to local
				int I = (i - 1) + Dom.global_i;
				int J = (j - 1) + Dom.global_j;

				// set the grid jumps
				double x = dx * I;
				double y = dy * J;

				// find the distance of a point from the splash
				double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));
				double h = 1;

				// if the distance is within splash range, set the vlaues to calculated h
				if (dist < r_splash)
				{
					// calculate wave height
					h = 5.0 * (cos(dist / r_splash * M_PI) + 1.0);

					// set the grid values
					Dom.current_grid_2d[i][j] = h;
					Dom.old_grid_2d[i][j] = h;
				}

			}
		}

}


void fillSubDomainsWithProcessId(Domain<double>& Dom)
{
	for (int i = 0; i < Dom.sub_imax; i++)
	{
		for (int j = 0; j < Dom.sub_jmax; j++)
		{
			Dom.current_grid_2d[i][j] = id;
		}
	}
}




// ------------------------------------------------------------------------------------------------|
// Function to do the communications
void doMPICommunications(Domain<double> &domain_instance,  MPI_Data<double> &MData, int count, MPI_Request *request_list)
{

	//--------------------------------------------------------------------------------------------|
	//		SEND FROM REGULAR BOUNDARIES AND RECEIVE INTO GHOST BOUNDARIES
	// -------------------------------------------------------------------------------------------|
	// Send the regualar cells of the sub-domain to the other processes

	if (domain_instance.neighbours[0] != -1)
	{
		// Receive left data in ghost cells
		MData.recv_data(domain_instance, MData.Left_Ghost_Type, domain_instance.neighbours[0], &request_list[count]);
		count++;

		// Send left Data
		MData.send_data(domain_instance, MData.Left_Type, domain_instance.neighbours[0], &request_list[count]);
		count++;
	}


	if (domain_instance.neighbours[2] != -1)
	{
		// Receive right data in ghost cells
		MData.recv_data(domain_instance, MData.Right_Ghost_Type, domain_instance.neighbours[2], &request_list[count]);
		count++;

		// Send right data
		MData.send_data(domain_instance, MData.Right_Type, domain_instance.neighbours[2], &request_list[count]);
		count++;
	}

		if (domain_instance.neighbours[1] != -1)
		{
			// Receive top data in ghost cells
			MData.recv_data(domain_instance, MData.Top_Ghost_Type, domain_instance.neighbours[1], &request_list[count]);
			count++;

			// Send top data
			MData.send_data(domain_instance, MData.Top_Type, domain_instance.neighbours[1], &request_list[count]);
			count++;
		}


		if (domain_instance.neighbours[3] != -1)
		{
			// Receive bottom data in ghost cells
			MData.recv_data(domain_instance, MData.Bottom_Ghost_Type, domain_instance.neighbours[3], &request_list[count]);
			count++;

			// Send bottom data
			MData.send_data(domain_instance, MData.Bottom_Type, domain_instance.neighbours[3], &request_list[count]);
			count++;
		}

	if(boundary == 2) 
	{		
			//--------------------------------------------------------------------------------------------------------------|
			//  ADDITIONAL PERIODIC COMMUNICATIONS ON THE TOP AND BOTTOM
			//--------------------------------------------------------------------------------------------------------------|
			// Send top data
			MData.send_data(domain_instance, MData.Top_Type, domain_instance.neighbours[1], &request_list[count]);
			count++;

			// Receive bottom data in ghost cells
			MData.recv_data(domain_instance, MData.Bottom_Ghost_Type, domain_instance.neighbours[1], &request_list[count]);
			count++;

			//--------------------------------------------------------------------------------------------------------------|

			// Receive top data in ghost cells
			MData.recv_data(domain_instance, MData.Top_Ghost_Type, domain_instance.neighbours[1], &request_list[count]);
			count++;

			// Send bottom data
			MData.send_data(domain_instance, MData.Bottom_Type, domain_instance.neighbours[3], &request_list[count]);
			count++;

	}
	MPI_Waitall(count, request_list, MPI_STATUSES_IGNORE);

}


// ------------------------------------------------------------------------------------------------|
// Solve inner domain
void doInnerIterations(Domain<double> &Dom) 
{
	// Iterate over the innder domain in both i and j directions
	for (int i = 2; i < Dom.sub_imax - 1; i++)
	{
		for (int j = 2; j < Dom.sub_jmax - 1; j++)
		{
			// Descretization of equation of the wave for time stepping in the inner domain
			Dom.new_grid_2d[i][j] = pow(dt * c, 2.0) 
				* ((Dom.current_grid_2d[i + 1][j] 
						- 2.0 * Dom.current_grid_2d[i][j] 
						+ Dom.current_grid_2d[i - 1][j]) / pow(dx, 2.0)

				+ (Dom.current_grid_2d[i][j + 1] 
					- 2.0 * Dom.current_grid_2d[i][j] 
					+ Dom.current_grid_2d[i][j - 1]) / pow(dy, 2.0)) 

				+ 2.0 * Dom.current_grid_2d[i][j] - Dom.old_grid_2d[i][j];
		}
	}
}

// Solve outer domain
void doBoundaryIterations(Domain<double>& Dom)
{
	// LEFT BOUNDARY
	if (Dom.neighbours[0] != -1)
	{
		for (int i = 1; i < Dom.sub_imax-1; i++)
		{
			Dom.new_grid_2d[i][1] = pow(dt * c, 2.0) * 
				((Dom.current_grid_2d[i + 1][1] 
					- 2.0 * Dom.current_grid_2d[i][1] 
					+ Dom.current_grid_2d[i - 1][1]) / pow(dx, 2.0)

					+ (Dom.current_grid_2d[i][2] 
						- 2.0 * Dom.current_grid_2d[i][1] 
						+ Dom.current_grid_2d[i][0]) / pow(dy, 2.0)) 

				+ 2.0 * Dom.current_grid_2d[i][1] - Dom.old_grid_2d[i][1];
		}
	}

	// TOP BOUNDARY
	if (Dom.neighbours[1] != -1)
	{
		for (int j = 1; j < Dom.sub_jmax -1; j++)
		{
			Dom.new_grid_2d[1][j] = pow(dt * c, 2.0) 
				* ((Dom.current_grid_2d[2][j] 
					- 2.0 * Dom.current_grid_2d[1][j] 
					+ Dom.current_grid_2d[0][j]) / pow(dx, 2.0)

				+ (Dom.current_grid_2d[1][j + 1] 
					- 2.0 * Dom.current_grid_2d[1][j] 
					+ Dom.current_grid_2d[1][j - 1]) / pow(dy, 2.0)) 

				+ 2.0 * Dom.current_grid_2d[1][j] - Dom.old_grid_2d[1][j];
		}
	}

	// RIGHT BOUNDARY
	if (Dom.neighbours[2] != -1)
	{
		for (int i = 1; i < Dom.sub_imax - 1; i++)
		{
			Dom.new_grid_2d[i][Dom.sub_jmax-1] = pow(dt * c, 2.0) 
				* ((Dom.current_grid_2d[i + 1][Dom.sub_jmax-1] 
					- 2.0 * Dom.current_grid_2d[i][Dom.sub_jmax - 1]
					+ Dom.current_grid_2d[i - 1][Dom.sub_jmax - 1]) / pow(dx, 2.0) 

				+ (Dom.current_grid_2d[i][Dom.sub_jmax - 1] 
					- 2.0 * Dom.current_grid_2d[i][Dom.sub_jmax - 1]
					+ Dom.current_grid_2d[i][Dom.sub_jmax - 2]) / pow(dy, 2.0)) 

				+ 2.0 * Dom.current_grid_2d[i][Dom.sub_jmax - 1] - Dom.old_grid_2d[i][Dom.sub_jmax - 1];
		}

	}

	// BOTTOM BOUNDARY
	if (Dom.neighbours[3] != -1)
	{
		for (int j = 1; j < Dom.sub_jmax-1; j++)
		{
			Dom.new_grid_2d[Dom.sub_imax-1][j] = pow(dt * c, 2.0) 
				* ((Dom.current_grid_2d[Dom.sub_imax - 1][j] 
					- 2.0 * Dom.current_grid_2d[Dom.sub_imax - 1][j]
					+ Dom.current_grid_2d[Dom.sub_imax - 2][j]) / pow(dx, 2.0) 

				+ (Dom.current_grid_2d[Dom.sub_imax - 1][j + 1]
					- 2.0 * Dom.current_grid_2d[Dom.sub_imax - 1][j]
					+ Dom.current_grid_2d[Dom.sub_imax - 1][j - 1]) / pow(dy, 2.0)) 

				+ 2.0 * Dom.current_grid_2d[Dom.sub_imax - 1][j] - Dom.old_grid_2d[Dom.sub_imax - 1][j];
		}
	}

}


// Output files after iterations for postprocessing
void grid_to_file(Domain<double>& Dom, int out, int boundary)
{

	// Function that uses fstream to write variable to the file
	stringstream fname;
	fstream f1;

	// Outputting into folders according to the boundary condition and number
	if (boundary == 1)
	{
		fname << "./output/neumann/"<< p <<"_processors/"<<"output_" << out << "_id_" << id  << ".dat";
	}
	else if (boundary == 0)
	{
		fname << "./output/dirichlit/" << p << "_processors/" << "output_" << out << "_id_" << id << ".dat";
	}
	else if (boundary == 2)
	{
		fname << "./output/periodic/" << p << "_processors/" << "output_" << out << "_id_" << id << ".dat";
	}

	// Open and write contents of each processor to the respective files
	f1.open(fname.str().c_str(), ios_base::out);

	for (int i = 1; i < Dom.sub_imax - 1; i++)
	{
		// Writing contents of processsor grid to file
		for (int j = 1; j < Dom.sub_jmax - 1; j++)
			f1 << Dom.current_grid_2d[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
}