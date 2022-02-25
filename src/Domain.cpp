
#define _USE_MATH_DEFINES

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <cmath>

#include "Domain.h"

using namespace std;


template <class T>
Domain<T>::Domain()
{
	//------------------------------------------------------------------------------|
	//			CREATE DOMAIN
	//------------------------------------------------------------------------------|
	// assign pointers to null
	// 1-D Arrays 
	T* old_grid_1d = nullptr;
	T* grid_1d = nullptr;
	T* new_grid_1d = nullptr;
	// 2-D Arrays
	T** old_2d = nullptr;
	T** grid_2d = nullptr;
	T** new_2d = nullptr;


	// Subdomain Divisions
	this->start_row = new int[p];
	this->num_rows = new int[p];
	this->start_column = new int[p];
	this->num_columns = new int[p];
}

template <class T>
Domain<T>::Domain(int n_i, int n_j)
{
	allocatePointers();
}

template <class T>
 Domain<T>::~Domain()
{
	// 1-D Arrays 
	delete[] old_grid_1d;
	delete[] current_grid_1d;
	delete[] new_grid_1d;

	// 2-D Arrays
	delete[] old_grid_2d;
	delete[] current_grid_2d;
	delete[] new_grid_2d;
}
 

 template <class T>
void Domain<T>::allocatePointers()
{
	// Allocate memory to the pointers and initialize to zero 
	old_grid_1d = new T[sub_imax * sub_jmax]{0};	 // set to total elements
	current_grid_1d = new T[sub_imax * sub_jmax]{0};
	new_grid_1d = new T[sub_imax * sub_jmax]{0};

	old_grid_2d = new T * [sub_imax];			// set to total number of rows
	current_grid_2d = new T * [sub_imax];
	new_grid_2d = new T * [sub_imax];



	// point the 2D array to the 1D array for the domain
	for (int i = 0; i < sub_imax; i++)
	{
		old_grid_2d[i] = &old_grid_1d[i * sub_jmax];			// set to total number of rows
		current_grid_2d[i] = &current_grid_1d[i * sub_jmax];
		new_grid_2d[i] = &new_grid_1d[i * sub_jmax];
	}
	
}


template <class T>
void Domain<T>::swapGridPointers(T**&subDomain2D_1, T**&subDomain2D_2, T**& subDomain2D_3, T * & subDomain1D_1, T*& subDomain1D_2, T*& subDomain1D_3)
{
	// swap the grid pointers of a 2-D grid using a swap variable
	T** swap1 = subDomain2D_1;
	subDomain2D_1 = subDomain2D_2;
	subDomain2D_2 = subDomain2D_3;
	subDomain2D_3 = swap1;

	// swap the grid pointers of a 1-D grid using a swap variable
	T* swap2 = subDomain1D_1;
	subDomain1D_1 = subDomain1D_2;
	subDomain1D_2 = subDomain1D_3;
	subDomain1D_3 = swap2;

}




//------------------------------------------------------------------------------------------------------
//								FIND SUBDOMAINS 
//------------------------------------------------------------------------------------------------------
template <class T>
void Domain<T>::createSubDomains(int p)
{
	int min_gap = p;
	int top = sqrt(p) + 1;
	for (int i = 1; i <= top; i++)
	{
		if (p % i == 0)
		{
			int gap = abs(p / i - i);
			if (gap < min_gap)
			{
				min_gap = gap;
				this->sub_rows = i;
				this->sub_cols = p / i;
			}
		}
		
	}


}

template <class T>
void Domain<T>::setSubDomainSize()
{
	// default processor variables
	num_rows = new int[sub_rows]{};
	start_row = new int[sub_rows]{};
	num_columns = new int[sub_cols]{};
	start_column = new int[sub_cols]{};

	int start_r = 0;
	for (int k = 0; k < this->sub_rows; k++)
	{
		// Dividing the rows among the processes
		num_rows[k] = (this->global_jmax - start_r) / (this->sub_rows-k);
		start_row[k] = start_r; // set the first start point as 0

		// incrementing start point to the next start
		start_r += num_rows[k];
	}


	int start_c = 0;
	for (int k = 0; k < this->sub_cols; k++)
	{
		// Dividing the columns among the processes
		num_columns[k] =(this->global_imax - start_c) / (this->sub_cols-k);
		start_column[k] = start_c;

		start_c += num_columns[k];
	}
	
	// With padding the maximum value of the subdomain is incremented by 2
	this->sub_imax = this->num_rows[id / sub_cols] + 2;
	this->sub_jmax = this->num_columns[id % sub_cols] + 2;

	// Set subdomain variables to processor indices
	this->sub_i = id / sub_cols;
	this->sub_j = id % sub_cols;


	
	// Set Global domain parameters 
	this->global_i = start_row[id / sub_cols];
	this->global_j = start_column[id % sub_cols];


}


template <class T>
void Domain<T>::findNeighbours()
{

	// set the neighbour array to -1
	for(int i = 0; i < 4; i++)
	{
		neighbours[i] = -1;
	}

	
		// Loop on left and right of the current processor position
		for (int i = -1; i <= 1; i++)
		{
			// Loop on top and bottom of the current processor position
			for (int j = -1; j <= 1; j++)
			{
				//------------------------------------------------------------------------------|
				// PERIODIC NEIGHBOURS
				//------------------------------------------------------------------------------|
				if (boundary == 2) // If the boubdary is periodic
				{
					// so we dont go negative we add sub_cols and sub_rows to the numbers
					int neigh_i = (this->sub_j + i + this->sub_cols) % this->sub_cols;
					int neigh_j = (this->sub_i + j + this->sub_rows) % this->sub_rows;
					int neigh_id = neigh_i + neigh_j * this->sub_cols;

					if (neigh_id != id)
					{
						if (i == -1 && j == 0) // left neighbour
						{
							this->neighbours[0] = neigh_id;
						}
						else if (i == 0 && j == -1) // top neighbour
						{
							this->neighbours[1] = neigh_id;
						}
						else if (i == 1 && j == 0) // right neighbour
						{
							this->neighbours[2] = neigh_id;
						}
						else if (i == 0 && j == 1) // bottom neighbour
						{
							this->neighbours[3] = neigh_id;
						}
					}
				}
				else
				{
					//------------------------------------------------------------------------------|
					// NON PERIODIC NEIGHBOURS
					//------------------------------------------------------------------------------|
					// Find potential neighbours
					int neigh_i = this->sub_j + i; // Find  left and right neighbours
					int neigh_j = this->sub_i + j; // Find  top and bottom neighbours

					// If found neighbours exist (ie they are non-negative), find the neightbour id 
					if (neigh_i >= 0 && neigh_i < this->sub_cols && neigh_j >= 0 && neigh_j < this->sub_rows)
					{
						int neigh_id = neigh_i + neigh_j * this->sub_cols;
						if (neigh_id != id)
						{
							if (i == -1 && j == 0) // left neighbour
							{
								this->neighbours[0] = neigh_id;
							}
							else if (i == 0 && j == -1) // top neighbour
							{
								this->neighbours[1] = neigh_id;
							}
							else if (i == 1 && j == 0) // right neighbour
							{
								this->neighbours[2] = neigh_id;
							}
							else if (i == 0 && j == 1) // bottom neighbour
							{
								this->neighbours[3] = neigh_id;
							}

						}
					}
				}
			}
		}

}


template <class T>
void Domain<T>:: id_to_index(int id, int& id_row, int& id_column) 
{
	id_column = id % sub_cols;
	id_row = id / sub_cols;
}

template <class T>
int Domain<T>:: id_from_index(int id_row, int id_column)
{
	if (id_row >= sub_rows || id_row < 0)
		return -1;
	if (id_column >= sub_cols || id_column < 0)
		return -1;

	return id_row * sub_cols + id_column;
}


template <class T>
void Domain<T>::printsubDomain(Domain<T>* Dom)
{
	cout << " Dom->sub_imax" << Dom->sub_imax << " Dom->sub_jmax "<< Dom->sub_jmax << endl;
	cout << endl;
	cout << " CURRENT GRID" << endl;
	for (int i = 0; i < Dom->sub_imax; i++)
	{
		for (int j = 0; j < Dom->sub_jmax; j++)
		{
			cout << Dom->current_grid_2d[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	cout << endl;
	cout << " NEW GRID" << endl;
	for (int i = 0; i < Dom->sub_imax; i++)
	{
		for (int j = 0; j < Dom->sub_jmax; j++)
		{
			cout << Dom->new_grid_2d[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;

	cout << " OLD GRID" << endl;
	for (int i = 0; i < Dom->sub_imax; i++)
	{
		for (int j = 0; j < Dom->sub_jmax; j++)
		{
			cout << Dom->new_grid_2d[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


