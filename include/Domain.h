#pragma once
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <vector>

#include "Parameters.h"


template <class T>
class Domain
{
public:

	//----------------------------------------------------------------------------------------------|
	// CONSTRUCTORS AND DESTRUCTORS																	|
	// ---------------------------------------------------------------------------------------------|
	// constructor that assigns pointers
	Domain();
	// constructor that allocates data
	Domain(int n_i, int n_j);
	// destructor that prevents memory leaks in the matrix
	~Domain();


	//----------------------------------------------------------------------------------------------|
	// LOCAL SUBDOMAIN GRID VARIABLES																|
	// ---------------------------------------------------------------------------------------------|
	// Class variables that are pointers to 1-D and 2_D arrays of 3 required sub domains
	// 1-D Arrays 
	T* old_grid_1d;
	T* current_grid_1d;
	T* new_grid_1d;
	// 2-D Arrays
	T** old_grid_2d;
	T** current_grid_2d;
	T** new_grid_2d;

	// Sub Domain Parameters------------------------------------------------------------------------|
	// Variables that define maximum size of sub domain
	int sub_imax;
	int sub_jmax;
	// Varibales that define the indexing on the sub domain
	int sub_i;
	int sub_j;
	// Parametris that define the subdivision of the grid
	int sub_rows;
	int sub_cols;


	//----------------------------------------------------------------------------------------------|
	// GLOBAL SUBDOMAIN GRID VARIABLES																|
	// ---------------------------------------------------------------------------------------------|
	// Global Domain Parameters  --------------------------------
	// Variables that define maximum size of global domain
	int global_imax;
	int global_jmax;
	// Varibales that define the indexing on the global domain
	int global_i;
	int global_j;

	// Allocate the number of start rows and num-rows 
	int* start_row;
	int* num_rows ;
	int* start_column ;
	int* num_columns;

	//----------------------------------------------------------------------------------------------|
	// GRID FUNCTIONS TO CREATE AND SET DOMAIN PARAMERETS											|
	// ---------------------------------------------------------------------------------------------|

	// class method to allocate the data in the domain and sub domains
	void allocatePointers();

	// Create Sub Domains
	void createSubDomains(int p);

	// Class function that sets the subdomain sizes
	void setSubDomainSize();

	//----------------------------------------------------------------------------------------------|
	// GLOBAL SUBDOMAIN GRID VARIABLES																|
	// ---------------------------------------------------------------------------------------------|
	// Array to store the neighbours for this process
	int neighbours[4];

	// Class function that finds neighbours of the processes
	void findNeighbours();

	// Class functions to convert indices from local to global domains and vice-versa
	void id_to_index(int id, int& id_row, int& id_column);
	int id_from_index(int id_row, int id_column);

	// Function to swap pointers of grid - used after each time step
	void swapGridPointers(T**& subDomain2D_1, T**& subDomain2D_2, T**& subDomain2D_3, T*& subDomain1D_1, T*& subDomain1D_2, T*& subDomain1D_3);

	// Functions to print the indiviual process's domain
	void printsubDomain(Domain<T>* Dom);
};

template class Domain<double>;
template class Domain<int>;