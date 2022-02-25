
#pragma once
#include <mpi.h>
#include <vector>
#include "Domain.h"


template <class T>
class MPI_Data
{
public:
	//---------------------------------------------------------------------------------------------|
    //  CONSTRUCTORS AND DESTRUCTORS															   |
    // --------------------------------------------------------------------------------------------|
	// Constructor 
	MPI_Data();
	// Destructor to avoid explicitly deleting the memory
	~MPI_Data();


	//---------------------------------------------------------------------------------------------|
	//  MPI DATATYPES FOR COMMUNICATIONS														   |
	// --------------------------------------------------------------------------------------------|
	// MPI TYPES for the regualr rows
	MPI_Datatype Left_Type, Right_Type, Top_Type, Bottom_Type;
	// MPI Types for the Ghost rows
	MPI_Datatype Left_Ghost_Type, Right_Ghost_Type, Top_Ghost_Type, Bottom_Ghost_Type;


	//----------------------------------------------------------------------------------------------|
	//  CREATE MPI DATA TYPES																		|
	// ---------------------------------------------------------------------------------------------|
	// function to create all data types
	void CreateMPI_Types(Domain<T>& Dom);

	//----------------------------------------------------------------------------------------------|
	// FUNCTIONS FOR COMMUNICATIONS																	|
	// ---------------------------------------------------------------------------------------------|
	void send_data(Domain<T>& Dom, MPI_Datatype MDataType, int neigh_to_send, MPI_Request* request);
	void recv_data(Domain<T>& Dom, MPI_Datatype MDataType, int neigh_to_receive, MPI_Request* request);

};

template class MPI_Data<double>;
template class MPI_Data<int>;