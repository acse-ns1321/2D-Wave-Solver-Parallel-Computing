
#include <vector>
#include <mpi.h>

#include "MPIDatatype.h"

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------|

template <class T>
MPI_Data<T>::MPI_Data() 
{
	// Set MPI Types for the  cells to null
	MPI_Datatype Left_Type = MPI_DATATYPE_NULL;
	MPI_Datatype Right_Type = MPI_DATATYPE_NULL;
	MPI_Datatype Top_Type = MPI_DATATYPE_NULL;
	MPI_Datatype Bottom_Type = MPI_DATATYPE_NULL;

	// Set MPI Types for the Ghost cells to null
	MPI_Datatype Left_Ghost_Type = MPI_DATATYPE_NULL; 
	MPI_Datatype Right_Ghost_Type = MPI_DATATYPE_NULL;
	MPI_Datatype Top_Ghost_Type = MPI_DATATYPE_NULL;
	MPI_Datatype Bottom_Ghost_Type = MPI_DATATYPE_NULL;
	
}


template <class T>
MPI_Data<T>::~MPI_Data()
{
	// deleting regular datatypes
	MPI_Type_free(&Left_Type);
	MPI_Type_free(&Right_Type);
	MPI_Type_free(&Top_Type);
	MPI_Type_free(&Bottom_Type);

	// deleting ghost datatypes
	MPI_Type_free(&Left_Ghost_Type);
	MPI_Type_free(&Right_Ghost_Type);
	MPI_Type_free(&Top_Ghost_Type);
	MPI_Type_free(&Bottom_Ghost_Type);
}

// -------------------------------------------------------------------------------------------------------------------------|

template <class T>
void MPI_Data<T>::send_data(Domain<T> &Dom, MPI_Datatype MDataType, int send_to, MPI_Request* request)
{
	// Function that facilitates the ease of sends the MPI Data Type
	MPI_Isend(Dom.current_grid_1d, 1, MDataType, send_to, tag_num, MPI_COMM_WORLD, request);
}


template <class T>
void MPI_Data<T>::recv_data(Domain<T>& Dom, MPI_Datatype MDataType, int receive_from, MPI_Request* request)
{
	// Function that facilitates the ease of receives the MPI Data Type
	MPI_Irecv(Dom.current_grid_1d, 1, MDataType, receive_from, tag_num, MPI_COMM_WORLD, request);
}


// -------------------------------------------------------------------------------------------------------------------------|

template <class T>
void MPI_Data<T>::CreateMPI_Types(Domain<T>& Dom)
{

	// Use vectors because we dont explicitly know the lengths
	vector<int> block_lengths;
	vector<MPI_Aint> displacements;
	MPI_Aint add_start;
	vector<MPI_Datatype> typelist;


	// Resize all lengths of vector to Dom.sub_imax-2
	block_lengths.resize(Dom.sub_imax - 2);
	displacements.resize(Dom.sub_imax - 2);
	typelist.resize(Dom.sub_imax - 2);

	// Get the address and set the object to the start of the address
	// Explanation: This is done so that we can send the data by referencing the object
	MPI_Get_address(Dom.current_grid_1d, &add_start);



	// ------------------------------------------------------------------------------------------------------------------------|
	//  MPI_DATATYPE 1 : LEFT_TYPE MPI - Regular
	// ------------------------------------------------------------------------------------------------------------------------|
	// Set each of the dataType variables
	// to the correct addresses depending on the location
	// in the array

	// Initialize the variables
	for (int i = 0; i < Dom.sub_imax - 2; i++)
	{
		typelist[i] = MPI_DOUBLE;
		block_lengths[i] = 1;
		MPI_Get_address(&Dom.current_grid_2d[i + 1][1], &displacements[i]);
		displacements[i] = displacements[i] - add_start;

	}
	// Create the dataType
	MPI_Type_create_struct((Dom.sub_imax - 2), block_lengths.data(), displacements.data(), typelist.data(), &Left_Type);
	// Commit to use
	MPI_Type_commit(&Left_Type);



	// ------------------------------------------------------------------------------------------------------------------------|
	//  MPI_DATATYPE 2 : LEFT_TYPE MPI - GHOST
	// ------------------------------------------------------------------------------------------------------------------------|
	for (int i = 0; i < Dom.sub_imax - 2; i++)
	{
		typelist[i] = MPI_DOUBLE;
		block_lengths[i] = 1;
		MPI_Get_address(&Dom.current_grid_2d[i + 1][0], &displacements[i]);
		displacements[i] = displacements[i] - add_start;
	}
	// Create the dataType
	MPI_Type_create_struct(Dom.sub_imax - 2, block_lengths.data(), displacements.data(), typelist.data(), &Left_Ghost_Type);
	// Commit to use
	MPI_Type_commit(&Left_Ghost_Type);



	// ------------------------------------------------------------------------------------------------------------------------|
	//  MPI_DATATYPE 3 :RIGHT_TYPE MPI - Regular
	// ------------------------------------------------------------------------------------------------------------------------|
	// The only difference is that the address to store in the 
	// 2D array is [i][Dom.sub_jmax-1]

	// In+ialize the variables for all the grids
	for (int i = 0; i < Dom.sub_imax - 2; i++)
	{
		typelist[i] = MPI_DOUBLE;
		block_lengths[i] = 1;
		MPI_Get_address(&Dom.current_grid_2d[i + 1][Dom.sub_jmax - 2], &displacements[i]);
		displacements[i] = displacements[i] - add_start;
	}
	// Create the dataType
	// Note: We use .data() because they are vectors
	MPI_Type_create_struct((Dom.sub_imax - 2), block_lengths.data(), displacements.data(), typelist.data(), &Right_Type);
	// Commit to use
	MPI_Type_commit(&Right_Type);



	// ------------------------------------------------------------------------------------------------------------------------|
	//  MPI_DATATYPE 4 :RIGHT_TYPE MPI - Ghost
	// ------------------------------------------------------------------------------------------------------------------------|
	// Initialize the variables for the ghost right grids
	for (int i = 0; i < Dom.sub_imax - 2; i++)
	{
		typelist[i] = MPI_DOUBLE;
		block_lengths[i] = 1;
		MPI_Get_address(&Dom.current_grid_2d[i + 1][Dom.sub_jmax - 1], &displacements[i]);
		displacements[i] = displacements[i] - add_start;
	}
	// Create the dataType
	// Note: We use .data() because they are vectors
	MPI_Type_create_struct((Dom.sub_imax - 2), block_lengths.data(), displacements.data(), typelist.data(), &Right_Ghost_Type);
	// Commit to use
	MPI_Type_commit(&Right_Ghost_Type);



	// ------------------------------------------------------------------------------------------------------------------------|
	//  TOP AND BOTTOM
	// ------------------------------------------------------------------------------------------------------------------------|
	// Note: the top and bottom are stored in contiguous memory so dont need 
	// to be sent as different types
	// ------------------------------------------------------------------------------------------------------------------------|



	// ------------------------------------------------------------------------------------------------------------------------|
	// MPI_DATATYPE 5 : TOP - REGULAR
	// ------------------------------------------------------------------------------------------------------------------------|

	// Initialize variables
	typelist[0] = MPI_DOUBLE;
	block_lengths[0] = Dom.sub_jmax - 2;
	MPI_Get_address(&Dom.current_grid_2d[1][1], &displacements[0]);
	displacements[0] = displacements[0] - add_start;
	// Create the dataType
	MPI_Type_create_struct(1, block_lengths.data(), displacements.data(), typelist.data(), &Top_Type);
	// Commit to use
	MPI_Type_commit(&Top_Type);


	// ------------------------------------------------------------------------------------------------------------------------|
	// MPI_DATATYPE 6 : TOP - GHOST
	// ------------------------------------------------------------------------------------------------------------------------|
	// Initialize variables
	typelist[0] = MPI_DOUBLE;
	block_lengths[0] = Dom.sub_jmax - 2;
	MPI_Get_address(&Dom.current_grid_2d[0][1], &displacements[0]);
	displacements[0] = displacements[0] - add_start;
	// Create the dataType
	MPI_Type_create_struct(1, block_lengths.data(), displacements.data(), typelist.data(), &Top_Ghost_Type);
	// Commit to use
	MPI_Type_commit(&Top_Ghost_Type);



	// ------------------------------------------------------------------------------------------------------------------------|
	// MPI_DATATYPE 7 : BOTTOM - REGULAR
	// ------------------------------------------------------------------------------------------------------------------------|
	// Initialize variables
	typelist[0] = MPI_DOUBLE;
	block_lengths[0] = Dom.sub_jmax - 2;
	MPI_Get_address(&Dom.current_grid_2d[Dom.sub_imax - 2][1], &displacements[0]);
	displacements[0] = displacements[0] - add_start;
	// Create the dataType
	MPI_Type_create_struct(1, block_lengths.data(), displacements.data(), typelist.data(), &Bottom_Type);
	// Commit to use
	MPI_Type_commit(&Bottom_Type);




	// ------------------------------------------------------------------------------------------------------------------------|
	// MPI_DATATYPE 8 : BOTTOM - GHOST
	// ------------------------------------------------------------------------------------------------------------------------|
	// Initialize variables
	typelist[0] = MPI_DOUBLE;
	block_lengths[0] = Dom.sub_jmax - 2;
	MPI_Get_address(&Dom.current_grid_2d[Dom.sub_imax - 1][1], &displacements[0]);
	displacements[0] = displacements[0] - add_start;
	// Create the dataType
	MPI_Type_create_struct(1, block_lengths.data(), displacements.data(), typelist.data(), &Bottom_Ghost_Type);
	// Commit to use
	MPI_Type_commit(&Bottom_Ghost_Type);

}
