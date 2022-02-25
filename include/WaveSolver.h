#pragma once
#define _USE_MATH_DEFINES
#include "MPIDatatype.h"
#include "Domain.h"



//--------------------------------------------------------------------------------|
//  FUNCTIONS DEALING WITH BOUNDRY CONDITIONS
// -------------------------------------------------------------------------------|
void selectBoundary(Domain<double>& Dom, int boundary);

void dirichlitBoundary(Domain<double>& Dom);

void neumannBouandry(Domain<double>& Dom);



//--------------------------------------------------------------------------------|
//  FUNCTIONS DEALING WITH COMMUNICATIONS
// -------------------------------------------------------------------------------|

// Dirichlit and Neumann Communications-------------------------------------------|
void doMPICommunications(Domain<double>& domain_instance, 
	MPI_Data<double>& MData, int count, MPI_Request* request_list);


//--------------------------------------------------------------------------------|
// FUNCTION TO SET INITIAL CONDITION
// -------------------------------------------------------------------------------|
void setInitialDisturbance(Domain<double> &Dom);
// Array pointer to store thenumber of splashes


//--------------------------------------------------------------------------------|
// FUNCTIONS DEALTING WITH TIME STEPPING
// -------------------------------------------------------------------------------|
// Solve inner domain
void doInnerIterations(Domain<double>&Dom);

// Solve outer domain
void doBoundaryIterations(Domain<double> &Dom);


//--------------------------------------------------------------------------------|
// TEST AND POSTPROCESSING FUNCTIONS 
// -------------------------------------------------------------------------------|
void fillSubDomainsWithProcessId(Domain<double>& Dom);

void grid_to_file(Domain<double>& Dom, int out, int boundary);