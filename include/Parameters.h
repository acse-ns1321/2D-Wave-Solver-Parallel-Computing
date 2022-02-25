#pragma once


// Processor Parameters
extern int id, p;
extern int tag_num;

// Domain Parameters - grid and spacing
extern int imax , jmax ;
extern double y_max , x_max , dx, dy;

// Time Stepping Parameters
extern double t_max;
extern double t, t_out, dt_out , dt;

// Wave Equation parameters
extern double c;

// Boundary Conditions
extern int boundary;

// Initial Conditions
extern int boundary;
