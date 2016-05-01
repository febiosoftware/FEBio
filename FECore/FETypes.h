#pragma once

//-----------------------------------------------------------------------------
struct FETimePoint
{
	FETimePoint()
	{
		t = dt = 0.0;
		alpha = 1.0;
		beta = 0.25;
		gamma = 0.5;
	}

	FETimePoint(double time, double tinc)
	{
		t = time;
		dt = tinc;
		alpha = 1.0;
		beta = 0.25;
		gamma = 0.5;
	}


	double	t;		// current time value
	double	dt;		// current time step (difference between this time and previous one)

	// HHT time integration parameters
	double	alpha;
	double	beta;
	double	gamma;
};
