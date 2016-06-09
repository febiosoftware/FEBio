#pragma once

//-----------------------------------------------------------------------------
struct FETimeInfo
{
	FETimeInfo()
	{
		currentTime = 0.0;
		timeIncrement = 0.0;
		alpha = 1.0;
		beta = 0.25;
		gamma = 0.5;
		currentIteration = 0;
	}

	FETimeInfo(double time, double tinc)
	{
		currentTime = time;
		timeIncrement = tinc;
		alpha = 1.0;
		beta = 0.25;
		gamma = 0.5;
		currentIteration = 0;
	}

	double	currentTime;	//!< current time value
	double	timeIncrement;		//!< current time step (difference between this time and previous one)

	// HHT time integration parameters
	double	alpha;
	double	beta;
	double	gamma;

	int		currentIteration;	//!< iteration number
};
