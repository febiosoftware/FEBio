#include "stdafx.h"
#include "FETimeInfo.h"
#include "DumpStream.h"

FETimeInfo::FETimeInfo()
{
	currentTime = 0.0;
	timeIncrement = 0.0;
	alpha = 1.0;
	beta = 0.25;
	gamma = 0.5;
	alphaf = 1.0;
	alpham = 1.0;
	currentIteration = 0;
}

FETimeInfo::FETimeInfo(double time, double tinc)
{
	currentTime = time;
	timeIncrement = tinc;
	alpha = 1.0;
	beta = 0.25;
	gamma = 0.5;
	alphaf = 1.0;
	alpham = 1.0;
	currentIteration = 0;
}

void FETimeInfo::Serialize(DumpStream& ar)
{
	ar & currentTime;
	ar & timeIncrement;
	ar & currentIteration;

	ar & alpha;
	ar & alphaf;
	ar & alpham;
	ar & beta;
	ar & gamma;
}
