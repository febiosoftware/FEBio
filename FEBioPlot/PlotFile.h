#pragma once
#include "FECore/FEMesh.h"
#include "FECore/FEPlotData.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class implements the facilities to write to a plot database. 
//!
class PlotFile
{
public:
	//! constructor
	PlotFile();

	//! descructor
	virtual ~PlotFile();

	//! close the plot database
	virtual void Close();

	//! Open the plot database
	virtual bool Open(FEModel& fem, const char* szfile) = 0;

	//! Open for appending
	virtual bool Append(FEModel& fem, const char* szfile) = 0;

	//! Write current FE state to plot database
	virtual bool Write(FEModel& fem, float ftime) = 0;

	//! see if the plot file is valid
	virtual bool IsValid() const = 0;

protected:
	FEModel*	m_pfem;		//!< pointer to FE model
};
