// PlotFile.h: interface for the PlotFile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PLOTFILE_H__6E7170ED_6C03_4720_96CF_C53411A7464E__INCLUDED_)
#define AFX_PLOTFILE_H__6E7170ED_6C03_4720_96CF_C53411A7464E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

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
	virtual bool Write(FEModel& fem) = 0;

	//! see if the plot file is valid
	virtual bool IsValid() const = 0;

protected:
	FEModel*	m_pfem;		//!< pointer to FE model
};

#endif // !defined(AFX_PLOTFILE_H__6E7170ED_6C03_4720_96CF_C53411A7464E__INCLUDED_)
