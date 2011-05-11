// PlotFile.h: interface for the PlotFile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PLOTFILE_H__6E7170ED_6C03_4720_96CF_C53411A7464E__INCLUDED_)
#define AFX_PLOTFILE_H__6E7170ED_6C03_4720_96CF_C53411A7464E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEMesh.h"
#include "FEPlotData.h"

class FEM;

// empty field
#define PLOT_NONE		0

// scalar fields
#define PLOT_FLUID_PRESSURE		1
#define PLOT_CONTACT_PRESSURE	2
#define PLOT_CONTACT_GAP		3
#define	PLOT_PLASTIC_STRAIN		4
#define PLOT_FIBER_STRAIN		5
#define PLOT_DEV_FIBER_STRAIN	6
#define PLOT_TEMPERATURE		7


// vector fields
#define PLOT_DISPLACEMENT		1
#define PLOT_VELOCITY			2
#define PLOT_ACCELERATION		3
#define	PLOT_FLUID_FLUX			4
#define PLOT_CONTACT_TRACTION	5
#define	PLOT_REACTION_FORCE		6
#define	PLOT_MATERIAL_FIBER		7
#define	PLOT_HEAT_FLUX			8

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
	void Close();

	//! Open the plot database
	virtual bool Open(FEM& fem, const char* szfile) = 0;

	//! Open for appending
	virtual bool Append(FEM& fem, const char* szfile) = 0;

	//! Write current FE state to plot database
	virtual bool Write(FEM& fem) = 0;

protected:
	FEM*	m_pfem;		//!< pointer to FE model
};

#endif // !defined(AFX_PLOTFILE_H__6E7170ED_6C03_4720_96CF_C53411A7464E__INCLUDED_)
