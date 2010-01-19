#pragma once
#include "PlotFile.h"

//-----------------------------------------------------------------------------
//! This class stores the results of the analysis in the LSDYNA database format.
//!
class LSDYNAPlotFile :	public PlotFile
{
protected:

	//-----------------------------------------------------------------------------
	//! This is the header of the plot database. 
	
	//! The header consists of 64 integers, however
	//! the first 10 integers are to be interpreted as characters (single bytes),
	//! since they contain the title of the problem. Note that the title is not
	//! zero-terminated.

	struct PLOTHEADER
	{
		char	Title[40];		//!< title of the problem
		int		UnUsed0[5];		//!< blanks (unused)
		int		ndim;			//!< number of dimensions
		int		nump;			//!< number of nodal points
		int		icode;			//!< code descriptor
		int		nglbv;			//!< number of global state variables
		int		flagT;			//!< state nodal temperatures included ?
		int		flagU;			//!< state nodal coordinates included ?
		int		flagV;			//!< state nodal velocities included ?
		int		flagA;			//!< state nodal accelerations included ?
		int		nel8;			//!< number of 8-node hexahedral elements
		int		nummat8;		//!< number of materials used by hexahedral elements
		int		UnUsed1[2];		//!< blanks (unused)
		int		nv3d;			//!< number of variables for hexahedral elements
		int		nel2;			//!< number of 2-node beam elements
		int		nummat2;		//!< number of materials used by beam elements
		int		nv1d;			//!< number of variables for beam elements
		int		nel4;			//!< number of 4-node shell elements
		int		nummat4;		//!< number of materials used by shell elements
		int		nv2d;			//!< number of variables for shell elements
		int		neiph;			//!< number of additional variables per solid element
		int		neips;			//!< number of additional variables per shell integration point
		int		maxint;			//!< number of integration points dumped for each shell
		int		UnUsed3[7];		//!< blank (unused)
		int		ioshl1;			//!< 6 stress component flag for shells
		int		ioshl2;			//!< plastic strain flag for shells
		int		ioshl3;			//!< shell force resultant flag
		int		ioshl4;			//!< shell thickness, energy + 2 more
		int		UnUsed4[16];	//!< blank (unused)
	};

public:
	LSDYNAPlotFile(void);
	~LSDYNAPlotFile(void);

	//! Open the plot database
	bool Open(FEM& fem, const char* szfile);

	//! Open for appending
	bool Append(FEM& fem, const char* szfile);

	//! Write current FE state to plot database
	bool Write(FEM& fem);

protected:
	// vector fields
	void write_displacements();
	void write_velocities();
	void write_accelerations();
	void write_fluid_flux();
	void write_contact_tractions();
	void write_reaction_forces();
	void write_material_fibers();
	void write_heat_flux();

	// scalar fields
	void write_fluid_pressures();
	void write_contact_pressures();
	void write_contact_gaps();
	void write_temperatures();

	// plastic stress fields
	float fiber_strain(FESolidElement& el, int j);
	float dev_fiber_strain(FESolidElement& el, int j);

public:
	bool	m_bsstrn;		//!< shell strain flag
	int		m_nfield[5];	//!< field maps

protected:
	PLOTHEADER	m_ph;	//!< The plot file header
};
