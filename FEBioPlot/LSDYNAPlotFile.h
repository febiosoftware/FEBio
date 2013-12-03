#pragma once
#include "FEBioPlot/PlotFile.h"

//-----------------------------------------------------------------------------
// Field tags

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
	bool Open(FEModel& fem, const char* szfile);

	//! close the plot database
	void Close();

	//! Open for appending
	bool Append(FEModel& fem, const char* szfile);

	//! Write current FE state to plot database
	bool Write(FEModel& fem);

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

	// tensor fields
	void write_solid_stress();
	void write_shell_stress();
	void write_truss_stress();

public:
	bool	m_bsstrn;		//!< shell strain flag
	int		m_nfield[5];	//!< field maps

protected:
	PLOTHEADER	m_ph;	//!< The plot file header
	FILE*		m_fp;	//!< the file pointer
};
