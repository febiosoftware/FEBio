#pragma once
#include "PlotFile.h"
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
//! This class implements the facilities to export FE data in the FEBio
//! plot file format.
//!
class FEBioPlotFile : public PlotFile
{
protected:
	// FEBio ID tag
	enum { FEBIO_TAG = 706966 };

	// variable types
	enum Var_Type { SCALAR, VEC3F, MAT3FS };

	// HEADER structure
	struct HEADER
	{
		int		nsize;			// sizeof(HEADER)
		int		nnodes;			// number of nodes
		int		n3d;			// number of solid elements
		int		n2d;			// number of shell elements
		int		n1d;			// number of beam elements
		int		nglv;			// number of global variables
		int		nnv;			// number of nodal variables
		int		nv3d;			// number of variables for solid elements
		int		nv2d;			// number of variables for shell elements
		int		nv1d;			// number of variables for beam elements
		int		nmat;			// number of parts
		
		int		nreserved[53];	// reverved for future use
	};

	// size of name variables
	enum { DI_NAME_SIZE = 64 };

	// Dictionary entry
	struct DICTIONARY_ITEM
	{
		FEPlotData*		m_psave;
		unsigned int	m_ntype;
		char			m_szname[DI_NAME_SIZE];
	};

	class Dictionary
	{
	public:
		void AddGlobalVariable(FEPlotData* ps, unsigned int ntype, const char* szname);
		void AddNodalVariable (FEPlotData* ps, unsigned int ntype, const char* szname);
		void AddSolidVariable (FEPlotData* ps, unsigned int ntype, const char* szname);
		void AddShellVariable (FEPlotData* ps, unsigned int ntype, const char* szname);
		void AddBeamVariable  (FEPlotData* ps, unsigned int ntype, const char* szname);

	protected:
		void Save(Archive& ar);
		
	protected:
		list<DICTIONARY_ITEM>	m_Glob;
		list<DICTIONARY_ITEM>	m_Node;
		list<DICTIONARY_ITEM>	m_Elem;
		list<DICTIONARY_ITEM>	m_Shell;
		list<DICTIONARY_ITEM>	m_Beam;

		friend class FEBioPlotFile;
	};

public:
	FEBioPlotFile(void);
	~FEBioPlotFile(void);

	//! Open the plot database
	bool Open(FEM& fem, const char* szfile);

	//! Open for appending
	bool Append(FEM& fem, const char* szfile);

	//! Write current FE state to plot database
	bool Write(FEM& fem);

protected:
	HEADER		m_hdr;	// plot file header
	Dictionary	m_dic;	// dictionary
};
