#pragma once
#include "PlotFile.h"
#include <list>
using namespace std;

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
		unsigned int	m_ntype;
		char			m_szname[DI_NAME_SIZE];
	};

	class Dictionary
	{
	public:
		void AddGlobalVariable(unsigned int ntype, const char* szname);
		void AddNodalVariable (unsigned int ntype, const char* szname);
		void AddSolidVariable (unsigned int ntype, const char* szname);
		void AddShellVariable (unsigned int ntype, const char* szname);
		void AddBeamVariable  (unsigned int ntype, const char* szname);

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
	void write_stresses();

protected:
	HEADER		m_hdr;	// plot file header
	Dictionary	m_dic;	// dictionary
};
