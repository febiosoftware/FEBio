#pragma once
#include "PlotFile.h"
#include "Archive.h"
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
//! This class implements the facilities to export FE data in the FEBio
//! plot file format.
//!
class FEBioPlotFile : public PlotFile
{
protected:
	// file version
	enum { PLT_VERSION = 0x0001 };

	// file tags
	enum { 
		PLT_ROOT					= 0x000AC996,		// = 706966 ('F', 'E', 'B')
		PLT_HEADER					= 0x00100000,
			PLT_HDR_VERSION			= 0x00100001,
		PLT_DICTIONARY				= 0x00200000,
			PLT_DIC_ITEM			= 0x00200001,
			PLT_DIC_ITEM_TYPE		= 0x00200002,
			PLT_DIC_ITEM_FMT		= 0x00200003,
			PLT_DIC_ITEM_NAME		= 0x00200004,
			PLT_DIC_GLOBAL			= 0x00201000,
			PLT_DIC_MATERIAL		= 0x00202000,
			PLT_DIC_NODAL			= 0x00203000,
			PLT_DIC_DOMAIN			= 0x00204000,
			PLT_DIC_SURFACE			= 0x00205000,
		PLT_MATERIALS				= 0x00300000,
			PLT_MATERIAL			= 0x00300001,
			PLT_MAT_ID				= 0x00300002,
			PLT_MAT_NAME			= 0x00300003,
		PLT_GEOMETRY				= 0x00400000,
			PLT_NODE_SECTION		= 0x00401000,
				PLT_NODE_COORDS		= 0x00401001,
			PLT_DOMAIN_SECTION		= 0x00402000,
				PLT_DOMAIN			= 0x00402001,
				PLT_DOM_ELEM_TYPE	= 0x00402002,
				PLT_DOM_MAT_ID		= 0x00402003,
				PLT_DOM_ELEM_LIST	= 0x00402004,
			PLT_SURFACE_SECTION		= 0x00403000,
		PLT_STATE					= 0x00500000,
			PLT_STATE_HEADER		= 0x00501000,
				PLT_STATE_HDR_ID	= 0x00501001,
				PLT_STATE_HDR_TIME	= 0x00501002,
			PLT_STATE_DATA			= 0x00502000,
				PLT_STATE_VARIABLE	= 0x00502001,
				PLT_STATE_VAR_ID	= 0x00502002,
				PLT_STATE_VAR_DATA	= 0x00502003,
				PLT_GLOBAL_DATA		= 0x00502100,
				PLT_MATERIAL_DATA	= 0x00502200,
				PLT_NODE_DATA		= 0x00502300,
				PLT_ELEMENT_DATA	= 0x00502400,
				PLT_FACE_DATA		= 0x00502500
	};

	// --- element types ---
	enum Elem_Type { PLT_ELEM_HEX, PLT_ELEM_PENTA, PLT_ELEM_TET, PLT_ELEM_QUAD, PLT_ELEM_TRI, PLT_ELEM_TRUSS };

	// size of name variables
	enum { STR_SIZE = 64 };

	// Dictionary entry
	struct DICTIONARY_ITEM
	{
		FEPlotData*		m_psave;
		unsigned int	m_ntype;	// data type
		unsigned int	m_nfmt;		// storage format
		char			m_szname[STR_SIZE];
	};

	class Dictionary
	{
	public:
		void AddGlobalVariable  (FEPlotData* ps, const char* szname);
		void AddMaterialVariable(FEPlotData* ps, const char* szname);
		void AddNodalVariable   (FEPlotData* ps, const char* szname);
		void AddElementVariable (FEPlotData* ps, const char* szname);
		void AddFaceVariable    (FEPlotData* ps, const char* szname);

	protected:
		list<DICTIONARY_ITEM>	m_Glob;		// Global variables
		list<DICTIONARY_ITEM>	m_Mat;		// Material variables
		list<DICTIONARY_ITEM>	m_Node;		// Node variables
		list<DICTIONARY_ITEM>	m_Elem;		// Element variables
		list<DICTIONARY_ITEM>	m_Face;		// Face variables

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
	bool WriteHeader    (FEM& fem);
	bool WriteDictionary(FEM& fem);
	bool WriteMaterials (FEM& fem);
	bool WriteGeometry  (FEM& fem);

	void WriteDicList(list<DICTIONARY_ITEM>& dic);

	void WriteNodeSection   (FEMesh& m);
	void WriteDomainSection (FEMesh& m);
	void WriteSurfaceSection(FEMesh& m);

	void WriteSolidDomain(FESolidDomain& dom);
	void WriteShellDomain(FEShellDomain& dom);
	void WriteTrussDomain(FETrussDomain& dom);

	void WriteGlobalData  (FEM& fem);
	void WriteMaterialData(FEM& fem);
	void WriteNodeData    (FEM& fem);
	void WriteElementData (FEM& fem);
	void WriteFaceData    (FEM& fem);

protected:
	Dictionary	m_dic;	// dictionary
	Archive		m_ar;	// the data archive
};
