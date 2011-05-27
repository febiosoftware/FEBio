#pragma once
#include "PlotFile.h"
#include "FECore/Archive.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEShellDomain.h"
#include "FECore/FETrussDomain.h"
#include "FECore/FEDiscreteDomain.h"
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
		PLT_ROOT						= 0x01000000,
		PLT_HEADER						= 0x01010000,
			PLT_HDR_VERSION				= 0x01010001,
			PLT_HDR_NODES				= 0x01010002,
		PLT_DICTIONARY					= 0x01020000,
			PLT_DIC_ITEM				= 0x01020001,
			PLT_DIC_ITEM_TYPE			= 0x01020002,
			PLT_DIC_ITEM_FMT			= 0x01020003,
			PLT_DIC_ITEM_NAME			= 0x01020004,
			PLT_DIC_GLOBAL				= 0x01021000,
			PLT_DIC_MATERIAL			= 0x01022000,
			PLT_DIC_NODAL				= 0x01023000,
			PLT_DIC_DOMAIN				= 0x01024000,
			PLT_DIC_SURFACE				= 0x01025000,
		PLT_MATERIALS					= 0x01030000,
			PLT_MATERIAL				= 0x01030001,
			PLT_MAT_ID					= 0x01030002,
			PLT_MAT_NAME				= 0x01030003,
		PLT_GEOMETRY					= 0x01040000,
			PLT_NODE_SECTION			= 0x01041000,
				PLT_NODE_COORDS			= 0x01041001,
			PLT_DOMAIN_SECTION			= 0x01042000,
				PLT_DOMAIN				= 0x01042100,
				PLT_DOMAIN_HDR			= 0x01042101,
					PLT_DOM_ELEM_TYPE	= 0x01042102,
					PLT_DOM_MAT_ID		= 0x01042103,
					PLT_DOM_ELEMS		= 0x01032104,
				PLT_DOM_ELEM_LIST		= 0x01042200,
					PLT_ELEMENT			= 0x01042201,
			PLT_SURFACE_SECTION			= 0x01043000,
				PLT_SURFACE				= 0x01043100,
				PLT_SURFACE_HDR			= 0x01043101,
					PLT_SURFACE_ID		= 0x01043102,
					PLT_SURFACE_FACES	= 0x01043103,
				PLT_FACE_LIST			= 0x01043200,
					PLT_FACE			= 0x01043201,
		PLT_STATE						= 0x02000000,
			PLT_STATE_HEADER			= 0x02010000,
				PLT_STATE_HDR_ID		= 0x02010001,
				PLT_STATE_HDR_TIME		= 0x02010002,
			PLT_STATE_DATA				= 0x02020000,
				PLT_STATE_VARIABLE		= 0x02020001,
				PLT_STATE_VAR_ID		= 0x02020002,
				PLT_STATE_VAR_DATA		= 0x02020003,
				PLT_GLOBAL_DATA			= 0x02020100,
				PLT_MATERIAL_DATA		= 0x02020200,
				PLT_NODE_DATA			= 0x02020300,
				PLT_ELEMENT_DATA		= 0x02020400,
				PLT_FACE_DATA			= 0x02020500
	};

	// --- element types ---
	enum Elem_Type { PLT_ELEM_HEX, PLT_ELEM_PENTA, PLT_ELEM_TET, PLT_ELEM_QUAD, PLT_ELEM_TRI, PLT_ELEM_TRUSS };

	// size of name variables
	enum { STR_SIZE = 64 };

public:
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
		bool AddVariable(const char* szname);

		int NodalVariables  () { return m_Node.size(); }
		int DomainVarialbes () { return m_Elem.size(); }
		int SurfaceVariables() { return m_Face.size(); }

		void Defaults(FEM& fem);

		void Clear();

	public:
		const list<DICTIONARY_ITEM>& GlobalVariableList  () const { return m_Glob; }
		const list<DICTIONARY_ITEM>& MaterialVariableList() const { return m_Mat;  }
		const list<DICTIONARY_ITEM>& NodalVariableList   () const { return m_Node; }
		const list<DICTIONARY_ITEM>& DomainVariableList  () const { return m_Elem; }
		const list<DICTIONARY_ITEM>& SurfaceVariableList () const { return m_Face; }

	protected:
		bool AddGlobalVariable  (FEPlotData* ps, const char* szname);
		bool AddMaterialVariable(FEPlotData* ps, const char* szname);
		bool AddNodalVariable   (FEPlotData* ps, const char* szname);
		bool AddDomainVariable  (FEPlotData* ps, const char* szname);
		bool AddSurfaceVariable (FEPlotData* ps, const char* szname);

	protected:
		list<DICTIONARY_ITEM>	m_Glob;		// Global variables
		list<DICTIONARY_ITEM>	m_Mat;		// Material variables
		list<DICTIONARY_ITEM>	m_Node;		// Node variables
		list<DICTIONARY_ITEM>	m_Elem;		// Domain variables
		list<DICTIONARY_ITEM>	m_Face;		// Surface variables

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

	//! Add a variable to the dictionary
	bool AddVariable(const char* sz) { return m_dic.AddVariable(sz); }

public:
	const Dictionary& GetDictionary() const { return m_dic; }

protected:
	bool WriteHeader    (FEM& fem);
	bool WriteDictionary(FEM& fem);
	bool WriteMaterials (FEM& fem);
	bool WriteGeometry  (FEM& fem);

	void WriteDicList(list<DICTIONARY_ITEM>& dic);

	void WriteNodeSection   (FEMesh& m);
	void WriteDomainSection (FEMesh& m);
	void WriteSurfaceSection(FEMesh& m);

	void WriteSolidDomain   (FESolidDomain&    dom);
	void WriteShellDomain   (FEShellDomain&    dom);
	void WriteTrussDomain   (FETrussDomain&    dom);
	void WriteDiscreteDomain(FEDiscreteDomain& dom);

	void WriteGlobalData  (FEM& fem);
	void WriteMaterialData(FEM& fem);
	void WriteNodeData    (FEM& fem);
	void WriteDomainData  (FEM& fem);
	void WriteSurfaceData (FEM& fem);

protected:
	bool ReadDictionary();
	bool ReadDicList();

protected:
	Dictionary	m_dic;	// dictionary
	Archive		m_ar;	// the data archive
};
