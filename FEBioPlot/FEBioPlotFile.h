#pragma once
#include "PlotFile.h"
#include "PltArchive.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEShellDomain.h"
#include "FECore/FETrussDomain.h"
#include "FECore/FEDiscreteDomain.h"
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
// VERSION INFO
//==============
// -0.4: fixed bug with the PLT_MAX_FACET_NODES: Previous version was writing the 
//       surface facets using a buffer of with PLT_MAX_FACET_NODES+2, but was only
//       writing PLT_MAX_FACET_NODES entries, essentially cutting of two nodes for 
//       eight node facets.
// -0.3: added PLT_HDR_COMPRESSION flag and compression for state data
// -0.2: increased max surface nodes to 8 to accomodate tri6 and quad8 facets.
//       facets store nr of nodes for each face PLT_FACE
// -0.1: initial design

//-----------------------------------------------------------------------------
//! This class implements the facilities to export FE data in the FEBio
//! plot file format.
//!
class FEBioPlotFile : public PlotFile
{
public:
	// file version
	enum { PLT_VERSION = 0x0004 };

	// max nodes per facet
	enum { PLT_MAX_FACET_NODES = 9 };

	// file tags
	enum { 
		PLT_ROOT						= 0x01000000,
		PLT_HEADER						= 0x01010000,
			PLT_HDR_VERSION				= 0x01010001,
			PLT_HDR_NODES				= 0x01010002,
			PLT_HDR_MAX_FACET_NODES		= 0x01010003,
			PLT_HDR_COMPRESSION			= 0x01010004,	// added in version 0x03
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
					PLT_DOM_NAME		= 0x01032105,	// added in version 0x03
				PLT_DOM_ELEM_LIST		= 0x01042200,
					PLT_ELEMENT			= 0x01042201,
			PLT_SURFACE_SECTION			= 0x01043000,
				PLT_SURFACE				= 0x01043100,
				PLT_SURFACE_HDR			= 0x01043101,
					PLT_SURFACE_ID		= 0x01043102,
					PLT_SURFACE_FACES	= 0x01043103,
					PLT_SURFACE_NAME	= 0x01043104,	// added in version 0x03
				PLT_FACE_LIST			= 0x01043200,
					PLT_FACE			= 0x01043201, 
			PLT_NODESET_SECTION			= 0x01044000,	// added in version 0x03
				PLT_NODESET				= 0x01044100,	// added in version 0x03
				PLT_NODESET_HDR			= 0x01044101,	// added in version 0x03
					PLT_NODESET_ID		= 0x01044102,	// added in version 0x03
					PLT_NODESET_NAME	= 0x01044103,	// added in version 0x03
					PLT_NODESET_SIZE	= 0x01044104,	// added in version 0x03
				PLT_NODESET_LIST		= 0x01044200,	// added in version 0x03
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
	enum Elem_Type { 
		PLT_ELEM_HEX, 
		PLT_ELEM_PENTA, 
		PLT_ELEM_TET, 
		PLT_ELEM_QUAD, 
		PLT_ELEM_TRI, 
		PLT_ELEM_TRUSS, 
		PLT_ELEM_HEX20, 
		PLT_ELEM_TET10, 
		PLT_ELEM_TET15, 
		PLT_ELEM_HEX27 
	};

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
		bool AddVariable(FEModel* pfem, const char* szname, vector<int>& item);

		int NodalVariables  () { return m_Node.size(); }
		int DomainVarialbes () { return m_Elem.size(); }
		int SurfaceVariables() { return m_Face.size(); }

		void Defaults(FEModel& fem);

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
		bool AddNodalVariable   (FEPlotData* ps, const char* szname, vector<int>& item);
		bool AddDomainVariable  (FEPlotData* ps, const char* szname, vector<int>& item);
		bool AddSurfaceVariable (FEPlotData* ps, const char* szname, vector<int>& item);

	protected:
		list<DICTIONARY_ITEM>	m_Glob;		// Global variables
		list<DICTIONARY_ITEM>	m_Mat;		// Material variables
		list<DICTIONARY_ITEM>	m_Node;		// Node variables
		list<DICTIONARY_ITEM>	m_Elem;		// Domain variables
		list<DICTIONARY_ITEM>	m_Face;		// Surface variables

		friend class FEBioPlotFile;
	};

public:
	FEBioPlotFile(FEModel& fem);
	~FEBioPlotFile(void);

	//! Open the plot database
	bool Open(FEModel& fem, const char* szfile);

	//! Close the plot database
	void Close();

	//! Open for appending
	bool Append(FEModel& fem, const char* szfile);

	//! Write current FE state to plot database
	bool Write(FEModel& fem);

	//! Add a variable to the dictionary
	bool AddVariable(const char* sz, vector<int>& item) { return m_dic.AddVariable(&m_fem, sz, item); }

	//! Set the compression level
	void SetCompression(int n);

	//! see if the plot file is valid
	virtual bool IsValid() const;

public:
	const Dictionary& GetDictionary() const { return m_dic; }

protected:
	bool WriteRoot      (FEModel& fem);
	bool WriteHeader    (FEModel& fem);
	bool WriteDictionary(FEModel& fem);
	bool WriteMaterials (FEModel& fem);
	bool WriteGeometry  (FEModel& fem);

	void WriteDicList(list<DICTIONARY_ITEM>& dic);

	void WriteNodeSection   (FEMesh& m);
	void WriteDomainSection (FEMesh& m);
	void WriteSurfaceSection(FEMesh& m);
	void WriteNodeSetSection(FEMesh& m);

	void WriteSolidDomain   (FESolidDomain&    dom);
	void WriteShellDomain   (FEShellDomain&    dom);
	void WriteTrussDomain   (FETrussDomain&    dom);
	void WriteDiscreteDomain(FEDiscreteDomain& dom);

	void WriteGlobalData  (FEModel& fem);
	void WriteMaterialData(FEModel& fem);
	void WriteNodeData    (FEModel& fem);
	void WriteDomainData  (FEModel& fem);
	void WriteSurfaceData (FEModel& fem);

protected:
	bool ReadDictionary();
	bool ReadDicList();

protected:
	Dictionary	m_dic;	// dictionary
	PltArchive	m_ar;	// the data archive
	FEModel&	m_fem;
	int			m_ncompress;	// compression level
};
