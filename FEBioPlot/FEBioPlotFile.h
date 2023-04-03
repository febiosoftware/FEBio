/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "PlotFile.h"
#include "PltArchive.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEShellDomain.h"
#include "FECore/FEBeamDomain.h"
#include "FECore/FEDiscreteDomain.h"
#include "FECore/FEDomain2D.h"
#include <list>

//-----------------------------------------------------------------------------
//! This class implements the facilities to export FE data in the FEBio
//! plot file format (version 3).
//!
class FEBioPlotFile : public PlotFile
{
public:
	// file version
	// 32: added PLT_ELEMENTSET_SECTION
	enum { PLT_VERSION = 0x0032 };

	// file tags
	enum { 
		PLT_ROOT						= 0x01000000,
		PLT_HEADER						= 0x01010000,
			PLT_HDR_VERSION				= 0x01010001,
//			PLT_HDR_NODES				= 0x01010002,
//			PLT_HDR_MAX_FACET_NODES		= 0x01010003,	// removed (redefined in seach SURFACE section)
			PLT_HDR_COMPRESSION			= 0x01010004,
			PLT_HDR_AUTHOR				= 0x01010005,	// new in 2.0
			PLT_HDR_SOFTWARE			= 0x01010006,	// new in 2.0
			PLT_HDR_UNITS				= 0x01010007,	// new in 4.0
		PLT_DICTIONARY					= 0x01020000,
			PLT_DIC_ITEM				= 0x01020001,
			PLT_DIC_ITEM_TYPE			= 0x01020002,
			PLT_DIC_ITEM_FMT			= 0x01020003,
			PLT_DIC_ITEM_NAME			= 0x01020004,
			PLT_DIC_ITEM_ARRAYSIZE		= 0x01020005,	// added in version 0x05
			PLT_DIC_ITEM_ARRAYNAME		= 0x01020006,	// added in version 0x05
			PLT_DIC_ITEM_UNITS			= 0x01020007,	// added in version 4.0
			PLT_DIC_GLOBAL				= 0x01021000,
//			PLT_DIC_MATERIAL			= 0x01022000,	// this was removed
			PLT_DIC_NODAL				= 0x01023000,
			PLT_DIC_DOMAIN				= 0x01024000,
			PLT_DIC_SURFACE				= 0x01025000,
//		PLT_MATERIALS					= 0x01030000,		// This was removed
//			PLT_MATERIAL				= 0x01030001,
//			PLT_MAT_ID					= 0x01030002,
//			PLT_MAT_NAME				= 0x01030003,
		PLT_MESH						= 0x01040000,		// this was PLT_GEOMETRY
			PLT_NODE_SECTION			= 0x01041000,
				PLT_NODE_HEADER			= 0x01041100,		// new in 2.0
					PLT_NODE_SIZE		= 0x01041101,		// new in 2.0
					PLT_NODE_DIM		= 0x01041102,		// new in 2.0
					PLT_NODE_NAME		= 0x01041103,		// new in 2.0
				PLT_NODE_COORDS			= 0x01041200,		// new in 2.0
			PLT_DOMAIN_SECTION			= 0x01042000,
				PLT_DOMAIN				= 0x01042100,
				PLT_DOMAIN_HDR			= 0x01042101,
					PLT_DOM_ELEM_TYPE	= 0x01042102,
					PLT_DOM_PART_ID		= 0x01042103,		// this was PLT_DOM_MAT_ID
					PLT_DOM_ELEMS		= 0x01032104,
					PLT_DOM_NAME		= 0x01032105,
				PLT_DOM_ELEM_LIST		= 0x01042200,
					PLT_ELEMENT			= 0x01042201,
			PLT_SURFACE_SECTION			= 0x01043000,
				PLT_SURFACE				= 0x01043100,
				PLT_SURFACE_HDR			= 0x01043101,
					PLT_SURFACE_ID		= 0x01043102,
					PLT_SURFACE_FACES	= 0x01043103,
					PLT_SURFACE_NAME	= 0x01043104,
					PLT_SURFACE_MAX_FACET_NODES = 0x01043105,	// new in 2.0 (max number of nodes per facet)
				PLT_FACE_LIST			= 0x01043200,
					PLT_FACE			= 0x01043201,
			PLT_NODESET_SECTION			= 0x01044000,
				PLT_NODESET				= 0x01044100,
				PLT_NODESET_HDR			= 0x01044101,
					PLT_NODESET_ID		= 0x01044102,
					PLT_NODESET_NAME	= 0x01044103,
					PLT_NODESET_SIZE	= 0x01044104,
				PLT_NODESET_LIST		= 0x01044200,
			PLT_PARTS_SECTION			= 0x01045000,		// new in 2.0
				PLT_PART				= 0x01045100,
				PLT_PART_ID				= 0x01045101,
				PLT_PART_NAME			= 0x01045102,

			// element set section was added in 4.1
			PLT_ELEMENTSET_SECTION		= 0x01046000,
				PLT_ELEMENTSET			= 0x01046100,
				PLT_ELEMENTSET_HDR		= 0x01046101,
					PLT_ELEMENTSET_ID	= 0x01046102,
					PLT_ELEMENTSET_NAME	= 0x01046103,
					PLT_ELEMENTSET_SIZE	= 0x01046104,
				PLT_ELEMENTSET_LIST		= 0x01046200,

			// facet set section was added in 4.1
			PLT_FACETSET_SECTION			= 0x01047000,
				PLT_FACETSET				= 0x01047100,
				PLT_FACETSET_HDR			= 0x01047101,
					PLT_FACETSET_ID			= 0x01047102,
					PLT_FACETSET_NAME		= 0x01047103,
					PLT_FACETSET_SIZE		= 0x01047104,
					PLT_FACETSET_MAXNODES	= 0x01047105,
				PLT_FACETSET_LIST			= 0x01047200,
					PLT_FACET				= 0x01047201,

			// plot objects were added in 3.0
			PLT_OBJECTS_SECTION			= 0x01050000,
					PLT_OBJECT_ID		= 0x01050001,
					PLT_OBJECT_NAME		= 0x01050002,
					PLT_OBJECT_TAG		= 0x01050003,
					PLT_OBJECT_POS		= 0x01050004,
					PLT_OBJECT_ROT		= 0x01050005,
					PLT_OBJECT_DATA		= 0x01050006,
				PLT_POINT_OBJECT		= 0x01051000,
					PLT_POINT_COORD		= 0x01051001,
				PLT_LINE_OBJECT			= 0x01052000,
					PLT_LINE_COORDS		= 0x01052001,

		PLT_STATE						= 0x02000000,
			PLT_STATE_HEADER			= 0x02010000,
				PLT_STATE_HDR_ID		= 0x02010001,
				PLT_STATE_HDR_TIME		= 0x02010002,
				PLT_STATE_STATUS        = 0x02010003,	// new in 3.1
			PLT_STATE_DATA				= 0x02020000,
				PLT_STATE_VARIABLE		= 0x02020001,
				PLT_STATE_VAR_ID		= 0x02020002,
				PLT_STATE_VAR_DATA		= 0x02020003,
				PLT_GLOBAL_DATA			= 0x02020100,
//				PLT_MATERIAL_DATA		= 0x02020200,		// this was removed
				PLT_NODE_DATA			= 0x02020300,
				PLT_ELEMENT_DATA		= 0x02020400,
				PLT_FACE_DATA			= 0x02020500,
			PLT_MESH_STATE				= 0x02030000,
				PLT_ELEMENT_STATE		= 0x02030001,
			PLT_OBJECTS_STATE			= 0x02040000
	};
	// --- element types ---
	enum Elem_Type { 
		PLT_ELEM_HEX, 
		PLT_ELEM_PENTA, 
		PLT_ELEM_TET4, 
		PLT_ELEM_QUAD, 
		PLT_ELEM_TRI, 
		PLT_ELEM_LINE2, 
		PLT_ELEM_HEX20, 
		PLT_ELEM_TET10, 
		PLT_ELEM_TET15, 
		PLT_ELEM_HEX27,
        PLT_ELEM_TRI6,
        PLT_ELEM_QUAD8,
        PLT_ELEM_QUAD9,
        PLT_ELEM_PENTA15,
		PLT_ELEM_TET20,
		PLT_ELEM_TRI10,
		PLT_ELEM_PYRA5,
		PLT_ELEM_TET5,
        PLT_ELEM_PYRA13
    };

	// size of name variables
	enum { STR_SIZE = 64 };


public:
	// Dictionary entry
	class DICTIONARY_ITEM
	{
	public:
		DICTIONARY_ITEM();
		DICTIONARY_ITEM(const DICTIONARY_ITEM& item);

	public:
		FEPlotData*		m_psave;
		unsigned int	m_ntype;	// data type
		unsigned int	m_nfmt;		// storage format
		unsigned int	m_arraySize;	// size of arrays (only used by arrays)
		std::vector<string>	m_arrayNames;	// names of array components (optional)
		char			m_szname[STR_SIZE];
		char			m_szunit[STR_SIZE];
	};

	class Dictionary
	{
	public:
		bool AddVariable(FEModel* pfem, const char* szname, std::vector<int>& item, const char* szdom = "");

		int NodalVariables() { return (int)m_Node.size(); }
		int DomainVarialbes() { return (int)m_Elem.size(); }
		int SurfaceVariables() { return (int)m_Face.size(); }

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
		bool AddNodalVariable   (FEPlotData* ps, const char* szname, std::vector<int>& item);
		bool AddDomainVariable  (FEPlotData* ps, const char* szname, std::vector<int>& item);
		bool AddSurfaceVariable (FEPlotData* ps, const char* szname, std::vector<int>& item);

	protected:
		list<DICTIONARY_ITEM>	m_Glob;		// Global variables
		list<DICTIONARY_ITEM>	m_Mat;		// Material variables
		list<DICTIONARY_ITEM>	m_Node;		// Node variables
		list<DICTIONARY_ITEM>	m_Elem;		// Domain variables
		list<DICTIONARY_ITEM>	m_Face;		// Surface variables

		friend class FEBioPlotFile;
	};

	struct Surface
	{
		int			maxNodes;
		FESurface*	surf;
	};

	class PlotObject
	{
	public:
		PlotObject() {}

		void AddData(const char* szname, Var_Type type, FEPlotData* psave = nullptr);

	public:
		int		m_id;	// object ID
		int		m_tag;	// user tag

		vec3d	m_pos;	// object's position
		quatd	m_rot;	// object's orientation

		std::string	m_name;	// object's name

		list<DICTIONARY_ITEM>	m_data;
	};

	class PointObject : public PlotObject
	{
	public:
		PointObject() {}

	public:
		vec3d	m_r;	// point position
	};

	class LineObject : public PlotObject
	{
	public:
		LineObject() {}

	public:
		vec3d	m_r1;	// point 1
		vec3d	m_r2;	// point 2
	};

public:
	FEBioPlotFile(FEModel* fem);
	~FEBioPlotFile(void);

	//! Open the plot database
	bool Open(const char* szfile) override;

	//! Close the plot database
	void Close() override;

	//! Open for appending
	bool Append(const char* szfile) override;

	//! Write current FE state to plot database
	bool Write(float ftime, int flag = 0)  override;

	//! see if the plot file is valid
	bool IsValid() const override;

public:
	//! Add a variable to the dictionary
	bool AddVariable(FEPlotData* ps, const char* szname);
	bool AddVariable(const char* sz);
	bool AddVariable(const char* sz, std::vector<int>& item, const char* szdom = "");

	//! Set the compression level
	void SetCompression(int n);

	// Write a mesh section
	bool WriteMeshSection(FEModel& fem);

	//! set the software variable
	void SetSoftwareString(const std::string& softwareString);

public:
	int PointObjects();
	PointObject* GetPointObject(int i);
	PointObject* AddPointObject(const std::string& name);

	int LineObjects();
	LineObject* GetLineObject(int i);
	LineObject* AddLineObject(const std::string& name);

public:
	const Dictionary& GetDictionary() const { return m_dic; }

protected:
	bool WriteRoot      (FEModel& fem);
	bool WriteHeader    (FEModel& fem);
	bool WriteDictionary(FEModel& fem);

	void WriteDicList(list<DICTIONARY_ITEM>& dic);
	void WriteDictionaryItem(DICTIONARY_ITEM& it);

	void WriteNodeSection   (FEMesh& m);
	void WriteDomainSection (FEMesh& m);
	void WriteSurfaceSection(FEMesh& m);
	void WriteNodeSetSection(FEMesh& m);
	void WriteElementSetSection(FEMesh& m);
	void WriteFacetSetSection(FEMesh& m);
	void WritePartsSection  (FEModel& fem);
	void WriteObjectsSection();
	void WriteObject(PlotObject* po);

	void WriteSolidDomain   (FESolidDomain&    dom);
	void WriteShellDomain   (FEShellDomain&    dom);
	void WriteBeamDomain    (FEBeamDomain&    dom);
	void WriteDiscreteDomain(FEDiscreteDomain& dom);
    void WriteDomain2D      (FEDomain2D&       dom);

	void WriteGlobalData  (FEModel& fem);
	void WriteNodeData    (FEModel& fem);
	void WriteDomainData  (FEModel& fem);
	void WriteSurfaceData (FEModel& fem);
	void WriteObjectsState();
	void WriteObjectData(PlotObject* po);

	void WriteNodeDataField(FEModel& fem, FEPlotData* pd);
	void WriteDomainDataField(FEModel& fem, FEPlotData* pd);
	void WriteSurfaceDataField(FEModel& fem, FEPlotData* pd);

	void WriteMeshState(FEMesh& mesh);

protected:
	bool ReadDictionary();
	bool ReadDicList();
	void BuildSurfaceTable();

protected:
	Dictionary	m_dic;	// dictionary
	PltArchive	m_ar;	// the data archive
	int			m_ncompress;	// compression level
	int			m_meshesWritten;	// nr of meshes written
	string		m_softwareString;	// the software string
	bool		m_exportUnitsFlag;	// flag that indicates whether to write units

	std::vector<Surface>	m_Surf;

	std::vector<PointObject*>	m_Points;
	std::vector<LineObject*>		m_Lines;
};

//-----------------------------------------------------------------------------
class FEPlotObjectData : public FEPlotData
{
	FECORE_BASE_CLASS(FEPlotObjectData)

public:
	FEPlotObjectData(FEModel* fem) : FEPlotData(fem) {}

	virtual bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar) = 0;
};
