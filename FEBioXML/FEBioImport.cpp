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



#include "stdafx.h"
#include "FEBioImport.h"
#include "FEBioIncludeSection.h"
#include "FEBioModuleSection.h"
#include "FEBioControlSection.h"
#include "FEBioControlSection3.h"
#include "FEBioControlSection4.h"
#include "FEBioGlobalsSection.h"
#include "FEBioMaterialSection.h"
#include "FEBioGeometrySection.h"
#include "FEBioBoundarySection.h"
#include "FEBioBoundarySection3.h"
#include "FEBioCodeSection.h"
#include "FEBioLoadsSection.h"
#include "FEBioContactSection.h"
#include "FEBioConstraintsSection.h"
#include "FEBioInitialSection.h"
#include "FEBioInitialSection3.h"
#include "FEBioLoadDataSection.h"
#include "FEBioOutputSection.h"
#include "FEBioStepSection.h"
#include "FEBioStepSection4.h"
#include "FEBioDiscreteSection.h"
#include "FEBioMeshDataSection.h"
#include "FEBioCodeSection.h"
#include "FEBioRigidSection.h"
#include "FEBioRigidSection4.h"
#include "FEBioMeshAdaptorSection.h"
#include "FEBioMeshSection.h"
#include "FEBioMeshSection4.h"
#include "FEBioMeshDomainsSection4.h"
#include "FEBioStepSection3.h"
#include "FECore/DataStore.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FESurfaceMap.h>
#include <FECore/FEFunction1D.h>
#include <FECore/tens3d.h>
#include "FECore/DOFS.h"
#include <string.h>
#include <stdarg.h>
#include "xmltool.h"

FEBioFileSection::FEBioFileSection(FEBioImport* feb) : FEFileSection(feb) {}

FEBioImport* FEBioFileSection::GetFEBioImport() { return static_cast<FEBioImport*>(GetFileReader()); }

//-----------------------------------------------------------------------------
FEBioImport::InvalidVersion::InvalidVersion()
{
	SetErrorString("Invalid version");
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidMaterial::InvalidMaterial(int nel)
{
	SetErrorString("Element %d has an invalid material type", nel);
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidDomainType::InvalidDomainType()	
{
	SetErrorString("Invalid domain type");
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidDomainMaterial::InvalidDomainMaterial()
{
	SetErrorString("Invalid domain material");
}

//-----------------------------------------------------------------------------
FEBioImport::FailedCreatingDomain::FailedCreatingDomain()
{
	SetErrorString("Failed creating domain");
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidElementType::InvalidElementType()
{
	SetErrorString("Invalid element type\n");
}

//-----------------------------------------------------------------------------
FEBioImport::FailedLoadingPlugin::FailedLoadingPlugin(const char* szfile)
{
	SetErrorString("failed loading plugin %s\n", szfile);
}

//-----------------------------------------------------------------------------
FEBioImport::DuplicateMaterialSection::DuplicateMaterialSection()
{
	SetErrorString("Material section has already been defined");
}

//-----------------------------------------------------------------------------
FEBioImport::MissingProperty::MissingProperty(const std::string& matName, const char* szprop)
{
	SetErrorString("Component \"%s\" needs to have property \"%s\" defined", matName.c_str(), szprop);
}

//-----------------------------------------------------------------------------
FEBioImport::FailedAllocatingSolver::FailedAllocatingSolver(const char* sztype)
{
	SetErrorString("Failed allocating solver \"%s\"", sztype);
}

//-----------------------------------------------------------------------------
FEBioImport::DataGeneratorError::DataGeneratorError()
{
	SetErrorString("Error in data generation");
}

//-----------------------------------------------------------------------------
FEBioImport::FailedBuildingPart::FailedBuildingPart(const std::string& partName)
{
	SetErrorString("Failed building part %s", partName.c_str());
}

//-----------------------------------------------------------------------------
FEBioImport::MeshDataError::MeshDataError()
{
	SetErrorString("An error occurred processing mesh_data section.");
}

FEBioImport::RepeatedNodeSet::RepeatedNodeSet(const std::string& name)
{
	SetErrorString("A nodeset with name \"%s\" was already defined.", name.c_str());
}

FEBioImport::RepeatedSurface::RepeatedSurface(const std::string& name)
{
	SetErrorString("A surface with name \"%s\" was already defined.", name.c_str());
}

FEBioImport::RepeatedEdgeSet::RepeatedEdgeSet(const std::string& name)
{
	SetErrorString("An edge with name \"%s\" was already defined.", name.c_str());
}

FEBioImport::RepeatedElementSet::RepeatedElementSet(const std::string& name)
{
	SetErrorString("An element set with name \"%s\" was already defined.", name.c_str());
}

//-----------------------------------------------------------------------------
FEBioImport::FEBioImport()
{
}

//-----------------------------------------------------------------------------
FEBioImport::~FEBioImport()
{
}

//-----------------------------------------------------------------------------
// Build the file section map based on the version number
void FEBioImport::BuildFileSectionMap(int nversion)
{
	// define the file structure
	m_map["Module"     ] = new FEBioModuleSection     (this);
	m_map["Globals"    ] = new FEBioGlobalsSection    (this);
	m_map["Output"     ] = new FEBioOutputSection     (this);

	// older formats
	if (nversion < 0x0200)
	{
		m_map["Control"    ] = new FEBioControlSection      (this);
		m_map["Material"   ] = new FEBioMaterialSection     (this);
	    m_map["Geometry"   ] = new FEBioGeometrySection1x   (this);
		m_map["Boundary"   ] = new FEBioBoundarySection1x   (this);
		m_map["Loads"      ] = new FEBioLoadsSection1x      (this);
		m_map["Constraints"] = new FEBioConstraintsSection1x(this);
		m_map["Step"       ] = new FEBioStepSection         (this);
		m_map["Initial"    ] = new FEBioInitialSection      (this);
		m_map["LoadData"   ] = new FEBioLoadDataSection     (this);
	}

	// version 2.0
	if (nversion == 0x0200)
	{
		m_map["Control"    ] = new FEBioControlSection     (this);
		m_map["Material"   ] = new FEBioMaterialSection    (this);
	    m_map["Geometry"   ] = new FEBioGeometrySection2   (this);
		m_map["Initial"    ] = new FEBioInitialSection     (this);
		m_map["Boundary"   ] = new FEBioBoundarySection2   (this);
		m_map["Loads"      ] = new FEBioLoadsSection2      (this);
		m_map["Include"    ] = new FEBioIncludeSection     (this);
		m_map["Contact"    ] = new FEBioContactSection2    (this);
		m_map["Discrete"   ] = new FEBioDiscreteSection    (this);
		m_map["Code"       ] = new FEBioCodeSection        (this); // added in FEBio 2.4 (experimental feature!)
		m_map["Constraints"] = new FEBioConstraintsSection2(this);
		m_map["Step"       ] = new FEBioStepSection2       (this);
		m_map["LoadData"   ] = new FEBioLoadDataSection    (this);
	}

	// version 2.5
	if (nversion == 0x0205)
	{
		m_map["Control"    ] = new FEBioControlSection      (this);
		m_map["Material"   ] = new FEBioMaterialSection     (this);
	    m_map["Geometry"   ] = new FEBioGeometrySection25   (this);
		m_map["Include"    ] = new FEBioIncludeSection      (this);
		m_map["Initial"    ] = new FEBioInitialSection25    (this);
		m_map["Boundary"   ] = new FEBioBoundarySection25   (this);
		m_map["Loads"      ] = new FEBioLoadsSection25      (this);
		m_map["Contact"    ] = new FEBioContactSection25    (this);
		m_map["Discrete"   ] = new FEBioDiscreteSection25   (this);
		m_map["Constraints"] = new FEBioConstraintsSection25(this);
		m_map["Code"       ] = new FEBioCodeSection         (this); // added in FEBio 2.4 (experimental feature!)
		m_map["LoadData"   ] = new FEBioLoadDataSection     (this);
		m_map["MeshData"   ] = new FEBioMeshDataSection     (this);
		m_map["Step"       ] = new FEBioStepSection25       (this);
		m_map["MeshAdaptor"] = new FEBioMeshAdaptorSection  (this);	// added in FEBio 3.0
	}

	// version 3.0
	if (nversion == 0x0300)
	{
		// we no longer allow unknown attributes
		SetStopOnUnknownAttribute(true);

		m_map["Control"    ] = new FEBioControlSection3     (this);
		m_map["Material"   ] = new FEBioMaterialSection3    (this);
		m_map["Geometry"   ] = new FEBioGeometrySection3    (this);
		m_map["Mesh"       ] = new FEBioMeshSection         (this);
		m_map["MeshDomains"] = new FEBioMeshDomainsSection  (this);
		m_map["Include"    ] = new FEBioIncludeSection      (this);
		m_map["Initial"    ] = new FEBioInitialSection3     (this);
		m_map["Boundary"   ] = new FEBioBoundarySection3    (this);
		m_map["Loads"      ] = new FEBioLoadsSection3       (this);
		m_map["Contact"    ] = new FEBioContactSection25    (this);
		m_map["Discrete"   ] = new FEBioDiscreteSection25   (this);
		m_map["Constraints"] = new FEBioConstraintsSection25(this);
		m_map["Code"       ] = new FEBioCodeSection         (this); // added in FEBio 2.4 (experimental feature!)
		m_map["MeshData"   ] = new FEBioMeshDataSection3    (this);
		m_map["LoadData"   ] = new FEBioLoadDataSection3    (this);
		m_map["Rigid"      ] = new FEBioRigidSection        (this); // added in FEBio 3.0
		m_map["Step"       ] = new FEBioStepSection3        (this);
		m_map["MeshAdaptor"] = new FEBioMeshAdaptorSection  (this);	// added in FEBio 3.0
	}

	// version 4.0
	if (nversion == 0x0400)
	{
		// we no longer allow unknown attributes
		SetStopOnUnknownAttribute(true);

		m_map["Control"    ] = new FEBioControlSection4     (this);
		m_map["Material"   ] = new FEBioMaterialSection3    (this);
		m_map["Mesh"       ] = new FEBioMeshSection4        (this);
		m_map["MeshDomains"] = new FEBioMeshDomainsSection4 (this);
		m_map["Include"    ] = new FEBioIncludeSection      (this);
		m_map["Initial"    ] = new FEBioInitialSection3     (this);
		m_map["Boundary"   ] = new FEBioBoundarySection3    (this);
		m_map["Loads"      ] = new FEBioLoadsSection3       (this);
		m_map["Contact"    ] = new FEBioContactSection4     (this);
		m_map["Discrete"   ] = new FEBioDiscreteSection25   (this);
		m_map["Constraints"] = new FEBioConstraintsSection25(this);
		m_map["Code"       ] = new FEBioCodeSection         (this); // added in FEBio 2.4 (experimental feature!)
		m_map["MeshData"   ] = new FEBioMeshDataSection4    (this);	// added in febio4
		m_map["LoadData"   ] = new FEBioLoadDataSection3    (this);
		m_map["Rigid"      ] = new FEBioRigidSection4       (this); // added in FEBio 4.0
		m_map["Step"       ] = new FEBioStepSection4        (this);
		m_map["MeshAdaptor"] = new FEBioMeshAdaptorSection  (this);	// added in FEBio 3.0
	}
}

//-----------------------------------------------------------------------------
bool FEBioImport::Load(FEModel& fem, const char* szfile)
{
	if (m_builder == nullptr)
	{
		m_builder = new FEModelBuilder(fem);
	}

	// intialize some variables
	m_szdmp[0] = 0;
	m_szlog[0] = 0;
	m_szplt[0] = 0;

	m_data.clear();

	// extract the path
	strcpy(m_szpath, szfile);
	char* ch = strrchr(m_szpath, '\\');
	if (ch==0) ch = strrchr(m_szpath, '/');
	if (ch==0) m_szpath[0] = 0; else *(ch+1)=0;

	// clean up
	fem.GetMesh().ClearDataMaps();

	// read the file
	if (ReadFile(szfile) == false) return false;

	// finish building
	try {
		bool b = m_builder->Finish();
		if (b == false) return errf("FAILED building FEBio model.");
	}
	catch (std::exception e)
	{
		const char* szerr = e.what();
		if (szerr == nullptr) szerr = "(unknown exception)";
		return errf("%s", e.what());
	}	
	catch (...)
	{
		return errf("unknown exception.");
	}

	return true;
}

//-----------------------------------------------------------------------------
//! set a custom model builder 
void FEBioImport::SetModelBuilder(FEModelBuilder* modelBuilder)
{
	delete m_builder;
	m_builder = modelBuilder;
}

//-----------------------------------------------------------------------------
// This function parses the XML input file. The broot parameter is used to indicate
// if this is the main control file or an included file. 
bool FEBioImport::ReadFile(const char* szfile, bool broot)
{
	// Open the XML file
	XMLReader xml;
	if (xml.Open(szfile) == false) return errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);

	// Find the root element
	XMLTag tag;
	try
	{
		if (xml.FindTag("febio_spec", tag) == false) return errf("FATAL ERROR: febio_spec tag was not found. This is not a valid input file.\n\n");
	}
	catch (...)
	{
		return errf("An error occured while finding the febio_spec tag.\nIs this a valid FEBio input file?\n\n");
	}

	// parse the file
	try
	{
		// get the version number
		ParseVersion(tag);

		// FEBio4 only supports file version 1.2, 2.0, 2.5, 3.0, and 4.0
		int nversion = GetFileVersion();
		if ((nversion != 0x0102) && 
			(nversion != 0x0200) && 
			(nversion != 0x0205) && 
			(nversion != 0x0300) && 
			(nversion != 0x0400)) throw InvalidVersion();

		// build the file section map based on the version number
		BuildFileSectionMap(nversion);

		// For versions before 2.5 we need to allocate all the degrees of freedom beforehand. 
		// This is necessary because the Module section doesn't have to defined until a Control section appears.
		// That means that model components that depend on DOFs can be defined before the Module tag (e.g. in multi-step analyses) and this leads to problems.
		// In 2.5 this is solved by requiring that the Module tag is defined at the top of the file. 
		if (broot && (nversion < 0x0205))
		{
			// We need to define a default Module type since before 2.5 this tag is optional for structural mechanics model definitions.
			GetBuilder()->SetActiveModule("solid");

			// set default variables for older files.
			GetBuilder()->SetDefaultVariables();
		}

		// parse the file
		++tag;

		// From version 2.5 and up the first tag of the (main control) file has to be the Module tag.
		if (broot && (nversion >= 0x0205))
		{
			if (tag != "Module")
			{
				return errf("First tag must be the Module section.\n\n");
			}

			// try to find a section parser
			FEFileSectionMap::iterator is = m_map.find(tag.Name());

			// make sure we found a section reader
			if (is == m_map.end()) throw XMLReader::InvalidTag(tag);

			// parse the module tag
			is->second->Parse(tag);

			// Now that the Module tag is read in, we'll want to create an analysis step.
			// Creating an analysis step will allocate a solver class (based on the module) 
			// and this in turn will allocate the degrees of freedom.
			// TODO: This is kind of a round-about way and I really want to find a better solution.
			// NOTE: For version 4.0 we do not allocate the solver by default
			GetBuilder()->GetStep(nversion >= 0x0400 ? false : true);

			// let's get the next tag
			++tag;
		}

		do
		{
			// try to find a section parser
			FEFileSectionMap::iterator is = m_map.find(tag.Name());

			// make sure we found a section reader
			if (is == m_map.end()) throw XMLReader::InvalidTag(tag);

			// see if the file has the "from" attribute (for version 2.0 and up)
			if (nversion >= 0x0200)
			{
				const char* szinc = tag.AttributeValue("from", true);
				if (szinc)
				{
					// make sure this is a leaf
					if (tag.isleaf() == false) return errf("FATAL ERROR: included sections may not have child sections.\n\n");

					// read this section from an included file.
					XMLReader xml2;
					if (xml2.Open(szinc) == false) return errf("FATAL ERROR: failed opening input file %s\n\n", szinc);

					// find the febio_spec tag
					XMLTag tag2;
					if (xml2.FindTag("febio_spec", tag2) == false) return errf("FATAL ERROR: febio_spec tag was not found. This is not a valid input file.\n\n");

					// find the section we are looking for
					char sz[512] = {0};
					sprintf(sz, "febio_spec/%s", tag.Name());
					if (xml2.FindTag(sz, tag2) == false) return errf("FATAL ERROR: Couldn't find %s section in file %s.\n\n", tag.Name(), szinc);

					// parse the section
					is->second->Parse(tag2);
				}
				else is->second->Parse(tag);
			}
			else
			{
				// parse the section
				is->second->Parse(tag);
			}

			// go to the next tag
			++tag;
		}
		while (!tag.isend());
	}
	// --- XML Reader Exceptions ---
	catch (XMLReader::Error& e)
	{
		return errf("%s", e.what());
	}
	// --- FEBioImport Exceptions ---
	catch (FEFileException& e)
	{
		return errf("%s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
	}
	// --- Exception from DataStore ---
	catch (UnknownDataField& e)
	{
		return errf("\"%s\" is not a valid field variable name (line %d)\n", e.what(), xml.GetCurrentLine()-1);
	}
	// std::exception
	catch (std::exception e)
	{
		const char* szerr = e.what();
		if (szerr == nullptr) szerr = "(unknown exception)";
		return errf("%s", szerr);
	}
	// --- Unknown exceptions ---
	catch (...)
	{
		return errf("unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	// close the XML file
	xml.Close();

	// we're done!
	return true;
}

//-----------------------------------------------------------------------------
//! This function parses the febio_spec tag for the version number
void FEBioImport::ParseVersion(XMLTag &tag)
{
	const char* szv = tag.AttributeValue("version");
	assert(szv);
	int n1, n2;
	int nr = sscanf(szv, "%d.%d", &n1, &n2);
	if (nr != 2) throw InvalidVersion();
	if ((n1 < 1) || (n1 > 0xFF)) throw InvalidVersion();
	if ((n2 < 0) || (n2 > 0xFF)) throw InvalidVersion();
	int nversion = (n1 << 8) + n2;
	SetFileVerion(nversion);
}

//-----------------------------------------------------------------------------
void FEBioImport::SetDumpfileName(const char* sz) { sprintf(m_szdmp, "%s", sz); }
void FEBioImport::SetLogfileName (const char* sz) { sprintf(m_szlog, "%s", sz); }
void FEBioImport::SetPlotfileName(const char* sz) { sprintf(m_szplt, "%s", sz); }

//-----------------------------------------------------------------------------
void FEBioImport::AddDataRecord(DataRecord* pd)
{
	m_data.push_back(pd);
}

//-----------------------------------------------------------------------------
// This tag parses a node set.
FENodeSet* FEBioImport::ParseNodeSet(XMLTag& tag, const char* szatt)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FENodeSet* pns = 0;

	// see if the set attribute is defined
	const char* szset = tag.AttributeValue(szatt, true);
	if (szset)
	{
		// Make sure this is an empty tag
		if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

		// find the node set
		pns = mesh.FindNodeSet(szset);
		if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, szatt, szset);
	}
	else
	{
		// This defines a node set, so we need a name tag
		// (For now this name is optional)
		const char* szname = tag.AttributeValue("name", true);
		if (szname == 0) szname = "_unnamed";

		// create a new node set
		pns = new FENodeSet(GetFEModel());
		pns->SetName(szname);

		// add the nodeset to the mesh
		mesh.AddNodeSet(pns);

		// read the nodes
		if (tag.isleaf())
		{
			// This format is deprecated
			vector<int> l;
			fexml::readList(tag, l);
			for (int i=0; i<l.size(); ++i) pns->Add(GetBuilder()->FindNodeFromID(l[i]));
		}
		else
		{
			// read the nodes
			++tag;
			do
			{
				if ((tag == "n") || (tag == "node"))
				{
					int nid = -1;
					tag.AttributeValue("id", nid);

					nid = GetBuilder()->FindNodeFromID(nid);
					pns->Add(nid);
				}
				else if (tag == "NodeSet")
				{
					const char* szset = tag.AttributeValue(szatt);
					
					// Make sure this is an empty tag
					if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

					// find the node set
					FENodeSet* ps = mesh.FindNodeSet(szset);
					if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, szatt, szset);

					// add the node set
					pns->Add(ps->GetNodeList());
				}
				else if (tag == "node_list")
				{
					vector<int> nl;
					fexml::readList(tag, nl);
					for (int i = 0; i<nl.size(); ++i) pns->Add(GetBuilder()->FindNodeFromID(nl[i]));
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
	}

	return pns;
}

//-----------------------------------------------------------------------------
FESurface* FEBioImport::ParseSurface(XMLTag& tag, const char* szatt)
{
	FEMesh& m = GetFEModel()->GetMesh();

	// create new surface
	FESurface* psurf = fecore_alloc(FESurface, GetFEModel());

	// see if the surface is referenced by a set of defined explicitly
	const char* szset = tag.AttributeValue(szatt, true);
	if (szset)
	{
		// make sure this tag does not have any children
		if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);

		// see if we can find the facet set
		FEMesh& m = GetFEModel()->GetMesh();
		FEFacetSet* ps = m.FindFacetSet(szset);

		// create a surface from the facet set
		if (ps)
		{
			if (GetBuilder()->BuildSurface(*psurf, *ps) == false) throw XMLReader::InvalidTag(tag);
		}
		else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
	}
	else
	{
		// count how many pressure cards there are
		int npr = tag.children();
		psurf->Create(npr);

		FEModelBuilder* feb = GetBuilder();

		++tag;
		int nf[FEElement::MAX_NODES ], N;
		for (int i=0; i<npr; ++i)
		{
			FESurfaceElement& el = psurf->Element(i);

			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(feb->m_ntri3);
			else if (tag == "tri6" ) el.SetType(feb->m_ntri6);
			else if (tag == "tri7" ) el.SetType(feb->m_ntri7);
			else if (tag == "tri10") el.SetType(feb->m_ntri10);
			else if (tag == "quad8") el.SetType(FE_QUAD8G9);
			else if (tag == "quad9") el.SetType(FE_QUAD9G9);
			else throw XMLReader::InvalidTag(tag);

			N = el.Nodes();
			tag.value(nf, N);
			for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

			++tag;
		}
	}

	return psurf;
}
