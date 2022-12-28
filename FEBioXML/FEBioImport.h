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
#include "FileImport.h"
#include <XML/XMLReader.h>
#include "FECore/FEAnalysis.h"
#include "FECore/FESolver.h"
#include "FECore/DataStore.h"
#include <FECore/FEMesh.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/tens3d.h>
#include <string>

class FENodeSet;
class FEBioImport;

//-----------------------------------------------------------------------------
//! Base class for FEBio (feb) file sections.
class FEBioFileSection : public FEFileSection
{
public:
	FEBioFileSection(FEBioImport* feb);

	FEBioImport* GetFEBioImport();
};

//=============================================================================
//! Implements a class to import FEBio input files
//!
class FEBIOXML_API FEBioImport : public FEFileImport
{
public:
	// invalid version
	class InvalidVersion : public FEFileException
	{
	public: InvalidVersion();
	};

	// invalid material defintion
	class InvalidMaterial : public FEFileException
	{ 
	public: InvalidMaterial(int nel);
	};

	// invalid domain type
	class InvalidDomainType : public FEFileException
	{
	public: InvalidDomainType();
	};

	// cannot create domain
	class FailedCreatingDomain : public FEFileException
	{
	public: FailedCreatingDomain();
	};

	// invalid element type
	class InvalidElementType : public FEFileException
	{
	public: InvalidElementType();
	};

	// failed loading plugin
	class FailedLoadingPlugin : public FEFileException
	{
	public: FailedLoadingPlugin(const char* sz);
	};

	// duplicate material section
	class DuplicateMaterialSection : public FEFileException
	{
	public: DuplicateMaterialSection();
	};

	// invalid domain material
	class InvalidDomainMaterial : public FEFileException
	{ 
	public: InvalidDomainMaterial();
	};

	// missing property
	class MissingProperty : public FEFileException
	{
	public: MissingProperty(const std::string& matName, const char* szprop);
	};

	//! Failed allocating solver
	class FailedAllocatingSolver : public FEFileException
	{
	public: FailedAllocatingSolver(const char* sztype);
	};

	//! error in data generation
	class DataGeneratorError : public FEFileException
	{
	public: 
		DataGeneratorError();
	};

	//! failed building a part
	class FailedBuildingPart : public FEFileException
	{
	public:
		FailedBuildingPart(const std::string& partName);
	};

	// Error while reading mesh data section
	class MeshDataError : public FEFileException
	{
	public:
		MeshDataError();
	};

	// repeated node set
	class RepeatedNodeSet : public FEFileException
	{
	public: RepeatedNodeSet(const std::string& name);
	};

	// repeated surface
	class RepeatedSurface : public FEFileException
	{
	public: RepeatedSurface(const std::string& name);
	};

	// repeated edge set
	class RepeatedEdgeSet : public FEFileException
	{
	public: RepeatedEdgeSet(const std::string& name);
	};

	// repeated element set
	class RepeatedElementSet : public FEFileException
	{
	public: RepeatedElementSet(const std::string& name);
	};

public:
	//! constructor
	FEBioImport();

	//! destructor
	~FEBioImport();

	//! open the file
	bool Load(FEModel& fem, const char* szfile);

	//! read the contents of a file
	bool ReadFile(const char* szfile, bool broot = true);

	//! set a custom model builder (takes ownership of modelBuilder)
	void SetModelBuilder(FEModelBuilder* modelBuilder);

public:
	void SetDumpfileName(const char* sz);
	void SetLogfileName (const char* sz);
	void SetPlotfileName(const char* sz);

	void AddDataRecord(DataRecord* pd);

public:
	// Helper functions for reading node sets, surfaces, etc.
	FENodeSet* ParseNodeSet(XMLTag& tag, const char* szatt = "set");
	FESurface* ParseSurface(XMLTag& tag, const char* szatt = "surf");

protected:
	void ParseVersion(XMLTag& tag);

	void BuildFileSectionMap(int nversion);

public:
	char	m_szdmp[512];
	char	m_szlog[512];
	char	m_szplt[512];

public:
	std::vector<DataRecord*>		m_data;
};
