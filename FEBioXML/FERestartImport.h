#pragma once
#include "FileImport.h"
#include "XMLReader.h"

//-----------------------------------------------------------------------------
class FERestartControlSection : public FEFileSection
{
public:
	FERestartControlSection(FEFileImport* reader) : FEFileSection(reader) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
//! Restart input file reader.
class FEBIOXML_API FERestartImport : public FEFileImport
{
public:
	FERestartImport();
	virtual ~FERestartImport();

	bool Load(FEModel& fem, const char* szfile);

public:
	char		m_szdmp[256];	// user defined restart file name

protected:
	XMLReader	m_xml;			// the file reader
};
