// FERestartImport.h: interface for the FERestartImport class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERESTARTIMPORT_H__A5A88D72_026C_45F5_BECB_5B3C7B3C767C__INCLUDED_)
#define AFX_FERESTARTIMPORT_H__A5A88D72_026C_45F5_BECB_5B3C7B3C767C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

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

#endif // !defined(AFX_FERESTARTIMPORT_H__A5A88D72_026C_45F5_BECB_5B3C7B3C767C__INCLUDED_)
