// FEDiagnostic.h: interface for the FEDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEDIAGNOSTIC_H__75EB5A08_CE16_45BD_A223_7BD93BF0837A__INCLUDED_)
#define AFX_FEDIAGNOSTIC_H__75EB5A08_CE16_45BD_A223_7BD93BF0837A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEModel.h"
#include "FEBioXML/FEBioImport.h"

//-----------------------------------------------------------------------------
//! The FEDiagnostic class is a base class that can be used to create
//! diagnostic classes to test FEBio's performance.

class FEDiagnostic
{
public:
	//! constructor
	FEDiagnostic(FEModel& fem);

	//! destructor
	virtual ~FEDiagnostic();

	//! initialization
	virtual bool Init() { return true; }

	//! run the diagnostic. Returns true on pass, false on failure
	virtual bool Run() = 0;

	//! load data from file
	virtual bool ParseSection(XMLTag& tag) { return false; }

public:
	FEModel&	m_fem;	//!< the FEModel object the diagnostic is performed on
};

//-----------------------------------------------------------------------------
// Scenario Section parser
class FEBioScenarioSection : public FEBioFileSection
{
public:
	FEBioScenarioSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
//! The FEDiagnosticImport class creates a specific diagnostic test. Currently
//! the only way to create a diagnostic is to load a diagnostic from file

class FEDiagnosticImport : public FEFEBioImport
{
public:
	FEDiagnostic* LoadFile(FEModel& fem, const char* szfile);

protected:
	FEDiagnostic* m_pdia;

	friend class FEBioScenarioSection;
};

#endif // !defined(AFX_FEDIAGNOSTIC_H__75EB5A08_CE16_45BD_A223_7BD93BF0837A__INCLUDED_)
