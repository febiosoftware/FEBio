#include "stdafx.h"
#include "FEBioStepSection.h"
#include "FEBioModuleSection.h"
#include "FEBioControlSection.h"
#include "FEBioConstraintsSection.h"
#include "FEBioBoundarySection.h"
#include "FEBioLoadsSection.h"
#include "FEBioContactSection.h"
#include "FEBioInitialSection.h"

//-----------------------------------------------------------------------------
void FEBioStepSection::Parse(XMLTag& tag)
{
	// create next step
	FEBioImport* feb = GetFEBioImport();
	feb->GetBuilder()->NextStep();

	FEFileSectionMap Map;
	Map["Module"     ] = new FEBioModuleSection       (feb);
	Map["Control"    ] = new FEBioControlSection      (feb);
	Map["Constraints"] = new FEBioConstraintsSection1x(feb);
	Map["Boundary"   ] = new FEBioBoundarySection1x   (feb);
	Map["Loads"      ] = new FEBioLoadsSection        (feb);
	Map["Initial"    ] = new FEBioInitialSection      (feb);

	// parse the file sections
	Map.Parse(tag);
}

//-----------------------------------------------------------------------------
void FEBioStepSection2::Parse(XMLTag& tag)
{
	// create next step
	FEBioImport* feb = GetFEBioImport();
	feb->GetBuilder()->NextStep();

	FEFileSectionMap Map;
	Map["Module"     ] = new FEBioModuleSection      (feb);
	Map["Control"    ] = new FEBioControlSection     (feb);
	Map["Constraints"] = new FEBioConstraintsSection2(feb);
	Map["Boundary"   ] = new FEBioBoundarySection2   (feb);
	Map["Loads"      ] = new FEBioLoadsSection       (feb);
	Map["Initial"    ] = new FEBioInitialSection     (feb);
	Map["Contact"    ] = new FEBioContactSection     (feb);

	// parse the file sections
	Map.Parse(tag);
}

//-----------------------------------------------------------------------------
void FEBioStepSection25::Parse(XMLTag& tag)
{
	// get the (optional) type attribute
	const char* sztype = tag.AttributeValue("type", true);
	if (sztype)
	{
		GetBuilder()->SetModuleName(sztype);
	}

	// create next step
	GetBuilder()->NextStep();

	FEFileImport* imp = GetFileReader();
	FEFileSectionMap Map;
	Map["Control"    ] = new FEStepControlSection     (imp);
	Map["Constraints"] = new FEBioConstraintsSection25(imp);
	Map["Boundary"   ] = new FEBioBoundarySection25   (imp);
	Map["Loads"      ] = new FEBioLoadsSection        (imp);
	Map["Initial"    ] = new FEBioInitialSection25    (imp);
	Map["Contact"    ] = new FEBioContactSection      (imp);

	// parse the file sections
	Map.Parse(tag);
}
