#pragma once
#include <FEBioXML/XMLReader.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// IO exceptions

//! the variable name is not recognized
class InvalidVariableName
{
public:
	InvalidVariableName(const char* sz);
	char szname[256];
};

//! there is nothing to optimize
class NothingToOptimize{};

//! FEBio error terminated during the optimization
class FEErrorTermination{};

//-----------------------------------------------------------------------------
class FEOptimizeData;

//=============================================================================
//! Class that reads the optimization input file
class FEOptimizeInput
{
public:
	bool Input(const char* szfile, FEOptimizeData* pOpt);

protected:
	bool ParseOptions(XMLTag& tag, FEOptimizeData& opt);
	bool ParseTask(XMLTag& tag, FEOptimizeData& opt);
	bool ParseObjective(XMLTag& tag, FEOptimizeData& opt);
	bool ParseParameters(XMLTag& tag, FEOptimizeData& opt);
	bool ParseConstraints(XMLTag& tag, FEOptimizeData& opt);
	FEDataSource* ParseDataSource(XMLTag& tag, FEOptimizeData& opt);

protected:
	bool ReadParameter(XMLTag& tag, FEParameterList& pl);
};
