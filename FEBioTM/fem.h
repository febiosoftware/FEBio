#pragma once
#include "FEBioLib/FEBioModel.h"
#include "FEBioPlot/PlotFile.h"
#include "FEBioLib/Timer.h"
#include "FEBioLib/DataStore.h"

//-----------------------------------------------------------------------------
// The FE model class
class FEM : public FEBioModel
{
public:
	// constructor
	FEM();

public:
	virtual void PushState();
	virtual void PopState ();

protected:
	void ShallowCopy(FEM& fem);

public:
	//! Evaluate parameter list
	virtual void EvaluateParameterList(FEParameterList& pl) {}
	virtual void EvaluateMaterialParameters(FEMaterial* pm) {}

	//! return a pointer to the named variable
	virtual double* FindParameter(const char* szname) { return 0; }

	//! find a boundary condition from the ID
	virtual FEBoundaryCondition* FindBC(int nid) { return 0; }

	//! check for user interruption
	virtual void CheckInterruption() {}

	//! serialize data
	virtual bool Serialize(DumpFile& ar) { return false; }
};
