#pragma once
#include "FEBioLib/FEBioModel.h"
#include "FEBioPlot/PlotFile.h"
#include "FEBioLib/Timer.h"
#include "FEBioLib/DataStore.h"

class CTask;

//-----------------------------------------------------------------------------
// The FE model class
class FEM : public FEBioModel
{
public:
	// constructor
	FEM();
	FEM(CTask* pt);

public:
	virtual void PushState();
	virtual void PopState ();

protected:
	void ShallowCopy(FEM& fem);

public:
	//! check for user interruption
	virtual void CheckInterruption();

private:
	CTask*	m_pTask;
};
