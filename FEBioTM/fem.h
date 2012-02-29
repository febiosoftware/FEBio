#pragma once
#include "FEBioLib/FEBioModel.h"

//-----------------------------------------------------------------------------
// forward declaration of the task class.
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
	//! copy the model
	void ShallowCopy(FEM& fem);

public:
	//! check for user interruption
	virtual void CheckInterruption();

private:
	CTask*	m_pTask;
};
