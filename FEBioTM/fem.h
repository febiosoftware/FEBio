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

private:
	CTask*	m_pTask;
};
