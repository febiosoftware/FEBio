// FEException.h: interface for the FEException class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEEXCEPTION_H__A56BADF1_E3BA_4482_AE23_EA40A591ED88__INCLUDED_)
#define AFX_FEEXCEPTION_H__A56BADF1_E3BA_4482_AE23_EA40A591ED88__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "fecore_api.h"
#include <vector>
using namespace std;

class FEElement;

class FECORE_API FEException
{
public:
	FEException();
	virtual ~FEException();

};

class FECORE_API NegativeJacobian : public FEException
{
public:
	NegativeJacobian(int iel, int ng, double vol, FEElement* pe = 0);

	// print a message to the screen and log file
	void print();

	int		m_iel;	// element where the jacobian was negative
	int		m_ng;	// integration point
	double	m_vol;	// volume
	FEElement*	m_pel;	// pointer to element

public:
	static bool m_boutput;	//!< set to false to suppress output of negative jacobians
};

class FECORE_API ZeroDiagonal : public FEException
{
private:
	struct EQUATION
	{
		int	node;	// node
		int	dof;	// degree of node
	};

public:
	ZeroDiagonal(int node, int dof);

	char m_szerr[256];	// the error message
};

class FECORE_API EnergyDiverging : public FEException {};

class FECORE_API MaxStiffnessReformations : public FEException {};

class FECORE_API ZeroLinestepSize : public FEException {};

class FECORE_API ForceConversion {};

class FECORE_API IterationFailure {};

class FECORE_API MaxResidualError {};

class FECORE_API NANDetected {};

class FECORE_API FatalError {};

class FECORE_API FEMultiScaleException
{
public:
	FEMultiScaleException(int eid, int gpt) : elemId(eid), gptIndex(gpt) {}

public:
	int elemId;
	int gptIndex;
};

class FECORE_API LinearSolverFailed {};

class FECORE_API DoRunningRestart{};


#endif // !defined(AFX_FEEXCEPTION_H__A56BADF1_E3BA_4482_AE23_EA40A591ED88__INCLUDED_)
