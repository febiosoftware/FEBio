#pragma once
#include "FECore/FESurfacePairConstraint.h"

class FEModel;
class FESolver;
class FEGlobalMatrix;

// Macauley bracket
#define MBRACKET(x) ((x)>=0? (x): 0)

// Heavyside function
#define HEAVYSIDE(x) ((x)>=0?1:0)


//-----------------------------------------------------------------------------
//! This is the base class for contact interfaces
class FEContactInterface : public FESurfacePairConstraint
{
public:
	//! constructor
	FEContactInterface(FEModel* pfem);

	//! destructor
	~FEContactInterface();

	//! serialize data to archive
	void Serialize(DumpStream& ar);

protected:
	//! don't call the default constructor
	FEContactInterface() : FESurfacePairConstraint(0){}

	//! auto-penalty calculation
	double AutoPenalty(FESurfaceElement& el, FESurface& s);

public:
	bool	m_blaugon;	//!< augmented lagrangian flag
};
