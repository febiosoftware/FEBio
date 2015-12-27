// FEContactInterface.h: interface for the FEContactInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FECONTACTINTERFACE_H__DB89AB6C_87DA_4C47_923C_B9D488E04403__INCLUDED_)
#define AFX_FECONTACTINTERFACE_H__DB89AB6C_87DA_4C47_923C_B9D488E04403__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEMesh.h"
#include "FECore/DumpFile.h"
#include "FECore/FESurfacePairInteraction.h"

class FEModel;
class FESolver;
class FEStiffnessMatrix;

// Macauley bracket
#define MBRACKET(x) ((x)>=0? (x): 0)

// Heavyside function
#define HEAVYSIDE(x) ((x)>=0?1:0)


//-----------------------------------------------------------------------------
//! This is the base class for contact interfaces
//! \todo Merge the active interface with the boundary condition class
class FEContactInterface : public FESurfacePairInteraction
{
public:
	//! constructor
	FEContactInterface(FEModel* pfem);

	//! destructor
	virtual ~FEContactInterface();

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R) = 0;

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver) = 0;

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug) = 0;

	//! serialize data to archive
	virtual void Serialize(DumpFile& ar);

	//! build the matrix profile for use in the stiffness matrix
	virtual void BuildMatrixProfile(FEGlobalMatrix& K) = 0;

protected:
	//! don't call the default constructor
	FEContactInterface() : FESurfacePairInteraction(0){}

	//! auto-penalty calculation
	double AutoPenalty(FESurfaceElement& el, FESurface& s);

public:
	bool	m_blaugon;	//!< augmented lagrangian flag

};

#endif // !defined(AFX_FECONTACTINTERFACE_H__DB89AB6C_87DA_4C47_923C_B9D488E04403__INCLUDED_)
