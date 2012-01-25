#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FENLSolver.h"

//-----------------------------------------------------------------------------
class FEHeatDomain
{
public:
	virtual ~FEHeatDomain(){}
	virtual void ConductionMatrix (FENLSolver* pnls) = 0;
	virtual void CapacitanceMatrix(FENLSolver* pnls, double dt) = 0;
};

//-----------------------------------------------------------------------------
//! domain class for 3D heat elements
class FEHeatSolidDomain : public FESolidDomain, public FEHeatDomain
{
public:
	//! constructor
	FEHeatSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_HEAT_SOLID_DOMAIN, pm, pmat) {}

	//! Create a clone of this domain
	FEDomain* Clone();

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public: // overloaded from FEHeatDomain

	//! Calculate the conduction stiffness 
	void ConductionMatrix(FENLSolver* psolver);

	//! Calculate capacitance stiffness matrix
	void CapacitanceMatrix(FENLSolver* psolver, double dt);

protected:
	//! calculate the conductive element stiffness matrix
	void ElementConduction(FESolidElement& el, matrix& ke);

	//! calculate the capacitance element stiffness matrix
	void ElementCapacitance(FESolidElement& el, matrix& ke, double dt);
};
