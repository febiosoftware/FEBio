#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! domain class for uniform-deformation-gradient hex elements (UDG)
class FEUDGHexDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEUDGHexDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_UDGHEX_DOMAIN; }

	FEDomain* Clone()
	{
		FEUDGHexDomain* pd = new FEUDGHexDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

public:

	//! calculates the residual
	void InternalForces(FENLSolver* psolver, vector<double>& R);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver);

	// update stresses
	void UpdateStresses(FEModel& fem);

protected: // element residual contributions
	//! Calculates the internal stress vector for enhanced strain hex elements
	void UDGInternalForces(FEModel& fem, FESolidElement& el, vector<double>& fe);

	//! calculates hourglass forces for the UDG element
	void UDGHourglassForces(FEModel& fem, FESolidElement& el, vector<double>& fe);

protected: // element stiffness contributions
	//! hourglass stiffness for UDG hex elements
	void UDGHourglassStiffness(FEModel& fem, FESolidElement& el, matrix& ke);

	//! geometrical stiffness for UDG hex elements
	void UDGGeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness for UDG hex elements
	void UDGMaterialStiffness(FESolidElement& el, matrix& ke);

protected:
	void AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int state = 0);
	void AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8]);
	double HexVolume(FESolidElement& el, int state = 0);

public:
	static double	m_hg;	//!< hourglass parameter
};
