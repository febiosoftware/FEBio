#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! domain class for uniform-deformation-gradient hex elements (UDG)
class FEUDGHexDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEUDGHexDomain(FEModel* pfem);

	//! initialize class
	bool Init() override;

public:
	//! calculates the residual
	void InternalForces(FEGlobalVector& R) override;

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver) override;

	// update domain data
	void Update(const FETimeInfo& tp) override;

protected: // element residual contributions
	//! Calculates the internal stress vector for enhanced strain hex elements
	void UDGInternalForces(FESolidElement& el, vector<double>& fe);

	//! calculates hourglass forces for the UDG element
	void UDGHourglassForces(FESolidElement& el, vector<double>& fe);

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
	double	m_hg;	//!< hourglass parameter
};
