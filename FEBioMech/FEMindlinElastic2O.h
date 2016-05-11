#pragma once
#include "FEElasticMaterial2O.h"

class FEMindlinElastic2O : public FEElasticMaterial2O
{
public:
	FEMindlinElastic2O(FEModel* pfem);

	//! Calculate PK1 stress and higher order stress Q
	void Stress(FEMaterialPoint& mp, mat3d& P, tens3drs& Q);

	//! Calculate material tangents
	//! C = dP/dF
	//! L = dP/dG
	//! H = dQ/dF
	//! J = dQ/dG
	void Tangent(FEMaterialPoint& mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J);

public:
	double	m_lam;
	double	m_mu;
	double	m_a1;
	double	m_a2;
	double	m_a3;
	double	m_a4;
	double	m_a5;

	DECLARE_PARAMETER_LIST();
};
