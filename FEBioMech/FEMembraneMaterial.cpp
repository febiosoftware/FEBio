#include "stdafx.h"
#include "FEMembraneMaterial.h"

//-----------------------------------------------------------------------------
// convet the displacement gradient to a strain
void FEMembraneMaterialPoint::strain(double e[3])
{
	const double h1[6] = {1, 0, 0, 0, 0, 0};
	const double h2[6] = {0, 0, 0, 0, 1, 0};
	const double h3[6] = {0, 1, 0, 1, 0, 0};

	const double H1[6][6] = {
		{1, 0, 0, 0, 0, 0},
		{0, 1, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0}};

	const double H2[6][6] = {
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 1, 0},
		{0, 0, 0, 0, 0, 1}};

	const double H3[6][6] = {
		{0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 1, 0},
		{0, 0, 0, 0, 0, 1},
		{1, 0, 0, 0, 0, 0},
		{0, 1, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0}};

	e[0] = e[1] = e[2] = 0.0;
	for (int i=0; i<6; ++i)
	{
		e[0] += h1[i]*g[i];
		e[1] += h2[i]*g[i];
		e[2] += h3[i]*g[i];
		for (int j=0; j<6; ++j)
		{
			e[0] += 0.5*g[i]*H1[i][j]*g[j];
			e[1] += 0.5*g[i]*H2[i][j]*g[j];
			e[2] += 0.5*g[i]*H3[i][j]*g[j];
		}
	}
}

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEElasticMembrane, FEMembraneMaterial)
	ADD_PARAMETER2(m_E, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER2(m_v, FE_PARAM_DOUBLE, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEElasticMembrane::Stress(FEMaterialPoint &mp, double s[])
{
	FEMembraneMaterialPoint& pt = *mp.ExtractData<FEMembraneMaterialPoint>();

	// get the material point strain
	double e[3];
	pt.strain(e);

	// elasticity tangent
	double E = 1.0;
	double v = 0.35;
	double d = 1.0/(1.0 - v*v);
	double D[3][3] = {
		{d*E, d*v*E, 0},
		{d*v*E, d*E, 0},
		{0, 0, d*E*(1.0-v)*0.5}};

	// calculate 2D PK2-stress
	s[0] = D[0][0]*e[0] + D[0][1]*e[1] + D[0][2]*e[2];
	s[1] = D[1][0]*e[0] + D[1][1]*e[1] + D[1][2]*e[2];
	s[2] = D[2][0]*e[0] + D[2][1]*e[1] + D[2][2]*e[2];
}

//-----------------------------------------------------------------------------
void FEElasticMembrane::Tangent(FEMaterialPoint &mp, double D[][3])
{
	FEMembraneMaterialPoint& pt = *mp.ExtractData<FEMembraneMaterialPoint>();

	// elasticity tangent
	double E = 1.0;
	double v = 0.35;
	double d = 1.0/(1.0 - v*v);
	D[0][0] =   d*E; D[0][1] = d*v*E; D[0][2] = 0;
	D[1][0] = d*v*E; D[1][1] =   d*E; D[1][2] = 0;
	D[2][0] =     0; D[2][1] =     0; D[2][2] = d*E*(1.0-v)*0.5;
}
