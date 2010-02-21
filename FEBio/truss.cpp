#include "stdafx.h"
#include "FESolidSolver.h"
#include "FETrussMaterial.h"
#include "FEDomain.h"

//-----------------------------------------------------------------------------
void FETrussDomain::ElementStiffness(FEM& fem, FETrussElement& el, matrix& ke)
{
	// get the material
	FETrussMaterial* pm = dynamic_cast<FETrussMaterial*>(fem.GetMaterial(el.GetMatID()));
	assert(pm);

	// intial length
	double L = el.Length0();

	// current length
	double l = el.Length();

	// get the elastic tangent
	FEMaterialPoint& mp = *el.m_State[0];
	FETrussMaterialPoint& pt = *mp.ExtractData<FETrussMaterialPoint>();
	double E = pm->Tangent(pt);

	// element initial volume
	double V = el.Volume0();

	// Kirchhoff Stress
	double tau = pt.m_tau;

	// scalar stiffness
	double k = V / (l*l)*( E - 2*tau);

	// axial force T = s*a = t*V/l
	double T = tau*V/l;

	// element normal
	vec3d n = el.Normal();

	// calculate the tangent matrix
	ke.Create(6, 6);

	ke[0][0] = ke[3][3] = k*n.x*n.x + T/l;
	ke[1][1] = ke[4][4] = k*n.y*n.y + T/l;
	ke[2][2] = ke[5][5] = k*n.z*n.z + T/l;

	ke[0][1] = ke[1][0] = ke[3][4] = ke[4][3] = k*n.x*n.y;
	ke[1][2] = ke[2][1] = ke[4][5] = ke[5][4] = k*n.y*n.z;
	ke[0][2] = ke[2][0] = ke[3][5] = ke[5][3] = k*n.x*n.z;

	ke[0][3] = ke[3][0] = -ke[0][0]; ke[0][4] = ke[4][0] = -ke[0][1]; ke[0][5] = ke[5][0] = -ke[0][2];
	ke[1][3] = ke[3][1] = -ke[1][0]; ke[1][4] = ke[4][1] = -ke[1][1]; ke[1][5] = ke[5][1] = -ke[1][2];
	ke[2][3] = ke[3][2] = -ke[2][0]; ke[2][4] = ke[4][2] = -ke[2][1]; ke[2][5] = ke[5][2] = -ke[2][2];
}


//-----------------------------------------------------------------------------
void FETrussDomain::InternalForces(FETrussElement& el, vector<double>& fe)
{
	FEMaterialPoint& mp = *el.m_State[0];
	FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

	// get the element's normal
	vec3d n = el.Normal();

	// get the element's Kirchhoff stress
	double tau = pt.m_tau;

	// elements initial volume
	double V = el.Volume0();

	// current length
	double l = el.Length();

	// calculate nodal forces
	fe.create(6);
	fe[0] = tau*V/l*n.x;
	fe[1] = tau*V/l*n.y;
	fe[2] = tau*V/l*n.z;
	fe[3] = -fe[0];
	fe[4] = -fe[1];
	fe[5] = -fe[2];
}
