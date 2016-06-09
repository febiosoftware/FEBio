#include "stdafx.h"
#include "FEElasticTrussDomain.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
//! Constructor
FEElasticTrussDomain::FEElasticTrussDomain(FEModel* pfem) : FETrussDomain(&pfem->GetMesh()), FEElasticDomain(pfem)
{
	m_pMat = 0;
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FETrussMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		FETrussElement& el = m_Elem[i];
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->Init();
	}
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::UnpackLM(FEElement &el, vector<int>& lm)
{
	lm.resize(6);
	FENode& n1 = m_pMesh->Node(el.m_node[0]);
	FENode& n2 = m_pMesh->Node(el.m_node[1]);
	lm[0] = n1.m_ID[m_dofX];
	lm[1] = n1.m_ID[m_dofY];
	lm[2] = n1.m_ID[m_dofZ];
	lm[3] = n2.m_ID[m_dofX];
	lm[4] = n2.m_ID[m_dofY];
	lm[5] = n2.m_ID[m_dofZ];
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[m_dofX] = DOF_ACTIVE;
				node.m_ID[m_dofY] = DOF_ACTIVE;
				node.m_ID[m_dofZ] = DOF_ACTIVE;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FETrussElement& el = m_Elem[i];
		el.GetMaterialPoint(0)->Update(timeInfo);
	}
}

//-----------------------------------------------------------------------------

void FEElasticTrussDomain::StiffnessMatrix(FESolver* psolver)
{
	matrix ke;
	int NT = m_Elem.size();
	vector<int> lm;
	for (int iel =0; iel<NT; ++iel)
	{
		FETrussElement& el = m_Elem[iel];
		ElementStiffness(iel, ke);
		UnpackLM(el, lm);
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::ElementStiffness(int iel, matrix& ke)
{
	FETrussElement& el = Element(iel);

	// nodal coordinates
	vec3d r0[2], rt[2];
	for (int i=0; i<2; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	// intial length
	double L = (r0[1] - r0[0]).norm();

	// current length
	double l = (rt[1] - rt[0]).norm();

	// get the elastic tangent
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FETrussMaterialPoint& pt = *mp.ExtractData<FETrussMaterialPoint>();
	double E = m_pMat->Tangent(pt);

	// element initial volume
	double V = L*el.m_a0;

	// Kirchhoff Stress
	double tau = pt.m_tau;

	// scalar stiffness
	double k = V / (l*l)*( E - 2*tau);

	// axial force T = s*a = t*V/l
	double T = tau*V/l;

	// element normal
	vec3d n = TrussNormal(el);

	// calculate the tangent matrix
	ke.resize(6, 6);

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

//----------------------------------------------------------------------------
void FEElasticTrussDomain::InternalForces(FEGlobalVector& R)
{
	// element force vector
	vector<double> fe;
	vector<int> lm;
	int NT = m_Elem.size();
	for (int i=0; i<NT; ++i)
	{
		FETrussElement& el = m_Elem[i];
		ElementInternalForces(el, fe);
		UnpackLM(el, lm);
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::ElementInternalForces(FETrussElement& el, vector<double>& fe)
{
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

	// get the element's normal
	vec3d n = TrussNormal(el);

	// get the element's Kirchhoff stress
	double tau = pt.m_tau;

	// nodal coordinates
	vec3d r0[2], rt[2];
	for (int i=0; i<2; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	// initial length
	double L = (r0[1] - r0[0]).norm();

	// current length
	double l = (rt[1] - rt[0]).norm();

	// elements initial volume
	double V = L*el.m_a0;

	// calculate nodal forces
	fe.resize(6);
	fe[0] = tau*V/l*n.x;
	fe[1] = tau*V/l*n.y;
	fe[2] = tau*V/l*n.z;
	fe[3] = -fe[0];
	fe[4] = -fe[1];
	fe[5] = -fe[2];
}

//-----------------------------------------------------------------------------
//! Update the truss' stresses
void FEElasticTrussDomain::Update(const FETimeInfo& tp)
{
	// loop over all elements
	vec3d r0[2], rt[2];
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// unpack the element
		FETrussElement& el = m_Elem[i];

		// setup the material point
		FEMaterialPoint& mp = *(el.GetMaterialPoint(0));
		FETrussMaterialPoint& pt = *(mp.ExtractData<FETrussMaterialPoint>());

		// nodal coordinates
		for (int j=0; j<2; ++j)
		{
			r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
			rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
		}

		double l = (rt[1] - rt[0]).norm();
		double L = (r0[1] - r0[0]).norm();

		// calculate strain
		pt.m_l = l / L;

		// calculate stress
		pt.m_tau = m_pMat->Stress(pt);
	}
}
