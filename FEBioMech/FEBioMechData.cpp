#include "stdafx.h"
#include "FEBioMechData.h"
#include "FEElasticMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FEDamageMaterialPoint.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FERigidMaterial.h"
#include "FESolidSolver.h"
#include "FESolidSolver2.h"
#include "FERigidBody.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FERigidConnector.h"
#include "FEVolumeConstraint.h"

//-----------------------------------------------------------------------------
double FENodeXPos::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.x; 
}

//-----------------------------------------------------------------------------
double FENodeYPos::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.y; 
}

//-----------------------------------------------------------------------------
double FENodeZPos::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.z; 
}

//-----------------------------------------------------------------------------
double FENodeXDisp::value(int nnode) 
{
	const int dof_X = m_pfem->GetDOFIndex("x");
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.get(dof_X); 
}

//-----------------------------------------------------------------------------
double FENodeYDisp::value(int nnode) 
{
	const int dof_Y = m_pfem->GetDOFIndex("y");
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.get(dof_Y); 
}

//-----------------------------------------------------------------------------
double FENodeZDisp::value(int nnode) 
{
	const int dof_Z = m_pfem->GetDOFIndex("z");
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.get(dof_Z); 
}

//-----------------------------------------------------------------------------
double FENodeXVel::value(int nnode) 
{
	const int dof_VX = m_pfem->GetDOFIndex("vx");
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.get(dof_VX);
}

//-----------------------------------------------------------------------------
double FENodeYVel::value(int nnode) 
{
	const int dof_VY = m_pfem->GetDOFIndex("vy");
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.get(dof_VY); 
}

//-----------------------------------------------------------------------------
double FENodeZVel::value(int nnode) 
{
	const int dof_VZ = m_pfem->GetDOFIndex("vz");
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.get(dof_VZ);
}

//-----------------------------------------------------------------------------
double FENodeXAcc::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_at.x; 
}

//-----------------------------------------------------------------------------
double FENodeYAcc::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_at.y; 
}

//-----------------------------------------------------------------------------
double FENodeZAcc::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_at.z; 
}

//-----------------------------------------------------------------------------
double FENodeForceX::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FESolidSolver2* psolid_solver = dynamic_cast<FESolidSolver2*>(m_pfem->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		vector<int>& id = mesh.Node(nnode).m_ID;
		return (-id[0] - 2 >= 0 ? Fr[-id[0]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceY::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FESolidSolver2* psolid_solver = dynamic_cast<FESolidSolver2*>(m_pfem->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		vector<int>& id = mesh.Node(nnode).m_ID;
		return (-id[1] - 2 >= 0 ? Fr[-id[1]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceZ::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FESolidSolver2* psolid_solver = dynamic_cast<FESolidSolver2*>(m_pfem->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		vector<int>& id = mesh.Node(nnode).m_ID;
		return (-id[2] - 2 >= 0 ? Fr[-id[2]-2] : 0);
	}
	return 0;
}


//-----------------------------------------------------------------------------
double FELogElemPosX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_rt.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPosY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_rt.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPosZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_rt.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemJacobian::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_J;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xx();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.yy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.zz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.yz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain1::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain2::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[1];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain3::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[2];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.xx();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.yy();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.zz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.xy();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.yz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.xz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xx();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.yy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.zz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.yz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress1::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress2::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[1];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress3::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[2];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemElasticity_::value(FEElement& el, int n)
{
    FEElasticMaterial* pme = m_pfem->GetMaterial(el.GetMatID())->GetElasticMaterial();
    if ((pme == 0) || pme->IsRigid()) return 0;

    tens4ds c;
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
        c = pme->Tangent(pt);
		val += c.d[n];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainEnergyDensity::value(FEElement& el)
{
    FEElasticMaterial* pme = m_pfem->GetMaterial(el.GetMatID())->GetElasticMaterial();
    if ((pme == 0) || pme->IsRigid()) return 0;
    
    double sed;
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
        sed = pme->StrainEnergyDensity(pt);
		val += sed;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDevStrainEnergyDensity::value(FEElement& el)
{
    FEElasticMaterial* pme = m_pfem->GetMaterial(el.GetMatID())->GetElasticMaterial();
    FEUncoupledMaterial* pmu = dynamic_cast<FEUncoupledMaterial*>(pme);
    if ((pme == 0) || pme->IsRigid() || (pmu == 0)) return 0;
    
    double sed;
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
        sed = pmu->DevStrainEnergyDensity(pt);
		val += sed;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFiberStretch::value(FEElement& el)
{
	int n = el.GaussPoints();
	double l = 0.0;
	for (int j=0; j<n; ++j)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
		vec3d ri = pt.m_Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.norm();
	}
	l /= (double) n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogElemFiberVectorX::value(FEElement& el)
{
	int n = el.GaussPoints();
	double l = 0.0;
	for (int j = 0; j<n; ++j)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
		vec3d ri = pt.m_Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.x;
	}
	l /= (double)n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogElemFiberVectorY::value(FEElement& el)
{
	int n = el.GaussPoints();
	double l = 0.0;
	for (int j = 0; j<n; ++j)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
		vec3d ri = pt.m_Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.y;
	}
	l /= (double)n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogElemFiberVectorZ::value(FEElement& el)
{
	int n = el.GaussPoints();
	double l = 0.0;
	for (int j = 0; j<n; ++j)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
		vec3d ri = pt.m_Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.z;
	}
	l /= (double)n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogDamage::value(FEElement& el)
{
    int nint = el.GaussPoints();
    double D = 0;
    for (int j=0; j<nint; ++j)
    {
        FEMaterialPoint& pt = *el.GetMaterialPoint(j);
        FEDamageMaterialPoint* ppd = pt.ExtractData<FEDamageMaterialPoint>();
        FEElasticMixtureMaterialPoint* pem = pt.ExtractData<FEElasticMixtureMaterialPoint>();
        FEMultigenerationMaterialPoint* pmg = pt.ExtractData<FEMultigenerationMaterialPoint>();
        if (ppd) D += (float) ppd->m_D;
        else if (pem) {
            for (int k=0; k<pem->Components(); ++k)
            {
                FEDamageMaterialPoint* ppd = pem->GetPointData(k)->ExtractData<FEDamageMaterialPoint>();
                if (ppd) D += (float) ppd->m_D;
            }
        }
        else if (pmg) {
            for (int k=0; k<pmg->Components(); ++k)
            {
                FEDamageMaterialPoint* ppd = pt.GetPointData(k)->ExtractData<FEDamageMaterialPoint>();
                FEElasticMixtureMaterialPoint* pem = pt.GetPointData(k)->ExtractData<FEElasticMixtureMaterialPoint>();
                if (ppd) D += (float) ppd->m_D;
                else if (pem)
                {
                    for (int l=0; l<pem->Components(); ++l)
                    {
                        FEDamageMaterialPoint* ppd = pem->GetPointData(l)->ExtractData<FEDamageMaterialPoint>();
                        if (ppd) D += (float) ppd->m_D;
                    }
                }
            }
        }
    }
    D /= (double) nint;
    return D;
}

//-----------------------------------------------------------------------------
double FELogRigidBodyR11::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(0,0)); }
double FELogRigidBodyR12::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(0, 1)); }
double FELogRigidBodyR13::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(0, 2)); }
double FELogRigidBodyR21::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(1, 0)); }
double FELogRigidBodyR22::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(1, 1)); }
double FELogRigidBodyR23::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(1, 2)); }
double FELogRigidBodyR31::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(2, 0)); }
double FELogRigidBodyR32::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(2, 1)); }
double FELogRigidBodyR33::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(2, 2)); }

//-----------------------------------------------------------------------------
double FELogRigidBodyPosX::value(FERigidBody& rb) { return rb.m_rt.x; }
double FELogRigidBodyPosY::value(FERigidBody& rb) { return rb.m_rt.y; }
double FELogRigidBodyPosZ::value(FERigidBody& rb) { return rb.m_rt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyVelX::value(FERigidBody& rb) { return rb.m_vt.x; }
double FELogRigidBodyVelY::value(FERigidBody& rb) { return rb.m_vt.y; }
double FELogRigidBodyVelZ::value(FERigidBody& rb) { return rb.m_vt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAccX::value(FERigidBody& rb) { return rb.m_at.x; }
double FELogRigidBodyAccY::value(FERigidBody& rb) { return rb.m_at.y; }
double FELogRigidBodyAccZ::value(FERigidBody& rb) { return rb.m_at.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngPosX::value(FERigidBody& rb) { return ((rb.GetRotation().GetVector()).x*rb.GetRotation().GetAngle()); }
double FELogRigidBodyAngPosY::value(FERigidBody& rb) { return ((rb.GetRotation().GetVector()).y*rb.GetRotation().GetAngle()); }
double FELogRigidBodyAngPosZ::value(FERigidBody& rb) { return ((rb.GetRotation().GetVector()).z*rb.GetRotation().GetAngle()); }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngVelX::value(FERigidBody& rb) { return rb.m_wt.x; }
double FELogRigidBodyAngVelY::value(FERigidBody& rb) { return rb.m_wt.y; }
double FELogRigidBodyAngVelZ::value(FERigidBody& rb) { return rb.m_wt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngAccX::value(FERigidBody& rb) { return rb.m_alt.x; }
double FELogRigidBodyAngAccY::value(FERigidBody& rb) { return rb.m_alt.y; }
double FELogRigidBodyAngAccZ::value(FERigidBody& rb) { return rb.m_alt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyQuatX::value(FERigidBody& rb) { return rb.GetRotation().x; }
double FELogRigidBodyQuatY::value(FERigidBody& rb) { return rb.GetRotation().y; }
double FELogRigidBodyQuatZ::value(FERigidBody& rb) { return rb.GetRotation().z; }
double FELogRigidBodyQuatW::value(FERigidBody& rb) { return rb.GetRotation().w; }

//-----------------------------------------------------------------------------
double FELogRigidBodyForceX::value(FERigidBody& rb) { return rb.m_Fr.x; }
double FELogRigidBodyForceY::value(FERigidBody& rb) { return rb.m_Fr.y; }
double FELogRigidBodyForceZ::value(FERigidBody& rb) { return rb.m_Fr.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyTorqueX::value(FERigidBody& rb) { return rb.m_Mr.x; }
double FELogRigidBodyTorqueY::value(FERigidBody& rb) { return rb.m_Mr.y; }
double FELogRigidBodyTorqueZ::value(FERigidBody& rb) { return rb.m_Mr.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyKineticEnergy::value(FERigidBody& rb) {
    FERigidBody&rbl = static_cast<FERigidBody&>(rb);
    return (rbl.m_mass*(rbl.m_vt*rbl.m_vt) + rbl.m_wt*(rbl.m_moi*rbl.m_wt))/2;
}

//-----------------------------------------------------------------------------
double FELogRigidConnectorForceX::value(FENLConstraint& rc) 
{ 
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_F.x : 0); 
}

double FELogRigidConnectorForceY::value(FENLConstraint& rc)
{
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_F.y : 0);
}

double FELogRigidConnectorForceZ::value(FENLConstraint& rc)
{
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_F.z : 0);
}

//-----------------------------------------------------------------------------
double FELogRigidConnectorMomentX::value(FENLConstraint& rc)
{ 
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_M.x : 0);
}

double FELogRigidConnectorMomentY::value(FENLConstraint& rc)
{ 
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_M.y : 0);
}

double FELogRigidConnectorMomentZ::value(FENLConstraint& rc)
{
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_M.z : 0);
}

//-----------------------------------------------------------------------------
double FELogVolumeConstraint::value(FENLConstraint& rc)
{
    FEVolumeConstraint* prc = dynamic_cast<FEVolumeConstraint*>(&rc);
    return (prc ? prc->m_s.m_Vt : 0);
}

//-----------------------------------------------------------------------------
double FELogVolumePressure::value(FENLConstraint& rc)
{
    FEVolumeConstraint* prc = dynamic_cast<FEVolumeConstraint*>(&rc);
    return (prc ? prc->m_s.m_p : 0);
}
