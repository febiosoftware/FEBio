#include "stdafx.h"
#include "FEBioMechData.h"
#include "FEElasticMaterial.h"
#include "FERigidMaterial.h"
#include "FESolidSolver.h"
#include "FECore/FERigidBody.h"

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
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.x - node.m_r0.x; 
}

//-----------------------------------------------------------------------------
double FENodeYDisp::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.y - node.m_r0.y; 
}

//-----------------------------------------------------------------------------
double FENodeZDisp::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.z - node.m_r0.z; 
}

//-----------------------------------------------------------------------------
double FENodeXVel::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_vt.x; 
}

//-----------------------------------------------------------------------------
double FENodeYVel::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_vt.y; 
}

//-----------------------------------------------------------------------------
double FENodeZVel::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_vt.z; 
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
	FESolidSolver* psolid_solver = dynamic_cast<FESolidSolver*>(m_pfem->GetCurrentStep()->m_psolver);
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
	FESolidSolver* psolid_solver = dynamic_cast<FESolidSolver*>(m_pfem->GetCurrentStep()->m_psolver);
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
	FESolidSolver* psolid_solver = dynamic_cast<FESolidSolver*>(m_pfem->GetCurrentStep()->m_psolver);
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
double FELogRigidBodyR11::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(0,0)); }
double FELogRigidBodyR12::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(0,1)); }
double FELogRigidBodyR13::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(0,2)); }
double FELogRigidBodyR21::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(1,0)); }
double FELogRigidBodyR22::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(1,1)); }
double FELogRigidBodyR23::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(1,2)); }
double FELogRigidBodyR31::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(2,0)); }
double FELogRigidBodyR32::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(2,1)); }
double FELogRigidBodyR33::value(FEObject& rb) { FERigidBody& o = static_cast<FERigidBody&>(rb); return (o.m_qt.RotationMatrix()(2,2)); }

//-----------------------------------------------------------------------------
double FELogRigidBodyPosX::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_rt.x; }
double FELogRigidBodyPosY::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_rt.y; }
double FELogRigidBodyPosZ::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_rt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyVelX::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_vt.x; }
double FELogRigidBodyVelY::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_vt.y; }
double FELogRigidBodyVelZ::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_vt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAccX::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_at.x; }
double FELogRigidBodyAccY::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_at.y; }
double FELogRigidBodyAccZ::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_at.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngPosX::value(FEObject& rb) { FERigidBody& o = dynamic_cast<FERigidBody&>(rb); return ((o.m_qt.GetVector()).x*o.m_qt.GetAngle()); }
double FELogRigidBodyAngPosY::value(FEObject& rb) { FERigidBody& o = dynamic_cast<FERigidBody&>(rb); return ((o.m_qt.GetVector()).y*o.m_qt.GetAngle()); }
double FELogRigidBodyAngPosZ::value(FEObject& rb) { FERigidBody& o = dynamic_cast<FERigidBody&>(rb); return ((o.m_qt.GetVector()).z*o.m_qt.GetAngle()); }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngVelX::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_wt.x; }
double FELogRigidBodyAngVelY::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_wt.y; }
double FELogRigidBodyAngVelZ::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_wt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngAccX::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_alt.x; }
double FELogRigidBodyAngAccY::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_alt.y; }
double FELogRigidBodyAngAccZ::value(FEObject& rb) { return dynamic_cast<FERigidBody&>(rb).m_alt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyQuatX::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_qt.x; }
double FELogRigidBodyQuatY::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_qt.y; }
double FELogRigidBodyQuatZ::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_qt.z; }
double FELogRigidBodyQuatW::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_qt.w; }

//-----------------------------------------------------------------------------
double FELogRigidBodyForceX::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_Fr.x; }
double FELogRigidBodyForceY::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_Fr.y; }
double FELogRigidBodyForceZ::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_Fr.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyTorqueX::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_Mr.x; }
double FELogRigidBodyTorqueY::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_Mr.y; }
double FELogRigidBodyTorqueZ::value(FEObject& rb) { return static_cast<FERigidBody&>(rb).m_Mr.z; }
