#include "stdafx.h"
#include "FEMicroMaterial2O.h"
#include "FECore/log.h"
#include "FESolidSolver2.h"
#include "FEElasticSolidDomain.h"
#include "FECore/FEAnalysis.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FECore/BC.h"
#include "FECore/tens3d.h"
#include "FEPeriodicBoundary2O.h"

//-----------------------------------------------------------------------------
FEMicroMaterialPoint2O::FEMicroMaterialPoint2O(FEMaterialPoint* mp) : FEMaterialPoint(mp)
{
	m_G.zero();
	m_Qa.zero();
	m_Pa.zero();

	m_Ca.zero();
	m_La.zero();
	m_Ha.zero();
	m_Ja.zero();

/*	m_tau.zero();

	m_inf_str_grad.zero();
	m_E.zero();
	m_H.zero();
	m_e.zero();
	m_h.zero();
	m_macro_energy = 0.;
	m_micro_energy = 0.;
	m_energy_diff = 0.;
	
	m_macro_energy_inc = 0.;
	m_micro_energy_inc = 0.;
	
	m_Ca.zero();
	m_Da.zero();
	m_Ea.zero();
	
	m_G_prev.zero();
*/
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FEMicroMaterialPoint2O::Init(bool bflag)
{
	FEMaterialPoint::Init(bflag);
}

//-----------------------------------------------------------------------------
//! create a shallow copy
FEMaterialPoint* FEMicroMaterialPoint2O::Copy()
{
	FEMicroMaterialPoint2O* pt = new FEMicroMaterialPoint2O(m_pNext?m_pNext->Copy():0);
	pt->m_G = m_G;
	pt->m_Qa = m_Qa;
	pt->m_Pa = m_Pa;
	pt->m_Ca = m_Ca;
	pt->m_La = m_La;
	pt->m_Ha = m_Ha;
	pt->m_Ja = m_Ja;
	return pt;
}

//-----------------------------------------------------------------------------
//! serialize material point data
void FEMicroMaterialPoint2O::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_G  << m_Qa << m_Pa;
		ar << m_Ca << m_La << m_Ha << m_Ja;
	}
	else
	{
		ar >> m_G  >> m_Qa >> m_Pa;
		ar >> m_Ca >> m_La >> m_Ha >> m_Ja;
	}
}

//=============================================================================

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEMicroMaterial2O, FEElasticMaterial)
	ADD_PARAMETER(m_szrve    , FE_PARAM_STRING, "RVE"     );
	ADD_PARAMETER(m_szbc     , FE_PARAM_STRING, "bc_set"  );
	ADD_PARAMETER(m_bperiodic, FE_PARAM_BOOL  , "periodic");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEMicroMaterial2O::FEMicroMaterial2O(FEModel* pfem) : FEElasticMaterial2O(pfem)
{
	// initialize parameters
	m_szrve[0] = 0;
	m_szbc[0] = 0;
	m_bperiodic = false;

	AddProperty(&m_probe, "probe", false);
}

//-----------------------------------------------------------------------------
FEMicroMaterial2O::~FEMicroMaterial2O(void)
{
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEMicroMaterial2O::CreateMaterialPointData()
{
	return new FEMicroMaterialPoint2O(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
bool FEMicroMaterial2O::Init()
{
	// initialize base class first
	if (FEElasticMaterial::Init() == false) return false;

	// load the master RVE model
	FEBioImport fim;
	if (fim.Load(m_mrve, m_szrve) == false)
	{
		return MaterialError("An error occured trying to read the RVE model from file %s.", m_szrve);
	}

	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::NEVER);

	// initialize master RVE
	if (m_mrve.InitRVE(m_bperiodic, m_szbc) == false) return MaterialError("An error occurred preparing RVE model");

	// reset the logfile mode
	felog.SetMode(nmode);

	return true;
}

//-----------------------------------------------------------------------------
// The stiffness is evaluated at the same time the stress is evaluated so we 
// can just return it here. Note that this assumes that the stress function 
// is always called prior to the tangent function.
// LTE - Note that this function is not used in the second-order implemenetation
tens4ds FEMicroMaterial2O::Tangent(FEMaterialPoint &mp)
{
	assert(false);
	tens4ds c; c.zero();
	return c;
}

//-----------------------------------------------------------------------------
// LTE - Note that this function is not used in the second-order implemenetation
mat3ds FEMicroMaterial2O::Stress(FEMaterialPoint &mp)
{
	assert(false);
	mat3ds sa; sa.zero();
	return sa;
}

//-----------------------------------------------------------------------------
// The stiffness is evaluated at the same time the stress is evaluated so we 
// can just return it here. Note that this assumes that the stress function 
// is always called prior to the tangent function.
void  FEMicroMaterial2O::Tangent2O(FEMaterialPoint& mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J)
{
	FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();
	C = mmpt.m_Ca;
	L = mmpt.m_La;
	H = mmpt.m_Ha;
	J = mmpt.m_Ja;
}

//-----------------------------------------------------------------------------
void FEMicroMaterial2O::Stress2O(FEMaterialPoint &mp)
{
	// get the deformation gradient and its gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
	mat3d F = pt.m_F;
	tens3drs G = mmpt2O.m_G;

	// solve the RVE
	bool bret = mmpt2O.m_rve.Solve(F, G);

	// make sure it converged
	if (bret == false) throw FEMultiScaleException();

	// calculate the averaged Cauchy stress
	mmpt2O.m_rve.AveragedStress2O(mmpt2O.m_Pa, mmpt2O.m_Qa);
	
	// calculate the averaged PK1 stress
//	AveragedStress2OPK1(mmpt2O.m_rve, mp, mmpt2O.m_PK1, mmpt2O.m_QK1);
	
	// calculate the averaged PK2 stress
//	mmpt2O.m_rve.AveragedStress2OPK2(mp, mmpt2O.m_S, mmpt2O.m_T);

	// calculate the averaged stiffness
//	mmpt2O.m_rve.AveragedStiffness(mp, mmpt2O.m_Ca, mmpt2O.m_Da, mmpt2O.m_Ea);

	// calculate the difference between the macro and micro energy for Hill-Mandel condition
//	calc_energy_diff(mmpt2O.m_rve, mp);	
}

//-----------------------------------------------------------------------------
//! Calculate the energy difference between the RVE problem and the macro material point
void FEMicroMaterial2O::calc_energy_diff(FEModel& rve, FEMaterialPoint& mp)
{
	// get the deformation gradient and deformation hessian
/*	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
	mat3d F = pt.m_F;
	tens3drs G = mmpt2O.m_G;

	mat3d Ftrans = F.transpose();
	tens3dls Gtrans = G.transpose();

	mat3d Finv = F.inverse();
	mat3d Finvtrans = Finv.transpose();
	tens3drs Ginv = G; Ginv.contractleg2(Finv,1); Ginv.contractleg2(Finv,2); Ginv.contractleg2(Finv,3);
	tens3dls Ginvtrans = Ginv.transpose();
	
	// calculate infinitesimal strain
	tens3d inf_strain_grad_nosym;

	inf_strain_grad_nosym.d[0] =  G.d[0];
	inf_strain_grad_nosym.d[1] =  G.d[1];
	inf_strain_grad_nosym.d[2] =  G.d[2];
	inf_strain_grad_nosym.d[3] =  0.5*(G.d[1] + G.d[6]);
	inf_strain_grad_nosym.d[4] =  0.5*(G.d[3] + G.d[7]);
	inf_strain_grad_nosym.d[5] =  0.5*(G.d[4] + G.d[8]);
	inf_strain_grad_nosym.d[6] =  0.5*(G.d[2] + G.d[12]);
	inf_strain_grad_nosym.d[7] =  0.5*(G.d[4] + G.d[13]);
	inf_strain_grad_nosym.d[8] =  0.5*(G.d[5] + G.d[14]);
	
	inf_strain_grad_nosym.d[9] =  0.5*(G.d[6] + G.d[1]);
	inf_strain_grad_nosym.d[10] = 0.5*(G.d[7] + G.d[3]);
	inf_strain_grad_nosym.d[11] = 0.5*(G.d[8] + G.d[4]);
	inf_strain_grad_nosym.d[12] = G.d[7];
	inf_strain_grad_nosym.d[13] = G.d[9];
	inf_strain_grad_nosym.d[14] = G.d[10];
	inf_strain_grad_nosym.d[15] = 0.5*(G.d[8] + G.d[13]);
	inf_strain_grad_nosym.d[16] = 0.5*(G.d[10] + G.d[15]);
	inf_strain_grad_nosym.d[17] = 0.5*(G.d[11] + G.d[16]);
	
	inf_strain_grad_nosym.d[18] = 0.5*(G.d[12] + G.d[2]);
	inf_strain_grad_nosym.d[19] = 0.5*(G.d[13] + G.d[4]);
	inf_strain_grad_nosym.d[20] = 0.5*(G.d[14] + G.d[5]);
	inf_strain_grad_nosym.d[21] = 0.5*(G.d[13] + G.d[8]);
	inf_strain_grad_nosym.d[22] = 0.5*(G.d[15] + G.d[10]);
	inf_strain_grad_nosym.d[23] = 0.5*(G.d[16] + G.d[11]);
	inf_strain_grad_nosym.d[24] = G.d[14];
	inf_strain_grad_nosym.d[25] = G.d[16];
	inf_strain_grad_nosym.d[26] = G.d[17];

	mmpt2O.m_inf_str_grad = inf_strain_grad_nosym.symm();

	// calculate Green-Lagrange strain
	mmpt2O.m_E = ((Ftrans*F - mat3dd(1))*0.5).sym();
	mmpt2O.m_H = (((Gtrans*F) + (Ftrans*G))*0.5).symm();

	// calculate Euler-Almansi strain
	mmpt2O.m_e = ((mat3dd(1) - Finvtrans*Finv)*0.5).sym();
	mmpt2O.m_h = (((Ginvtrans*Finv) + (Finvtrans*Ginv))*-0.5).symm();
	
	// calculate the energy difference between macro point and RVE
	// to verify that we have satisfied the Hill-Mandel condition
	mat3d F_prev = mmpt2O.m_F_prev;
	tens3drs G_prev = mmpt2O.m_G_prev;

	// calculate the macroscopic strain energy increment according to PK1 stress
	mmpt2O.m_macro_energy_inc = mmpt2O.m_PK1.dotdot(pt.m_F - F_prev) + mmpt2O.m_QK1.tripledot(mmpt2O.m_G - mmpt2O.m_G_prev);

	//// calculate the macroscopic strain energy increment according to PK2 stress
	///*mat3ds E_prev = ((F_prev.transpose()*F_prev - mat3dd(1))*0.5).sym();
	//tens3ds H_prev = ((G_prev.transpose().multiply2right(F).LStoUnsym() + G_prev.multiply2left(Ftrans).RStoUnsym())*0.5).symm();
	//mmpt2O.m_macro_energy_inc = mmpt2O.m_S.dotdot(mmpt2O.m_E - E_prev) + mmpt2O.m_T.tripledot(mmpt2O.m_H - H_prev);
	//
	//// calculate the macroscopic strain energy increment according to Cauchy stress
	///*mat3d Finv_prev = F_prev.inverse();
	//mat3d Finvtrans_prev = Finv_prev.transpose();
	//mat3ds e_prev = ((mat3dd(1) - F_prev.transinv()*F_prev.inverse())*0.5).sym();
	//tens3drs Ginv_prev; Ginv_prev = G_prev; Ginv_prev.contractleg2(Finv_prev,1); Ginv_prev.contractleg2(Finv_prev,2); Ginv_prev.contractleg2(Finv_prev,3);
	//tens3dls Ginvtrans_prev = Ginv_prev.transpose();
	//tens3ds h_prev = ((Ginvtrans_prev.multiply2right(Finv_prev).LStoUnsym() + Ginv_prev.multiply2left(Finvtrans_prev).RStoUnsym())*-0.5).symm();
	//mmpt2O.m_macro_energy_inc = pt.m_s.dotdot(mmpt2O.m_e - e_prev) + mmpt2O.m_tau.tripledot(mmpt2O.m_h - h_prev);

	// calculate the microscopic strain energy increment
	double rve_energy_avg = 0.;
	int nint; 
	double* w;
	
	mat3d rve_F;		
	mat3d rve_F_prev;
	double J0 = 0.;
	double V0 = 0.;
	mat3d rve_PK1;
	
	mat3ds rve_S;
	mat3ds rve_E;
	mat3ds rve_E_prev;
		
	mat3ds rve_s;
	mat3ds rve_e;
	mat3ds rve_e_prev;
	double J = 0.;
	double v = 0.;

	FEMesh& m = mmpt2O.m_rve.GetMesh();

	for (int k=0; k<m.Domains(); ++k)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(m.Domain(k));

		for (int i=0; i<dom.Elements(); ++i)
		{
			FESolidElement& el = dom.Element(i);
			
			nint = el.GaussPoints();
			w = el.GaussWeights();
			
			for (int n=0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();

				rve_F = pt.m_F;
				rve_F_prev = mmpt.m_F_prev;
				rve_s = pt.m_s;

				//// calculate microscopic strain energy according to PK1 stress
				//rve_PK1 = rve_F.det()*rve_s*rve_F.transinv();
				//J0 = dom.detJ0(el, n);		
				//V0 += J0*w[n];
				//rve_energy_avg += rve_PK1.dotdot(rve_F - rve_F_prev)*J0*w[n];
				//				
				//// calculate microscopic strain energy according to PK2 stress
				///*rve_E = ((rve_F.transpose()*rve_F - mat3dd(1))*0.5).sym();
				//rve_E_prev = ((rve_F_prev.transpose()*rve_F_prev - mat3dd(1))*0.5).sym();
				//rve_PK1 = rve_F.det()*rve_s*rve_F.transinv();
				//rve_S = (rve_F.inverse()*rve_PK1).sym();
				//J0 = dom.detJ0(el, n);		
				//V0 += J0*w[n];
				//rve_energy_avg += rve_S.dotdot(rve_E - rve_E_prev)*J0*w[n];
				//
				//// calculate microscopic strain energy according to Cauchy stress
				///*rve_s = rve_pt.m_s;
				//rve_e = ((mat3dd(1) - rve_F.transinv()*rve_F.inverse())*0.5).sym();
				//rve_e_prev = ((mat3dd(1) - rve_F_prev.transinv()*rve_F_prev.inverse())*0.5).sym();
				//J = dom.detJt(el, n);		
				//v += J*w[n];		
				//rve_energy_avg += rve_s.dotdot(rve_e - rve_e_prev)*J*w[n];
			}
		}
	}

	mmpt2O.m_micro_energy_inc = rve_energy_avg/V0;
*/
}
