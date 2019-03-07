#include "stdafx.h"
#include "FEBioModel.h"
#include <FECore/FEAnalysis.h>
#include <FECore/tens3d.h>
#include <FEBioMech/FERigidJoint.h>
#include <FEBioMech/FERigidSphericalJoint.h>
#include <FEBioPlot/FEBioPlotFile.h>
#include <FEBioMech/FERigidSystem.h>
#include <FECore/FEMaterial.h>
#include <FEBioMech/FERigidBody.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FELoadCurve.h>
#include <FECore/log.h>
#include <string.h>

//-----------------------------------------------------------------------------
// helper function to print a parameter to the logfile
void FEBioModel::print_parameter(FEParam& p, int level)
{
	char sz[512] = {0};
	int l = (int)strlen(p.name()) + 2*level;
	sprintf(sz, "\t%*s %.*s", l, p.name(), 50-l, "..................................................");
	if (p.dim() == 1)
	{
		switch (p.type())
		{
		case FE_PARAM_DOUBLE : feLog("%s : %lg\n", sz, p.value<double>()); break;
		case FE_PARAM_INT    : feLog("%s : %d\n" , sz, p.value<int   >()); break;
		case FE_PARAM_BOOL   : feLog("%s : %s\n" , sz, (p.value<bool>() ? "true" : "false")); break;
		case FE_PARAM_STRING : feLog("%s : %s\n" , sz, p.cvalue()); break;
		case FE_PARAM_STD_STRING: feLog("%s : %s\n", sz, p.value<string>().c_str()); break;
		case FE_PARAM_VEC3D  :
			{
				vec3d v = p.value<vec3d>();
				feLog("%s : %lg,%lg,%lg\n", sz, v.x, v.y, v.z);
			}
			break;
		case FE_PARAM_MAT3DS :
			{
				mat3ds m = p.value<mat3ds>();
				feLog("%s : %lg,%lg,%lg,%lg,%lg,%lg\n", sz, m.xx(), m.yy(), m.zz(), m.xy(), m.yz(), m.xz());
			}
			break;
		case FE_PARAM_MAT3D:
			{
				mat3d m = p.value<mat3d>();
				feLog("%s : %lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n", sz, m(0,0), m(0,1), m(0,2), m(1,0), m(1,1), m(1,2), m(2,0), m(2,1), m(2,2));
			}
			break;
		case FE_PARAM_TENS3DRS:
			{
				tens3drs m = p.value<tens3drs>();
				double* d = m.d;
				feLog("%s : %lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n", sz,d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14],d[15],d[16],d[17]);
			}
			break;
		default:
			break;
		}
	}
	else
	{
		switch (p.type())
		{
		case FE_PARAM_INT:
		case FE_PARAM_DOUBLE:
			{
				int n = p.dim();
				feLog("%s : ", sz);
				for (int k=0; k<n; ++k)
				{
					switch (p.type())
					{
					case FE_PARAM_INT   : feLog("%d" , p.pvalue<int   >()[k]); break;
					case FE_PARAM_DOUBLE: feLog("%lg", p.pvalue<double>()[k]); break;
                    default: break;
					}
					if (k!=n-1) feLog(","); else feLog("\n");
				}
			}
			break;
		case FE_PARAM_DOUBLE_MAPPED:
		{
			int n = p.dim();
			feLog("%s : ", sz);
			for (int k = 0; k < n; ++k)
			{
				FEParamDouble& v = p.value<FEParamDouble>(k);
				if (v.isConst())
					feLog("%lg", v.constValue());
				else
					feLog("???");
				if (k != n - 1) feLog(","); else feLog("\n");
			}
		}
		break;
		default:
			break;
		}
	}
}

//-----------------------------------------------------------------------------
// print the parameter list to the log file
void FEBioModel::print_parameter_list(FEParameterList& pl, int level)
{
	int n = pl.Parameters();
	if (n > 0)
	{
		FEParamIterator it = pl.first();
		for (int j=0; j<n; ++j, ++it) print_parameter(*it, level);
	}
}

//------------------------------------------------------------------------------
void FEBioModel::print_parameter_list(FECoreBase* pc, int level)
{
	print_parameter_list(pc->GetParameterList(), level);
	for (int i=0; i<pc->PropertyClasses(); ++i)
	{
		FEProperty* prop = pc->PropertyClass(i);
		int n = prop->size();
		for (int j=0; j<n; ++j)
		{
			feLog("\t%s: ", prop->GetName());
			FECoreBase* pcj = prop->get(j);
			if (pcj)
			{
				feLog("(type: %s)\n", pcj->GetTypeStr());
				print_parameter_list(pcj, level + 1);
			}
			else feLog("(unspecified)\n");
		}
	}
}

//------------------------------------------------------------------------------
//! This function outputs the input data to the felog file.
void FEBioModel::echo_input()
{
	// get the analysis step
	FEAnalysis& step = *GetCurrentStep();

	// get the FE mesh
	FEBioModel& fem = *this;
	FEMesh& mesh = GetMesh();

	// print title
	feLog("%s\n\n", fem.GetTitle().c_str());

	// print file info
	feLog(" FILES USED\n");
	feLog("===========================================================================\n");
	feLog("\tInput file : %s\n", fem.GetInputFileName());
	feLog("\tPlot file  : %s\n", fem.GetPlotFileName());
	feLog("\tLog file   : %s\n", fem.GetLogfileName());
	feLog("\n\n");

	// print mesh info
	feLog(" MESH INFO\n");
	feLog("===========================================================================\n");
	feLog("\tNumber of materials ............................ : %d\n", fem.Materials());
	feLog("\tNumber of domains .............................. : %d\n", mesh.Domains());
	feLog("\tNumber of nodes ................................ : %d\n", mesh.Nodes());
	int nsolid = mesh.Elements(FE_DOMAIN_SOLID    ); if (nsolid > 0) feLog("\tNumber of solid elements ....................... : %d\n", nsolid);
	int nshell = mesh.Elements(FE_DOMAIN_SHELL    ); if (nshell > 0) feLog("\tNumber of shell elements ....................... : %d\n", nshell);
	int ntruss = mesh.Elements(FE_DOMAIN_TRUSS    ); if (ntruss > 0) feLog("\tNumber of truss elements ....................... : %d\n", ntruss);
	int nelm2d = mesh.Elements(FE_DOMAIN_2D       ); if (nelm2d > 0) feLog("\tNumber of 2D elements .......................... : %d\n", nelm2d);
	feLog("\n\n");

	// print control info
	feLog(" CONTROL DATA\n");
	feLog("===========================================================================\n");
	const char* szmod = step.GetFESolver()->GetTypeStr();
	if (szmod == 0) { szmod = "unknown"; assert(false); }
	feLog("\tModule type .................................... : %s\n", szmod);

	const char* szan = 0;
	switch (step.m_nanalysis)
	{
	case FE_STATIC      : szan = "quasi-static"; break;
	case FE_DYNAMIC     : szan = "dynamic"     ; break;
	case FE_STEADY_STATE: szan = "steady-state"; break;
	default:
		szan = "unknown";
		assert(false);
	}
	feLog("\tAnalysis type .................................. : %s\n", szan);

	if (step.m_ntime > 0)
		feLog("\tNumber of timesteps ............................ : %d\n", step.m_ntime);
	else
		feLog("\tFinal time ..................................... : %lg\n", step.m_final_time);

	feLog("\tTime step size ................................. : %lg\n", step.m_dt0);
	feLog("\tAuto time stepper activated .................... : %s\n", (step.m_bautostep ? "yes" : "no"));
	if (step.m_bautostep)
	{
		FETimeStepController& tc = step.m_timeController;
		feLog("\t  Optimal nr of iterations ..................... : %d\n", tc.m_iteopt);
		feLog("\t  Minimum allowable step size .................. : %lg\n", tc.m_dtmin);
		feLog("\t  Maximum allowable step size .................. : %lg\n", tc.m_dtmax);
	}
	feLog("\tNumber of load controllers ..................... : %d\n", fem.LoadControllers());

	feLog("\n\n");

	// output solver data
	feLog(" SOLVER PARAMETERS\n");
	feLog("===========================================================================\n");

	FESolver* psolver = step.GetFESolver();
	if (psolver)
	{
		print_parameter_list(psolver->GetParameterList());
		feLog("\n\n");
	}

	// print output data
	feLog(" OUTPUT DATA\n");
	feLog("===========================================================================\n");
	switch (step.m_nplot)
	{
	case FE_PLOT_NEVER      : feLog("\tplot level ................................ : never\n"); break;
	case FE_PLOT_MAJOR_ITRS : feLog("\tplot level ................................ : major iterations\n"); break;
	case FE_PLOT_MINOR_ITRS : feLog("\tplot level ................................ : minor iterations\n"); break;
	case FE_PLOT_MUST_POINTS: feLog("\tplot level ................................ : must points only\n"); break;
	case FE_PLOT_FINAL      : feLog("\tplot level ................................ : final state\n"); break;
	case FE_PLOT_STEP_FINAL : feLog("\tplot level ................................ : step final state\n"); break;
	}

	PlotFile* pplt = fem.GetPlotFile();
	if (dynamic_cast<FEBioPlotFile*>(pplt))
	{
		FEBioPlotFile* pf = dynamic_cast<FEBioPlotFile*>(pplt);
		feLog("\tplotfile format ........................... : FEBIO\n");

		const FEBioPlotFile::Dictionary& dic = pf->GetDictionary();

		for (int i=0; i<5; ++i)
		{
			const list<FEBioPlotFile::DICTIONARY_ITEM>* pl=0;
			const char* szn = 0;
			switch (i)
			{
			case 0: pl = &dic.GlobalVariableList  (); szn = "Global Variables"  ; break;
			case 1: pl = &dic.MaterialVariableList(); szn = "Material Variables"; break;
			case 2: pl = &dic.NodalVariableList   (); szn = "Nodal Variables"   ; break;
			case 3: pl = &dic.DomainVariableList  (); szn = "Domain Variables"  ; break;
			case 4: pl = &dic.SurfaceVariableList (); szn = "Surface Variables" ; break;
			}

			if (!pl->empty())
			{
				feLog("\t\t%s:\n", szn);
				list<FEBioPlotFile::DICTIONARY_ITEM>::const_iterator it;
				for (it = pl->begin(); it != pl->end(); ++it)
				{
					const char* szt = 0;
					switch (it->m_ntype)
					{
					case PLT_FLOAT  : szt = "float"  ; break;
					case PLT_VEC3F  : szt = "vec3f"  ; break;
					case PLT_MAT3FS : szt = "mat3fs" ; break;
					case PLT_MAT3FD : szt = "mat3fd" ; break;
                    case PLT_TENS4FS: szt = "tens4fs"; break;
					case PLT_MAT3F  : szt = "mat3f"  ; break;
					}

					const char* szf = 0;
					switch (it->m_nfmt)
					{
					case FMT_NODE: szf = "NODE"; break;
					case FMT_ITEM: szf = "ITEM"; break;
					case FMT_MULT: szf = "COMP"; break;
					case FMT_REGION: szf = "REGION"; break;
					}

					feLog("\t\t\t%-20s (type = %5s, format = %4s)\n", it->m_szname, szt, szf);
				}
			}
		}
	}

	// material data
	feLog("\n\n");
	feLog(" MATERIAL DATA\n");
	feLog("===========================================================================\n");
	for (int i=0; i<fem.Materials(); ++i)
	{
		if (i>0) feLog("---------------------------------------------------------------------------\n");
		feLog("%3d - ", i+1);

		// get the material
		FEMaterial* pmat = fem.GetMaterial(i);

		// get the material name and type string
		const char* szname = pmat->GetName().c_str();
		const char* sztype = pmat->GetTypeStr();
		if (szname[0] == 0) szname = 0;

		// print type and name
		feLog("%s", (szname?szname:"unknown"));
		feLog(" (type: %s)", sztype);
		feLog("\n");

		// print the parameter list
		print_parameter_list(pmat);
	}
	feLog("\n\n");

	FERigidSystem& rigid = *fem.GetRigidSystem();
	if (rigid.Objects())
	{
		feLog(" RIGID BODY DATA\n");
		feLog("===========================================================================\n");
		for (int i=0; i<rigid.Objects(); ++i)
		{
			FERigidBody& rb = *rigid.Object(i);
			if (i>0) feLog("---------------------------------------------------------------------------\n");
			feLog("Rigid Body %d:\n", rb.m_nID+1);
			feLog("\tmaterial id    : %d\n", rb.m_mat+1);
			feLog("\tcenter of mass : %lg, %lg, %lg\n", rb.m_r0.x, rb.m_r0.y, rb.m_r0.z);
            feLog("\tmass           : %lg\n", rb.m_mass);
            feLog("\tIxx Ixy Ixz    : %lg, %lg, %lg\n", rb.m_moi.xx(), rb.m_moi.xy(), rb.m_moi.xz());
            feLog("\tIxy Iyy Iyz    : %lg, %lg, %lg\n", rb.m_moi.xy(), rb.m_moi.yy(), rb.m_moi.yz());
            feLog("\tIxz Iyz Izz    : %lg, %lg, %lg\n", rb.m_moi.xz(), rb.m_moi.yz(), rb.m_moi.zz());
		}
		feLog("\n\n");
	}

	if (fem.BodyLoads() > 0)
	{
		feLog(" BODY LOAD DATA\n");
		feLog("===========================================================================\n");
		for (int i=0; i<fem.BodyLoads(); ++i)
		{
			if (i>0) feLog("---------------------------------------------------------------------------\n");
			feLog("%3d - ", i+1);

			// get the body load
			FEBodyLoad* pbl = fem.GetBodyLoad(i);

			// get the type string
			const char* sztype = pbl->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog(" Type: %s\n", sztype);

			// print the parameter list
			FEParameterList& pl = pbl->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	if (fem.SurfacePairConstraints() > 0)
	{
		feLog(" CONTACT INTERFACE DATA\n");
		feLog("===========================================================================\n");
		for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
		{
			if (i>0) feLog("---------------------------------------------------------------------------\n");

			FESurfacePairConstraint* pi = fem.SurfacePairConstraint(i);
			const char* sztype = pi->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog("contact interface %d - Type: %s\n", i+1, sztype);
			FEParameterList& pl = pi->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	if (fem.NonlinearConstraints() != 0)
	{
		feLog(" NONLINEAR CONSTRAINT DATA\n");
		feLog("===========================================================================\n");
		int NC = fem.NonlinearConstraints();
		for (int i = 0; i<NC; ++i)
		{
			FENLConstraint* plc = fem.NonlinearConstraint(i);
			if (dynamic_cast<FERigidJoint*>(plc))
			{
				FERigidJoint& rj = static_cast<FERigidJoint&>(*plc);
				FERigidBody& ra = *rigid.Object(rj.m_nRBa);
				FERigidBody& rb = *rigid.Object(rj.m_nRBb);
				feLog("rigid joint %d:\n", i+1);
				feLog("\tRigid body A                   : %d\n", ra.m_mat + 1);
				feLog("\tRigid body B                   : %d\n", rb.m_mat + 1);
				feLog("\tJoint                          : (%lg, %lg, %lg)\n", rj.m_q0.x, rj.m_q0.y, rj.m_q0.z);
				feLog("\tPenalty factor                 : %lg\n", rj.m_eps );
				feLog("\tAugmented Lagrangian tolerance : %lg\n", rj.m_atol);
				feLog("---------------------------------------------------------------------------\n");
			}
		}
		feLog("\n\n");
	}

	FECoreKernel& fecore = FECoreKernel::GetInstance();
	feLog(" LINEAR SOLVER DATA\n");
	feLog("===========================================================================\n");
	feLog("\tSolver type ....................................... : %s\n", fecore.GetLinearSolverType());
	feLog("\tMatrix format ..................................... : ");
	switch (step.GetFESolver()->MatrixSymmetryFlag())
	{
	case REAL_UNSYMMETRIC   : feLog("unsymmetric\n"); break;
	case REAL_SYMMETRIC     : feLog("symmetric\n"); break;
	case REAL_SYMM_STRUCTURE: feLog("symmetric structure\n"); break;
	}
	FECoreFactory* fac = fecore.FindFactoryClass(FELINEARSOLVER_ID, fecore.GetLinearSolverType());
	if (fac) print_parameter_list(fac->GetParameterList());
	feLog("\n");
}
