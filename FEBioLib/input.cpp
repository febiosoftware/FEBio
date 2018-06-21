#include "stdafx.h"
#include "FEBioModel.h"
#include <FECore/FEAnalysis.h>
#include <FECore/tens3d.h>
#include <FECore/FEMathValue.h>
#include "FEBioMech/FERigidJoint.h"
#include "FEBioMech/FERigidSphericalJoint.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FECore/FERigidSystem.h"
#include "FECore/FEMaterial.h"
#include "FECore/FERigidBody.h"
#include "FECore/BC.h"
#include "FECore/FEBodyLoad.h"
#include <FECore/FEDataLoadCurve.h>
#include "FECore/log.h"
#include <string.h>

//-----------------------------------------------------------------------------
// helper function to print a parameter to the logfile
void print_parameter(FEParam& p)
{
	char sz[512] = {0};
	int l = (int)strlen(p.name());
	sprintf(sz, "\t%-*s %.*s", l, p.name(), 50-l, "..................................................");
	if (p.dim() == 1)
	{
		switch (p.type())
		{
		case FE_PARAM_DOUBLE : felog.printf("%s : %lg\n", sz, p.value<double>()); break;
		case FE_PARAM_INT    : felog.printf("%s : %d\n" , sz, p.value<int   >()); break;
		case FE_PARAM_BOOL   : felog.printf("%s : %s\n" , sz, (p.value<bool>() ? "true" : "false")); break;
		case FE_PARAM_STRING : felog.printf("%s : %s\n" , sz, p.cvalue()); break;
		case FE_PARAM_VEC3D  :
			{
				vec3d v = p.value<vec3d>();
				felog.printf("%s : %lg,%lg,%lg\n", sz, v.x, v.y, v.z);
			}
			break;
		case FE_PARAM_MAT3DS :
			{
				mat3ds m = p.value<mat3ds>();
				felog.printf("%s : %lg,%lg,%lg,%lg,%lg,%lg\n", sz, m.xx(), m.yy(), m.zz(), m.xy(), m.yz(), m.xz());
			}
			break;
		case FE_PARAM_MAT3D:
			{
				mat3d m = p.value<mat3d>();
				felog.printf("%s : %lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n", sz, m(0,0), m(0,1), m(0,2), m(1,0), m(1,1), m(1,2), m(2,0), m(2,1), m(2,2));
			}
			break;
		case FE_PARAM_TENS3DRS:
			{
				tens3drs m = p.value<tens3drs>();
				double* d = m.d;
				felog.printf("%s : %lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n", sz,d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14],d[15],d[16],d[17]);
			}
			break;
		case FE_PARAM_MATH_DOUBLE:
			{
				FEMathDouble& mv = p.value<FEMathDouble>();
				std::string s = mv.getExpression();
				felog.printf("%s : %s\n", sz, s.c_str());
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
				felog.printf("%s : ", sz);
				for (int k=0; k<n; ++k)
				{
					switch (p.type())
					{
					case FE_PARAM_INT   : felog.printf("%d" , p.pvalue<int   >()[k]); break;
					case FE_PARAM_DOUBLE: felog.printf("%lg", p.pvalue<double>()[k]); break;
                    default: break;
					}
					if (k!=n-1) felog.printf(","); else felog.printf("\n");
				}
			}
			break;
		default:
			assert(false);
		}
	}
}

//-----------------------------------------------------------------------------
// print the parameter list to the log file
void print_parameter_list(FEParameterList& pl)
{
	int n = pl.Parameters();
	if (n > 0)
	{
		FEParamIterator it = pl.first();
		for (int j=0; j<n; ++j, ++it) print_parameter(*it);
	}
}

//------------------------------------------------------------------------------
//! This function outputs the input data to the felog file.
void echo_input(FEBioModel& fem)
{
	// get the analysis step
	FEAnalysis& step = *fem.GetCurrentStep();

	// get the FE mesh
	FEMesh& mesh = fem.GetMesh();

	// print title
	felog.printf("%s\n\n", fem.GetTitle());

	// print file info
	felog.printf(" FILES USED\n");
	felog.printf("===========================================================================\n");
	felog.printf("\tInput file : %s\n", fem.GetInputFileName());
	felog.printf("\tPlot file  : %s\n", fem.GetPlotFileName());
	felog.printf("\tLog file   : %s\n", fem.GetLogfileName());
	felog.printf("\n\n");

	// print mesh info
	felog.printf(" MESH INFO\n");
	felog.printf("===========================================================================\n");
	felog.printf("\tNumber of materials ............................ : %d\n", fem.Materials());
	felog.printf("\tNumber of domains .............................. : %d\n", mesh.Domains());
	felog.printf("\tNumber of nodes ................................ : %d\n", mesh.Nodes());
	int nsolid = mesh.Elements(FE_DOMAIN_SOLID    ); if (nsolid > 0) felog.printf("\tNumber of solid elements ....................... : %d\n", nsolid);
	int nshell = mesh.Elements(FE_DOMAIN_SHELL    ); if (nshell > 0) felog.printf("\tNumber of shell elements ....................... : %d\n", nshell);
	int ntruss = mesh.Elements(FE_DOMAIN_TRUSS    ); if (ntruss > 0) felog.printf("\tNumber of truss elements ....................... : %d\n", ntruss);
	int nelm2d = mesh.Elements(FE_DOMAIN_2D       ); if (nelm2d > 0) felog.printf("\tNumber of 2D elements .......................... : %d\n", nelm2d);
	felog.printf("\n\n");

	// print control info
	felog.printf(" CONTROL DATA\n");
	felog.printf("===========================================================================\n");
	const char* szmod = step.GetFESolver()->GetTypeStr();
	if (szmod == 0) { szmod = "unknown"; assert(false); }
	felog.printf("\tModule type .................................... : %s\n", szmod);

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
	felog.printf("\tAnalysis type .................................. : %s\n", szan);

	if (step.m_ntime > 0)
		felog.printf("\tNumber of timesteps ............................ : %d\n", step.m_ntime);
	else
		felog.printf("\tFinal time ..................................... : %lg\n", step.m_final_time);

	felog.printf("\tTime step size ................................. : %lg\n", step.m_dt0);
	felog.printf("\tAuto time stepper activated .................... : %s\n", (step.m_bautostep ? "yes" : "no"));
	if (step.m_bautostep)
	{
		FETimeStepController& tc = step.m_timeController;
		felog.printf("\t  Optimal nr of iterations ..................... : %d\n", tc.m_iteopt);
		felog.printf("\t  Minimum allowable step size .................. : %lg\n", tc.m_dtmin);
		felog.printf("\t  Maximum allowable step size .................. : %lg\n", tc.m_dtmax);
	}
	felog.printf("\tNumber of loadcurves ........................... : %d\n", fem.LoadCurves());

	felog.printf("\n\n");

	// output solver data
	felog.printf(" SOLVER PARAMETERS\n");
	felog.printf("===========================================================================\n");

	FESolver* psolver = step.GetFESolver();
	if (psolver)
	{
		print_parameter_list(psolver->GetParameterList());
		felog.printf("\n\n");
	}

	// print output data
	felog.printf(" OUTPUT DATA\n");
	felog.printf("===========================================================================\n");
	switch (step.m_nplot)
	{
	case FE_PLOT_NEVER      : felog.printf("\tplot level ................................ : never\n"); break;
	case FE_PLOT_MAJOR_ITRS : felog.printf("\tplot level ................................ : major iterations\n"); break;
	case FE_PLOT_MINOR_ITRS : felog.printf("\tplot level ................................ : minor iterations\n"); break;
	case FE_PLOT_MUST_POINTS: felog.printf("\tplot level ................................ : must points only\n"); break;
	case FE_PLOT_FINAL      : felog.printf("\tplot level ................................ : final state\n"); break;
	case FE_PLOT_STEP_FINAL : felog.printf("\tplot level ................................ : step final state\n"); break;
	}

	PlotFile* pplt = fem.GetPlotFile();
	if (dynamic_cast<FEBioPlotFile*>(pplt))
	{
		FEBioPlotFile* pf = dynamic_cast<FEBioPlotFile*>(pplt);
		felog.printf("\tplotfile format ........................... : FEBIO\n");

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
				felog.printf("\t\t%s:\n", szn);
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

					felog.printf("\t\t\t%-20s (type = %5s, format = %4s)\n", it->m_szname, szt, szf);
				}
			}
		}
	}

	// material data
	felog.printf("\n\n");
	felog.printf(" MATERIAL DATA\n");
	felog.printf("===========================================================================\n");
	for (int i=0; i<fem.Materials(); ++i)
	{
		if (i>0) felog.printf("---------------------------------------------------------------------------\n");
		felog.printf("%3d - ", i+1);

		// get the material
		FEMaterial* pmat = fem.GetMaterial(i);

		// get the material name and type string
		const char* szname = pmat->GetName().c_str();
		const char* sztype = pmat->GetTypeStr();
		if (szname[0] == 0) szname = 0;

		// print type and name
		felog.printf("%s", (szname?szname:"unknown"));
		felog.printf(" (type: %s)", sztype);
		felog.printf("\n");

		// print the parameter list
		FEParameterList& pl = pmat->GetParameterList();
		print_parameter_list(pl);
	}
	felog.printf("\n\n");

	FERigidSystem& rigid = *fem.GetRigidSystem();
	if (rigid.Objects())
	{
		felog.printf(" RIGID BODY DATA\n");
		felog.printf("===========================================================================\n");
		for (int i=0; i<rigid.Objects(); ++i)
		{
			FERigidBody& rb = *rigid.Object(i);
			if (i>0) felog.printf("---------------------------------------------------------------------------\n");
			felog.printf("Rigid Body %d:\n", rb.m_nID+1);
			felog.printf("\tmaterial id    : %d\n", rb.m_mat+1);
			felog.printf("\tcenter of mass : %lg, %lg, %lg\n", rb.m_r0.x, rb.m_r0.y, rb.m_r0.z);
            felog.printf("\tmass           : %lg\n", rb.m_mass);
            felog.printf("\tIxx Ixy Ixz    : %lg, %lg, %lg\n", rb.m_moi.xx(), rb.m_moi.xy(), rb.m_moi.xz());
            felog.printf("\tIxy Iyy Iyz    : %lg, %lg, %lg\n", rb.m_moi.xy(), rb.m_moi.yy(), rb.m_moi.yz());
            felog.printf("\tIxz Iyz Izz    : %lg, %lg, %lg\n", rb.m_moi.xz(), rb.m_moi.yz(), rb.m_moi.zz());
		}
		felog.printf("\n\n");
	}

	if (fem.BodyLoads() > 0)
	{
		felog.printf(" BODY LOAD DATA\n");
		felog.printf("===========================================================================\n");
		for (int i=0; i<fem.BodyLoads(); ++i)
		{
			if (i>0) felog.printf("---------------------------------------------------------------------------\n");
			felog.printf("%3d - ", i+1);

			// get the body load
			FEBodyLoad* pbl = fem.GetBodyLoad(i);

			// get the type string
			const char* sztype = pbl->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			felog.printf(" Type: %s\n", sztype);

			// print the parameter list
			FEParameterList& pl = pbl->GetParameterList();
			print_parameter_list(pl);
		}
		felog.printf("\n\n");
	}

	if (fem.SurfacePairConstraints() > 0)
	{
		felog.printf(" CONTACT INTERFACE DATA\n");
		felog.printf("===========================================================================\n");
		for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
		{
			if (i>0) felog.printf("---------------------------------------------------------------------------\n");

			FESurfacePairConstraint* pi = fem.SurfacePairConstraint(i);
			const char* sztype = pi->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			felog.printf("contact interface %d - Type: %s\n", i+1, sztype);
			FEParameterList& pl = pi->GetParameterList();
			print_parameter_list(pl);
		}
		felog.printf("\n\n");
	}

	if (fem.NonlinearConstraints() != 0)
	{
		felog.printf(" NONLINEAR CONSTRAINT DATA\n");
		felog.printf("===========================================================================\n");
		int NC = fem.NonlinearConstraints();
		for (int i = 0; i<NC; ++i)
		{
			FENLConstraint* plc = fem.NonlinearConstraint(i);
			if (dynamic_cast<FERigidJoint*>(plc))
			{
				FERigidJoint& rj = static_cast<FERigidJoint&>(*plc);
				FERigidBody& ra = *rigid.Object(rj.m_nRBa);
				FERigidBody& rb = *rigid.Object(rj.m_nRBb);
				felog.printf("rigid joint %d:\n", i+1);
				felog.printf("\tRigid body A                   : %d\n", ra.m_mat + 1);
				felog.printf("\tRigid body B                   : %d\n", rb.m_mat + 1);
				felog.printf("\tJoint                          : (%lg, %lg, %lg)\n", rj.m_q0.x, rj.m_q0.y, rj.m_q0.z);
				felog.printf("\tPenalty factor                 : %lg\n", rj.m_eps );
				felog.printf("\tAugmented Lagrangian tolerance : %lg\n", rj.m_atol);
				felog.printf("---------------------------------------------------------------------------\n");
			}
		}
		felog.printf("\n\n");
	}

	felog.printf(" LOADCURVE DATA\n");
	felog.printf("===========================================================================\n");
	for (int i = 0; i<fem.LoadCurves(); ++i)
	{
		if (i>0) felog.printf("---------------------------------------------------------------------------\n");
		felog.printf("%3d\n", i+1);
		FEDataLoadCurve* plc = dynamic_cast<FEDataLoadCurve*>(fem.GetLoadCurve(i));
		if (plc)
		{
			for (int j = 0; j<plc->Points(); ++j)
			{
				LOADPOINT pt = plc->LoadPoint(j);
				felog.printf("%10lg%10lg\n", pt.time, pt.value);
			}
		}
	}
	felog.printf("\n\n");

	felog.printf(" LINEAR SOLVER DATA\n");
	felog.printf("===========================================================================\n");
	felog.printf("\tSolver type ............................... : ");
	switch (fem.GetLinearSolverType())
	{
	case SKYLINE_SOLVER     : felog.printf("Skyline\n"           ); break;
	case PSLDLT_SOLVER      : felog.printf("PSLDLT\n"            ); break;
	case SUPERLU_SOLVER     : felog.printf("SuperLU\n"           ); break;
	case SUPERLU_MT_SOLVER  : felog.printf("SuperLU_MT\n"        ); break;
	case PARDISO_SOLVER     : felog.printf("Pardiso\n"           ); break;
	case WSMP_SOLVER        : felog.printf("WSMP\n"              ); break;
	case LU_SOLVER          : felog.printf("LUSolver\n"          ); break;
	case CG_ITERATIVE_SOLVER: felog.printf("Conjugate gradient\n"); break;
	case RCICG_SOLVER       : felog.printf("RCICG\n"             ); break;
	case FGMRES_SOLVER      : felog.printf("FGMRES\n"            ); break;
	case FGMRES_ILUT_SOLVER : felog.printf("FGMRES_ILUT\n"       ); break;
	case FGMRES_ILU0_SOLVER : felog.printf("FGMRES_ILU0\n"       ); break;
	case BIPN_SOLVER        : felog.printf("BIPN\n"              ); break;
	case HYPRE_GMRES        : felog.printf("HYPRE GMRES\n"       ); break;
	case STOKES_SOLVER      : felog.printf("Stokes\n"            ); break;
	default:
		assert(false);
		felog.printf("Unknown solver\n");
	}
	felog.printf("\tMatrix format ............................. : ");
	if (step.GetFESolver()->m_bsymm) felog.printf("symmetric\n");
	else felog.printf("unsymmetric\n");
	felog.printf("\n\n");
}
