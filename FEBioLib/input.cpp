// input module

#include "stdafx.h"
#include "FEBioModel.h"
#include "FEBioXML/FileImport.h"
#include "FEBioXML/FEBioImport.h"
#include "FESlidingInterface.h"
#include "FETiedInterface.h"
#include "FERigidWallInterface.h"
#include "FEPeriodicBoundary.h"
#include "FESurfaceConstraint.h"
#include "FEFacet2FacetSliding.h"
#include "FEFacet2FacetTied.h"
#include "FERigidJoint.h"
#include "FEBiphasicSolver.h"
#include "FEBiphasicSoluteSolver.h"
#include <FECore/FERigidBody.h>
#include "FEBioPlot/LSDYNAPlotFile.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FECore/log.h"
#include <string.h>

//-----------------------------------------------------------------------------
// helper function to print a parameter to the logfile
void print_parameter(FEParam& p)
{
	char sz[512] = {0};
	int l = strlen(p.m_szname);
	sprintf(sz, "\t%-*s %.*s", l, p.m_szname, 50-l, "..................................................");
	switch (p.m_itype)
	{
	case FE_PARAM_DOUBLE : clog.printf("%s : %lg\n", sz, p.value<double>()); break;
	case FE_PARAM_INT    : clog.printf("%s : %d\n" , sz, p.value<int   >()); break;
	case FE_PARAM_BOOL   : clog.printf("%s : %d\n" , sz, (int) p.value<bool  >()); break;
	case FE_PARAM_STRING : clog.printf("%s : %s\n" , sz, p.cvalue()); break;
	case FE_PARAM_VEC3D  :
		{
			vec3d v = p.value<vec3d>();
			clog.printf("%s : %lg,%lg,%lg\n", sz, v.x, v.y, v.z);
		}
		break;
	case FE_PARAM_INTV   :
	case FE_PARAM_DOUBLEV:
		{
			int n = p.m_ndim;
			clog.printf("%s : ", sz);
			for (int k=0; k<n; ++k)
			{
				switch (p.m_itype)
				{
				case FE_PARAM_INTV   : clog.printf("%d", p.pvalue<int   >()[k]); break;
				case FE_PARAM_DOUBLEV: clog.printf("%lg", p.pvalue<double>()[k]); break;
				}
				if (k!=n-1) clog.printf(","); else clog.printf("\n");
			}
		}
		break;
	default:
		assert(false);
	}
}

//-----------------------------------------------------------------------------
// print the parameter list to the log file
void print_parameter_list(FEParameterList& pl)
{
	int n = pl.Parameters();
	if (n > 0)
	{
		list<FEParam>::iterator it = pl.first();
		for (int j=0; j<n; ++j, ++it) print_parameter(*it);
	}
}

//------------------------------------------------------------------------------
//! This function outputs the input data to the clog file.
void echo_input(FEBioModel& fem)
{
	// echo input
	int i, j;

	// we only output this data to the clog file and not the screen
	Logfile::MODE old_mode = clog.SetMode(Logfile::FILE_ONLY);

	// if for some reason the old_mode was set to NEVER, we should not output anything
	if (old_mode == Logfile::NEVER)
	{
		clog.SetMode(old_mode);
		return;
	}

	// get the analysis step
	FEAnalysisStep& step = dynamic_cast<FEAnalysisStep&>(*fem.GetCurrentStep());

	// get the FE mesh
	FEMesh& mesh = fem.GetMesh();

	// print title
	clog.printf("%s\n\n", fem.GetTitle());

	// print file info
	clog.printf(" FILES USED\n");
	clog.printf("===========================================================================\n");
	clog.printf("\tInput file : %s\n", fem.GetInputFileName());
	clog.printf("\tPlot file  : %s\n", fem.GetPlotFileName());
	clog.printf("\tLog file   : %s\n", fem.GetLogfileName());
	clog.printf("\n\n");

	// print control info
	clog.printf(" CONTROL DATA\n");
	clog.printf("===========================================================================\n");
	const char* szmod = 0;
	switch (step.GetType())
	{
	case FE_SOLID         : szmod = "solid mechanics"; break;
	case FE_EXPLICIT_SOLID: szmod = "explicit solid mechanics"; break;
	case FE_BIPHASIC      : szmod = "poroelastic"    ; break;
	case FE_HEAT          : szmod = "heat transfer"  ; break;
	case FE_POROSOLUTE    : szmod = "biphasic-solute"; break;
	case FE_LINEAR_SOLID  : szmod = "linear solid"   ; break;
	case FE_HEAT_SOLID    : szmod = "heat solid"     ; break;
	default:
		szmod = "unknown";
		assert(false);
	}
	clog.printf("\tModule type .................................... : %s\n", szmod);

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
	clog.printf("\tAnalysis type .................................. : %s\n", szan);

	clog.printf("\tPlane strain mode .............................. : %s\n", (fem.m_nplane_strain != -1? "yes" : "no"));
	clog.printf("\tNumber of materials ............................ : %d\n", fem.Materials());
	clog.printf("\tNumber of nodes ................................ : %d\n", mesh.Nodes() );
	clog.printf("\tNumber of solid elements ....................... : %d\n", mesh.SolidElements());
	clog.printf("\tNumber of shell elements ....................... : %d\n", mesh.ShellElements());
	clog.printf("\tNumber of truss elements ....................... : %d\n", mesh.TrussElements());
	if (step.m_ntime > 0)
		clog.printf("\tNumber of timesteps ............................ : %d\n", step.m_ntime);
	else
		clog.printf("\tFinal time ..................................... : %lg\n", step.m_final_time);

	clog.printf("\tTime step size ................................. : %lg\n", step.m_dt0);
	clog.printf("\tAuto time stepper activated .................... : %s\n", (step.m_bautostep ? "yes" : "no"));
	if (step.m_bautostep)
	{
		clog.printf("\t  Optimal nr of iterations ..................... : %d\n", step.m_iteopt);
		clog.printf("\t  Minimum allowable step size .................. : %lg\n", step.m_dtmin);
		clog.printf("\t  Maximum allowable step size .................. : %lg\n", step.m_dtmax);
	}
	clog.printf("\tNumber of loadcurves ........................... : %d\n", fem.LoadCurves());
	clog.printf("\tNumber of displacement boundary conditions ..... : %d\n", fem.PrescribedBCs());
//	clog.printf("\tNumber of pressure boundary cards .............. : %d\n", (fem.m_psurf ? fem.m_psurf->Surface().Elements() : 0));
//	clog.printf("\tNumber of constant traction boundary cards ..... : %d\n", (fem.m_ptrac ? fem.m_ptrac->Surface().Elements() : 0));
//	clog.printf("\tNumber of fluid flux boundary cards .............: %d\n", (fem.m_fsurf ? fem.m_fsurf->Surface().Elements() : 0));
	clog.printf("\tNumber of concentrated nodal forces ............ : %d\n", fem.NodalLoads());
	clog.printf("\tMax nr of stiffness reformations ............... : %d\n", step.m_psolver->m_bfgs.m_maxref);
	clog.printf("\tper time steps\n");
	clog.printf("\tMax nr of Quasi-Newton iterations .............. : %d\n", step.m_psolver->m_bfgs.m_maxups);
	clog.printf("\tbetween stiffness matrix reformations\n");
	FESolidSolver* ps = dynamic_cast<FESolidSolver*>(step.m_psolver);
	if (ps)
	{
		clog.printf("\tDisplacement convergence tolerance ............. : %lg\n", ps->m_Dtol);
		clog.printf("\tEnergy convergence tolerance ................... : %lg\n", ps->m_Etol);
		clog.printf("\tResidual convergence tolerance ................. : %lg\n", ps->m_Rtol);
		clog.printf("\tMinimal residual value ......................... : %lg\n", ps->m_Rmin);
	}
	FEBiphasicSolver* pps = dynamic_cast<FEBiphasicSolver*>(step.m_psolver);
	if (pps) clog.printf("\tFluid pressure convergence tolerance ........... : %lg\n", pps->m_Ptol);
	FEBiphasicSoluteSolver* pss = dynamic_cast<FEBiphasicSoluteSolver*>(step.m_psolver);
	if (pss) clog.printf("\tSolute concentration convergence tolerance ..... : %lg\n", pss->m_Ctol);
	clog.printf("\tLinesearch convergence tolerance ............... : %lg\n", step.m_psolver->m_bfgs.m_LStol );
	clog.printf("\tMinimum line search size ....................... : %lg\n", step.m_psolver->m_bfgs.m_LSmin );
	clog.printf("\tMaximum number of line search iterations ....... : %d\n" , step.m_psolver->m_bfgs.m_LSiter);
	clog.printf("\tMax condition number ........................... : %lg\n", step.m_psolver->m_bfgs.m_cmax  );
	clog.printf("\n\n");

	// print output data
	clog.printf(" OUTPUT DATA\n");
	clog.printf("===========================================================================\n");
	switch (step.m_nplot)
	{
	case FE_PLOT_NEVER      : clog.printf("\tplot level ................................ : never\n"); break;
	case FE_PLOT_MAJOR_ITRS : clog.printf("\tplot level ................................ : major iterations\n"); break;
	case FE_PLOT_MINOR_ITRS : clog.printf("\tplot level ................................ : minor iterations\n"); break;
	case FE_PLOT_MUST_POINTS: clog.printf("\tplot level ................................ : must points only\n"); break;
	case FE_PLOT_FINAL      : clog.printf("\tplot level ................................ : final state\n"); break;
	}

	PlotFile* pplt = fem.GetPlotFile();
	if (dynamic_cast<LSDYNAPlotFile*>(pplt))
	{
		LSDYNAPlotFile& plt = *dynamic_cast<LSDYNAPlotFile*>(pplt);

		clog.printf("\tplotfile format ........................... : LSDYNA\n");
		clog.printf("\tshell strains included .................... : %s\n", (plt.m_bsstrn? "Yes" : "No"));

		clog.printf("\tplot file field data:\n");
		clog.printf("\t\tdisplacement ...................... : ");
		switch (plt.m_nfield[0])
		{
		case PLOT_DISPLACEMENT: clog.printf("displacement\n"); break;
		default: clog.printf("???\n");
		}
		clog.printf("\t\tvelocity .......................... : ");
		switch (plt.m_nfield[1])
		{
		case PLOT_NONE            : clog.printf("none\n"); break;
		case PLOT_VELOCITY        : clog.printf("velocity\n"); break;
		case PLOT_FLUID_FLUX      : clog.printf("fluid flux\n"); break;
		case PLOT_CONTACT_TRACTION: clog.printf("contact traction\n"); break;
		case PLOT_REACTION_FORCE  : clog.printf("reaction force\n"); break;
		default: clog.printf("???\n");
		}
		clog.printf("\t\tacceleration ...................... : ");
		switch (plt.m_nfield[2])
		{
		case PLOT_NONE            : clog.printf("none\n"); break;
		case PLOT_ACCELERATION    : clog.printf("acceleration\n"); break;
		case PLOT_FLUID_FLUX      : clog.printf("fluid flux\n"); break;
		case PLOT_CONTACT_TRACTION: clog.printf("contact traction\n"); break;
		case PLOT_REACTION_FORCE  : clog.printf("reaction force\n"); break;
		default: clog.printf("???\n");
		}
		clog.printf("\t\ttemperature........................ : ");
		switch (plt.m_nfield[3])
		{
		case PLOT_NONE            : clog.printf("none\n"); break;
		case PLOT_FLUID_PRESSURE  : clog.printf("fluid pressure\n"); break;
		case PLOT_CONTACT_PRESSURE: clog.printf("contact pressure\n"); break;
		case PLOT_CONTACT_GAP     : clog.printf("contact gap\n"); break;
		default: clog.printf("???\n");
		}
	}

	if (dynamic_cast<FEBioPlotFile*>(pplt))
	{
		FEBioPlotFile* pf = dynamic_cast<FEBioPlotFile*>(pplt);
		clog.printf("\tplotfile format ........................... : FEBIO\n");

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
				clog.printf("\t\t%s:\n", szn);
				list<FEBioPlotFile::DICTIONARY_ITEM>::const_iterator it;
				for (it = pl->begin(); it != pl->end(); ++it)
				{
					const char* szt = 0;
					switch (it->m_ntype)
					{
					case PLT_FLOAT : szt = "float"; break;
					case PLT_VEC3F : szt = "vec3f"; break;
					case PLT_MAT3FS: szt = "mat3f"; break;
					}

					const char* szf = 0;
					switch (it->m_nfmt)
					{
					case FMT_NODE: szf = "NODE"; break;
					case FMT_ITEM: szf = "ITEM"; break;
					case FMT_MULT: szf = "COMP"; break;
					}

					clog.printf("\t\t\t%-20s (type = %5s, format = %4s)\n", it->m_szname, szt, szf);
				}
			}
		}
	}

	FEBioKernel& febio = FEBioKernel::GetInstance();

	// material data
	clog.printf("\n\n");
	clog.printf(" MATERIAL DATA\n");
	clog.printf("===========================================================================\n");
	for (i=0; i<fem.Materials(); ++i)
	{
		if (i>0) clog.printf("---------------------------------------------------------------------------\n");
		clog.printf("%3d - ", i+1);

		// get the material
		FEMaterial* pmat = fem.GetMaterial(i);

		// get the material name and type string
		const char* szname = pmat->GetName();
		const char* sztype = febio.GetTypeStr<FEMaterial>(pmat);
		if (szname[0] == 0) szname = 0;

		// print type and name
		clog.printf("%s", (szname?szname:"unknown"));
		clog.printf(" (type: %s)", sztype);
		clog.printf("\n");

		// print the parameter list
		FEParameterList& pl = pmat->GetParameterList();
		print_parameter_list(pl);
	}
	clog.printf("\n\n");

	if (fem.HasBodyLoads())
	{
		clog.printf(" BODY LOAD DATA\n");
		clog.printf("===========================================================================\n");
		for (i=0; i<fem.BodyLoads(); ++i)
		{
			if (i>0) clog.printf("---------------------------------------------------------------------------\n");
			clog.printf("%3d - ", i+1);

			// get the body load
			FEBodyLoad* pbl = fem.GetBodyLoad(i);

			// get the type string
			const char* sztype = febio.GetTypeStr<FEBodyLoad>(pbl);
			if (sztype == 0) sztype = "unknown";
			clog.printf(" Type: %s\n", sztype);

			// print the parameter list
			FEParameterList& pl = pbl->GetParameterList();
			print_parameter_list(pl);
		}
		clog.printf("\n\n");
	}

	if (fem.ContactInterfaces() > 0)
	{
		clog.printf(" CONTACT INTERFACE DATA\n");
		clog.printf("===========================================================================\n");
		for (i=0; i<fem.ContactInterfaces(); ++i)
		{
			if (i>0) clog.printf("---------------------------------------------------------------------------\n");

			FEContactInterface* pi = fem.ContactInterface(i);
			const char* sztype = febio.GetTypeStr<FEContactInterface>(pi);
			if (sztype == 0) sztype = "unknown";
			clog.printf("contact interface %d - Type: %s\n", i+1, sztype);
			FEParameterList& pl = pi->GetParameterList();
			print_parameter_list(pl);
		}
		clog.printf("\n\n");
	}

	if (fem.NonlinearConstraints() != 0)
	{
		clog.printf(" NONLINEAR CONSTRAINT DATA\n");
		clog.printf("===========================================================================\n");
		int NC = fem.NonlinearConstraints();
		for (i=0; i<NC; ++i)
		{
			FENLConstraint* plc = fem.NonlinearConstraint(i);
			if (dynamic_cast<FERigidJoint*>(plc))
			{
				FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*plc);
				FERigidBody& ra = dynamic_cast<FERigidBody&>(*fem.Object(rj.m_nRBa));
				FERigidBody& rb = dynamic_cast<FERigidBody&>(*fem.Object(rj.m_nRBb));
				clog.printf("rigid joint %d:\n", i+1);
				clog.printf("\tRigid body A                   : %d\n", ra.m_mat + 1);
				clog.printf("\tRigid body B                   : %d\n", rb.m_mat + 1);
				clog.printf("\tJoint                          : (%lg, %lg, %lg)\n", rj.m_q0.x, rj.m_q0.y, rj.m_q0.z);
				clog.printf("\tPenalty factor                 : %lg\n", rj.m_eps );
				clog.printf("\tAugmented Lagrangian tolerance : %lg\n", rj.m_atol);
				clog.printf("---------------------------------------------------------------------------\n");
			}
		}
		clog.printf("\n\n");
	}

	clog.printf(" LOADCURVE DATA\n");
	clog.printf("===========================================================================\n");
	for (i=0; i<fem.LoadCurves(); ++i)
	{
		if (i>0) clog.printf("---------------------------------------------------------------------------\n");
		clog.printf("%3d\n", i+1);
		FELoadCurve* plc = fem.GetLoadCurve(i);
		for (j=0; j<plc->Points(); ++j)
		{
			LOADPOINT& pt = plc->LoadPoint(j);
			clog.printf("%10lg%10lg\n", pt.time, pt.value);
		}
	}
	clog.printf("\n\n");

	clog.printf(" LINEAR SOLVER DATA\n");
	clog.printf("===========================================================================\n");
	clog.printf("\tSolver type ............................... : ");
	switch (fem.m_nsolver)
	{
	case SKYLINE_SOLVER     : clog.printf("Skyline\n"           ); break;
	case PSLDLT_SOLVER      : clog.printf("PSLDLT\n"            ); break;
	case SUPERLU_SOLVER     : clog.printf("SuperLU\n"           ); break;
	case SUPERLU_MT_SOLVER  : clog.printf("SuperLU_MT\n"        ); break;
	case PARDISO_SOLVER     : clog.printf("Pardiso\n"           ); break;
	case WSMP_SOLVER        : clog.printf("WSMP\n"              ); break;
	case LU_SOLVER          : clog.printf("LUSolver\n"          ); break;
	case CG_ITERATIVE_SOLVER: clog.printf("Conjugate gradient\n"); break;
	case RCICG_SOLVER       : clog.printf("RCICG\n"             ); break;
	default:
		assert(false);
		clog.printf("Unknown solver\n");
	}
	clog.printf("\n\n");

	// reset clog mode
	clog.SetMode(old_mode);
}
