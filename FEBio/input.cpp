// input module

#include "stdafx.h"
#include "fem.h"
#include "FileImport.h"
#include "NikeImport.h"
#include "FEBioImport.h"
#include "FEPeriodicBoundary.h"
#include "FESurfaceConstraint.h"
#include "FEFacet2FacetSliding.h"
#include "log.h"
#include "LSDYNAPlotFile.h"
#include <string.h>

// Forward declarations
void Hello(FILE* fp);

//-----------------------------------------------------------------------------
//! This routine reads in an input file and performs some initialization stuff.
//! The rest of the initialization is done in Init

bool FEM::Input(const char* szfile)
{
	// create file readers
	FEFileImport* pin = 0;
	FENIKEImport  in_nike;
	FEFEBioImport in_febio;

	// check the extension of the filename
	// to determine what file reader we need.
	const char* szext = strrchr(szfile, '.');

	if (szext && ((strcmp(szext, ".xml")== 0) ||
		          (strcmp(szext, ".XML")== 0) ||
				  (strcmp(szext, ".feb")== 0) ))
	{
		// It's an FEBio file
		pin = &in_febio;
	}
	else
	{
		// It's a NIKE file
		pin = &in_nike;
	}

	// Load the file
	if (pin->Load(*this, szfile) == false)
	{
		char szerr[256];
		pin->GetErrorMessage(szerr);
		fprintf(stderr, szerr);

		return false;
	}

	// we're done reading
	return true;
}

//------------------------------------------------------------------------------
//! This function outputs the input data to the log file.

void FEM::EchoInput()
{
	// echo input
	int i, j;

	// get the logfile
	Logfile& log = GetLogfile();

	// we only output this data to the log file and not the screen
	Logfile::MODE old_mode = log.SetMode(Logfile::FILE_ONLY);

	// if for some reason the old_mode was set to NEVER, we should not output anything
	if (old_mode == Logfile::NEVER)
	{
		log.SetMode(old_mode);
		return;
	}

	bool bporo = m_pStep->m_nModule == FE_POROELASTIC;

	log.printf("%s\n\n", m_sztitle);

	log.printf(" FILES USED\n");
	log.printf("===========================================================================\n");
	log.printf("\tInput file : %s\n", m_szfile);
	log.printf("\tPlot file  : %s\n", m_szplot);
	log.printf("\tLog file   : %s\n", m_szlog);
	log.printf("\n\n");

	log.printf(" CONTROL DATA\n");
	log.printf("===========================================================================\n");
	log.printf("Module type ...................................... : %d\n", m_pStep->m_nModule);
	log.printf("\t   eq.%2d: solid mechanics\n", FE_SOLID);
	log.printf("\t   eq.%2d: poroelastic\n", FE_POROELASTIC);
	log.printf("\t   eq.%2d: heat transfer\n", FE_HEAT);
	log.printf("\tAnalysis type .................................. : %d\n", m_pStep->m_nanalysis);
	log.printf("\t   eq.%2d: quasi-static\n", FE_STATIC);
	log.printf("\t   eq.%2d: dynamic\n", FE_DYNAMIC);
	log.printf("\tPlane strain mode .............................. : %s\n", (m_bplane_strain? "yes" : "no"));
	log.printf("\tNumber of materials ............................ : %d\n", Materials());
	log.printf("\tNumber of nodes ................................ : %d\n", m_mesh.Nodes() );
	log.printf("\tNumber of solid elements ....................... : %d\n", m_mesh.SolidElements());
	log.printf("\tNumber of shell elements ....................... : %d\n", m_mesh.ShellElements());
	log.printf("\tNumber of truss elements ....................... : %d\n", m_mesh.TrussElements());
	log.printf("\tNumber of timesteps ............................ : %d\n", m_pStep->m_ntime);
	log.printf("\tTime step size ................................. : %lg\n", m_pStep->m_dt0);
	log.printf("\tAuto time stepper activated .................... : %s\n", (m_pStep->m_bautostep ? "yes" : "no"));
	if (m_pStep->m_bautostep)
	{
		log.printf("\t  Optimal nr of iterations ..................... : %d\n", m_pStep->m_iteopt);
		log.printf("\t  Minimum allowable step size .................. : %lg\n", m_pStep->m_dtmin);
		log.printf("\t  Maximum allowable step size .................. : %lg\n", m_pStep->m_dtmax);
	}
	log.printf("\tNumber of loadcurves ........................... : %d\n", LoadCurves());
	log.printf("\tNumber of displacement boundary conditions ..... : %d\n", m_DC.size());
	log.printf("\tNumber of pressure boundary cards .............. : %d\n", (m_psurf ? m_psurf->Elements() : 0));
	log.printf("\tNumber of constant traction boundary cards ..... : %d\n", (m_ptrac ? m_ptrac->Elements() : 0));
	log.printf("\tNumber of concentrated nodal forces ............ : %d\n", m_FC.size());
	log.printf("\tMax nr of stiffness reformations ............... : %d\n", m_pStep->m_psolver->m_bfgs.m_maxref);
	log.printf("\tper time steps\n");
	log.printf("\tMax nr of Quasi-Newton iterations .............. : %d\n", m_pStep->m_psolver->m_bfgs.m_maxups);
	log.printf("\tbetween stiffness matrix reformations\n");
	log.printf("\tDisplacement convergence tolerance ............. : %lg\n", m_pStep->m_psolver->m_Dtol);
	log.printf("\tEnergy convergence tolerance ................... : %lg\n", m_pStep->m_psolver->m_Etol);
	log.printf("\tResidual convergence tolerance ................. : %lg\n", m_pStep->m_psolver->m_Rtol);
	if (bporo) log.printf("\tFluid pressure convergence tolernace ........... : %lg\n", m_pStep->m_psolver->m_Ptol);
	log.printf("\tLinesearch convergence tolerance ............... : %lg\n", m_pStep->m_psolver->m_LStol);
	log.printf("\tMinimum line search size ....................... : %lg\n", m_pStep->m_psolver->m_LSmin);
	log.printf("\tMaximum number of line search iterations ....... : %d\n", m_pStep->m_psolver->m_LSiter);
	log.printf("\tMax condition number ........................... : %lg\n", m_pStep->m_psolver->m_bfgs.m_cmax);
	log.printf("\n\n");

	log.printf(" OUTPUT DATA\n");
	log.printf("===========================================================================\n");
	switch (m_pStep->m_nplot)
	{
	case FE_PLOT_NEVER      : log.printf("\tplot level ................................ : never\n"); break;
	case FE_PLOT_MAJOR_ITRS : log.printf("\tplot level ................................ : major iterations\n"); break;
	case FE_PLOT_MINOR_ITRS : log.printf("\tplot level ................................ : minor iterations\n"); break;
	case FE_PLOT_MUST_POINTS: log.printf("\tplot level ................................ : must points only\n"); break;
	}

	if (dynamic_cast<LSDYNAPlotFile*>(m_plot))
	{
		LSDYNAPlotFile& plt = *dynamic_cast<LSDYNAPlotFile*>(m_plot);
		log.printf("\tshell strains included .................... : %s\n", (plt.m_bsstrn? "Yes" : "No"));

		log.printf("\tplot file field data:\n");
		log.printf("\t\tdisplacement ...................... : ");
		switch (plt.m_nfield[0])
		{
		case PLOT_DISPLACEMENT: log.printf("displacement\n"); break;
		default: log.printf("???\n");
		}
		log.printf("\t\tvelocity .......................... : ");
		switch (plt.m_nfield[1])
		{
		case PLOT_NONE            : log.printf("none\n"); break;
		case PLOT_VELOCITY        : log.printf("velocity\n"); break;
		case PLOT_FLUID_FLUX      : log.printf("fluid flux\n"); break;
		case PLOT_CONTACT_TRACTION: log.printf("contact traction\n"); break;
		case PLOT_REACTION_FORCE  : log.printf("reaction force\n"); break;
		default: log.printf("???\n");
		}
		log.printf("\t\tacceleration ...................... : ");
		switch (plt.m_nfield[2])
		{
		case PLOT_NONE            : log.printf("none\n"); break;
		case PLOT_ACCELERATION    : log.printf("acceleration\n"); break;
		case PLOT_FLUID_FLUX      : log.printf("fluid flux\n"); break;
		case PLOT_CONTACT_TRACTION: log.printf("contact traction\n"); break;
		case PLOT_REACTION_FORCE  : log.printf("reaction force\n"); break;
		default: log.printf("???\n");
		}
		log.printf("\t\ttemperature........................ : ");
		switch (plt.m_nfield[3])
		{
		case PLOT_NONE            : log.printf("none\n"); break;
		case PLOT_FLUID_PRESSURE  : log.printf("fluid pressure\n"); break;
		case PLOT_CONTACT_PRESSURE: log.printf("contact pressure\n"); break;
		case PLOT_CONTACT_GAP     : log.printf("contact gap\n"); break;
		default: log.printf("???\n");
		}
	}

	log.printf("\n\n");
	log.printf(" MATERIAL DATA\n");
	log.printf("===========================================================================\n");
	for (i=0; i<Materials(); ++i)
	{
		if (i>0) log.printf("---------------------------------------------------------------------------\n");
		log.printf("%3d - ", i+1);

		// get the material
		FEMaterial* pmat = GetMaterial(i);

		// get the material name and type string
		const char* szname = pmat->GetName();
		const char* sztype = pmat->GetTypeString();
		if (szname[0] == 0) szname = 0;

		// print type and name
		log.printf("%s", sztype);
		if (szname) log.printf(" (%s)", szname);
		log.printf("\n");

		// print the parameter list
		auto_ptr<FEParameterList> pl(pmat->GetParameterList());
		int n = pl->Parameters();
		if (n > 0)
		{
			char sz[256];
			list<FEParam>::iterator it = pl->first();
			for (int j=0; j<n; ++j, ++it)
			{
				int l = strlen(it->m_szname);
				sprintf(sz, "\t%-*s %.*s : ", l, it->m_szname, 50-l, "..................................................");
				switch (it->m_itype)
				{
				case FE_PARAM_DOUBLE : log.printf("%s : %lg\n", sz, it->value<double>()); break;
				case FE_PARAM_INT    : log.printf("%s : %d\n" , sz, it->value<int   >()); break;
				case FE_PARAM_BOOL   : log.printf("%s : %d\n" , sz, it->value<bool  >()); break;
				case FE_PARAM_STRING : log.printf("%s : %s\n" , sz, it->cvalue()); break;
				case FE_PARAM_INTV   :
				case FE_PARAM_DOUBLEV:
					{
						int n = it->m_ndim;
						log.printf("%s : ", sz);
						for (int k=0; k<n; ++k)
						{
							switch (it->m_itype)
							{
							case FE_PARAM_INTV   : log.printf("%d", it->pvalue<int   >()[k]); break;
							case FE_PARAM_DOUBLEV: log.printf("%lg", it->pvalue<double>()[k]); break;
							}
							if (k!=n-1) log.printf(","); else log.printf("\n");
						}
					}
					break;
				default:
					assert(false);
				}
			}
		}
	}
	log.printf("\n\n");

	log.printf(" LOADCURVE DATA\n");
	log.printf("===========================================================================\n");
	for (i=0; i<LoadCurves(); ++i)
	{
		if (i>0) log.printf("---------------------------------------------------------------------------\n");
		log.printf("%3d\n", i+1);
		FELoadCurve* plc = GetLoadCurve(i);
		for (j=0; j<plc->Points(); ++j)
		{
			LOADPOINT& pt = plc->LoadPoint(j);
			log.printf("%10lg%10lg\n", pt.time, pt.value);
		}
	}
	log.printf("\n\n");

	if (m_CI.size() > 0)
	{
		log.printf(" CONTACT INTERFACE DATA\n");
		log.printf("===========================================================================\n");
		for (i=0; i<m_CI.size(); ++i)
		{
			if (i>0) log.printf("---------------------------------------------------------------------------\n");

			FESlidingInterface *psi = dynamic_cast<FESlidingInterface*>(&m_CI[i]);
			if (psi)
			{
				log.printf("contact interface %d:\n", i+1);
				log.printf("\tType                           : sliding with gaps\n");
				log.printf("\tPenalty factor                 : %lg\n", psi->m_eps);
				log.printf("\tAuto-penalty                   : %s\n", (psi->m_nautopen==2? "on" : "off"));
				log.printf("\tTwo-pass algorithm             : %s\n", (psi->m_npass==1? "off":"on"));
				log.printf("\tAugmented Lagrangian           : %s\n", (psi->m_blaugon? "on" : "off"));
				if (psi->m_blaugon)
					log.printf("\tAugmented Lagrangian tolerance : %lg\n", psi->m_atol);
				log.printf("\tmaster segments                : %d\n", (psi->m_ms.Elements()));
				log.printf("\tslave segments                 : %d\n", (psi->m_ss.Elements()));
			}

			FEFacet2FacetSliding* pf2f = dynamic_cast<FEFacet2FacetSliding*>(&m_CI[i]);
			if (pf2f)
			{
				log.printf("contact interface %d:\n", i+1);
				log.printf("\tType                           : facet-to-facet sliding\n");
				log.printf("\tPenalty factor                 : %lg\n", pf2f->m_epsn);
				log.printf("\tAuto-penalty                   : %s\n", (pf2f->m_bautopen? "on" : "off"));
				log.printf("\tTwo-pass algorithm             : %s\n", (pf2f->m_npass==1? "off": "on" ));
				log.printf("\tAugmented Lagrangian           : %s\n", (pf2f->m_blaugon ? "on" : "off"));
				if (pf2f->m_blaugon)
					log.printf("\tAugmented Lagrangian tolerance : %lg\n", pf2f->m_atol);
				log.printf("\tmaster segments                : %d\n", (pf2f->m_ms.Elements()));
				log.printf("\tslave segments                 : %d\n", (pf2f->m_ss.Elements()));
			}

			FETiedInterface *pti = dynamic_cast<FETiedInterface*>(&m_CI[i]);
			if (pti)
			{
				log.printf("contact interface %d:\n", i+1);
				log.printf("\tType                           : tied\n");
				log.printf("\tPenalty factor                 : %lg\n", pti->m_eps);
				log.printf("\tAugmented Lagrangian tolerance : %lg\n", pti->m_atol);
			}

			FEPeriodicBoundary *pbi = dynamic_cast<FEPeriodicBoundary*>(&m_CI[i]);
			if (pbi)
			{
				log.printf("contact interface %d:\n", i+1);
				log.printf("\tType                           : periodic\n");
				log.printf("\tPenalty factor                 : %lg\n", pbi->m_eps);
				log.printf("\tAugmented Lagrangian tolerance : %lg\n", pbi->m_atol);
			}

			FESurfaceConstraint *psc = dynamic_cast<FESurfaceConstraint*>(&m_CI[i]);
			if (psc)
			{
				log.printf("contact interface %d:\n", i+1);
				log.printf("\tType                           : surface constraint\n");
				log.printf("\tPenalty factor                 : %lg\n", psc->m_eps);
				log.printf("\tAugmented Lagrangian tolerance : %lg\n", psc->m_atol);
			}

			FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(&m_CI[i]);
			if (pri)
			{
				log.printf("contact interface %d:\n", i+1);
				log.printf("\tType                           : rigid wall\n");
				log.printf("\tPenalty factor                 : %lg\n", pri->m_eps);
				log.printf("\tAugmented Lagrangian tolerance : %lg\n", pri->m_atol);
				if (dynamic_cast<FEPlane*>(pri->m_mp))
				{
					FEPlane* pp = dynamic_cast<FEPlane*>(pri->m_mp);
					log.printf("\tPlane equation:\n");
					double* a = pp->GetEquation();
					log.printf("\t\ta = %lg\n", a[0]);
					log.printf("\t\tb = %lg\n", a[1]);
					log.printf("\t\tc = %lg\n", a[2]);
					log.printf("\t\td = %lg\n", a[3]);
				}
				else if (dynamic_cast<FERigidSphere*>(pri->m_mp))
				{
					FERigidSphere* ps = dynamic_cast<FERigidSphere*>(pri->m_mp);
					log.printf("\trigid sphere:\n");
					log.printf("\t\tx ............... : %lg\n", ps->m_rc.x);
					log.printf("\t\ty ............... : %lg\n", ps->m_rc.y);
					log.printf("\t\tz ............... : %lg\n", ps->m_rc.z);
					log.printf("\t\tR ............... : %lg\n", ps->m_R);
				}
			}
		}
		log.printf("\n\n");
	}

	if (m_nrj > 0)
	{
		log.printf(" RIGID JOINT DATA\n");
		log.printf("===========================================================================\n");
		for (i=0; i<m_nrj; ++i)
		{
			if (i>0) log.printf("---------------------------------------------------------------------------\n");

			FERigidJoint& rj = m_RJ[i];
			log.printf("rigid joint %d:\n", i+1);
			log.printf("\tRigid body A                   : %d\n", m_RB[rj.m_nRBa].m_mat + 1);
			log.printf("\tRigid body B                   : %d\n", m_RB[rj.m_nRBb].m_mat + 1);
			log.printf("\tJoint                          : (%lg, %lg, %lg)\n", rj.m_q0.x, rj.m_q0.y, rj.m_q0.z);
			log.printf("\tPenalty factor                 : %lg\n", rj.m_eps );
			log.printf("\tAugmented Lagrangian tolerance : %lg\n", rj.m_atol);
		}
		log.printf("\n\n");
	}

	if (m_DE.size())
	{
		log.printf(" DISCRETE ELEMENT DATA\n");
		log.printf("===========================================================================\n");
		log.printf(" Nr of discrete elements : %d\n", m_DE.size());
		for (i=0; i<m_DE.size(); ++i)
		{
			FE_DISCRETE_ELEMENT& de = m_DE[i];
			log.printf(" discrete element %d:\n", i+1);
			log.printf("\tnodes : %d, %d\n", de.n1+1, de.n2+1);
			log.printf("\tE     : %lg\n", de.E);
		}
		log.printf("\n\n");
	}

	log.printf(" LINEAR SOLVER DATA\n");
	log.printf("===========================================================================\n");
	log.printf("\tSolver type ............................... : ");
	if (m_nsolver == SKYLINE_SOLVER     ) log.printf("Skyline\n");
	if (m_nsolver == PSLDLT_SOLVER      ) log.printf("PSLDLT\n");
	if (m_nsolver == SUPERLU_SOLVER     ) log.printf("SuperLU\n");
	if (m_nsolver == SUPERLU_MT_SOLVER  ) log.printf("SuperLU_MT\n");
	if (m_nsolver == PARDISO_SOLVER     ) log.printf("Pardiso\n");
	if (m_nsolver == WSMP_SOLVER        ) log.printf("WSMP\n");
	if (m_nsolver == LU_SOLVER          ) log.printf("LUSolver\n");
	if (m_nsolver == CG_ITERATIVE_SOLVER) log.printf("Conjugate gradient\n");
	log.printf("\n\n");

	// reset log mode
	log.SetMode(old_mode);
}
