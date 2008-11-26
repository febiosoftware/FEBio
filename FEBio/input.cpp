// input module

#include "stdafx.h"
#include "fem.h"
#include "FileImport.h"
#include "NikeImport.h"
#include "FEBioImport.h"
#include <string.h>

// Forward declarations
void Hello(FILE* fp);

//-----------------------------------------------------------------------------
//! This routine reads in an input file and performs some initialization stuff.
//! The rest of the initialization is done int Init

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

	// we only output this data to the log file and not the screen
	Logfile::MODE old_mode = m_log.SetMode(Logfile::FILE_ONLY);

	m_log.printf("%s\n\n", m_sztitle);

	m_log.printf(" FILES USED\n");
	m_log.printf("===========================================================================\n");
	m_log.printf("\tInput file : %s\n", m_szfile);
	m_log.printf("\tPlot file  : %s\n", m_szplot);
	m_log.printf("\tLog file   : %s\n", m_szlog);
	m_log.printf("\n\n");

	m_log.printf(" CONTROL DATA\n");
	m_log.printf("===========================================================================\n");
	m_log.printf("\tAnalysis type .................................. : %d\n", m_pStep->m_itype);
	m_log.printf("\t   eq.%2d: quasi-static\n", FE_STATIC);
	m_log.printf("\t   eq.%2d: dynamic\n", FE_DYNAMIC);
	m_log.printf("\t   eq.%2d: quasi-static + poro-elasticity\n", FE_STATIC_PORO);
	m_log.printf("\tNumber of materials ............................ : %d\n", Materials());
	m_log.printf("\tNumber of nodes ................................ : %d\n", m_mesh.Nodes() );
	m_log.printf("\tNumber of solid elements ....................... : %d\n", m_mesh.SolidElements());
	m_log.printf("\tNumber of shell elements ....................... : %d\n", m_mesh.ShellElements());
	m_log.printf("\tNumber of timesteps ............................ : %d\n", m_pStep->m_ntime);
	m_log.printf("\tTime step size ................................. : %lg\n", m_pStep->m_dt0);
	m_log.printf("\tAuto time stepper activated .................... : %s\n", (m_pStep->m_bautostep ? "yes" : "no"));
	if (m_pStep->m_bautostep)
	{
		m_log.printf("\t  Optimal nr of iterations ..................... : %d\n", m_pStep->m_iteopt);
		m_log.printf("\t  Minimum allowable step size .................. : %lg\n", m_pStep->m_dtmin);
		m_log.printf("\t  Maximum allowable step size .................. : %lg\n", m_pStep->m_dtmax);
	}
	m_log.printf("\tNumber of loadcurves ........................... : %d\n", LoadCurves());
	m_log.printf("\tNumber of displacement boundary conditions ..... : %d\n", m_DC.size());
	m_log.printf("\tNumber of pressure boundary cards .............. : %d\n", m_PC.size());
	m_log.printf("\tNumber of concentrated nodal forces ............ : %d\n", m_FC.size());
	m_log.printf("\tMax nr of stiffness reformations ............... : %d\n", m_pStep->m_psolver->m_maxref);
	m_log.printf("\tper time steps\n");
	m_log.printf("\tMax nr of Quasi-Newton iterations .............. : %d\n", m_pStep->m_psolver->m_maxups);
	m_log.printf("\tbetween stiffness matrix reformations\n"); 
	m_log.printf("\tDisplacement convergence tolerance ............. : %lg\n", m_pStep->m_psolver->m_Dtol);
	m_log.printf("\tEnergy convergence tolerance ................... : %lg\n", m_pStep->m_psolver->m_Etol);
	m_log.printf("\tResidual convergence tolerance ................. : %lg\n", m_pStep->m_psolver->m_Rtol);
	m_log.printf("\tLinesearch convergence tolerance ............... : %lg\n", m_pStep->m_psolver->m_LStol);
	m_log.printf("\tMinimum line search size ....................... : %lg\n", m_pStep->m_psolver->m_LSmin);
	m_log.printf("\tMaximum number of line search iterations ....... : %d\n", m_pStep->m_psolver->m_LSiter);
	m_log.printf("\tMax condition number ........................... : %lg\n", m_pStep->m_psolver->m_cmax);
	m_log.printf("\n\n");

	m_log.printf(" MATERIAL DATA\n");
	m_log.printf("===========================================================================\n");
	for (i=0; i<Materials(); ++i)
	{
		if (i>0) m_log.printf("---------------------------------------------------------------------------\n");
		m_log.printf("%3d - ", i+1);

		// get the material
		FEMaterial* pmat = GetMaterial(i);

		// get the material name and type string
		const char* szname = pmat->GetName();
		const char* sztype = pmat->GetTypeString();
		if (szname[0] == 0) szname = 0;

		// print type and name
		m_log.printf("%s", sztype);
		if (szname) m_log.printf(" (%s)", szname);
		m_log.printf("\n");

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
				case FE_PARAM_DOUBLE : m_log.printf("%s : %lg\n", sz, it->value<double>()); break;
				case FE_PARAM_INT    : m_log.printf("%s : %d\n" , sz, it->value<int   >()); break;
				case FE_PARAM_BOOL   : m_log.printf("%s : %d\n" , sz, it->value<bool  >()); break;
				case FE_PARAM_INTV   : 
				case FE_PARAM_DOUBLEV:
					{
						int n = it->m_ndim;
						m_log.printf("%s : ", sz);
						for (int k=0; k<n; ++k)
						{
							switch (it->m_itype)
							{
							case FE_PARAM_INTV   : m_log.printf("%d", it->pvalue<int   >()[k]); break;
							case FE_PARAM_DOUBLEV: m_log.printf("%lg", it->pvalue<double>()[k]); break;
							}
							if (k!=n-1) m_log.printf(","); else m_log.printf("\n");
						}
					}
					break;
				default:
					assert(false);
				}
			}
		}
	}
	m_log.printf("\n\n");

	m_log.printf(" LOADCURVE DATA\n");
	m_log.printf("===========================================================================\n");
	for (i=0; i<LoadCurves(); ++i)
	{
		if (i>0) m_log.printf("---------------------------------------------------------------------------\n");
		m_log.printf("%3d\n", i+1);
		FELoadCurve* plc = GetLoadCurve(i);
		for (j=0; j<plc->Points(); ++j)
		{
			LOADPOINT& pt = plc->LoadPoint(j);
			m_log.printf("%10lg%10lg\n", pt.time, pt.value);
		}
	}
	m_log.printf("\n\n");

	if (m_CI.size() > 0)
	{
		m_log.printf(" CONTACT INTERFACE DATA\n");
		m_log.printf("===========================================================================\n");
		for (i=0; i<m_CI.size(); ++i)
		{
			if (i>0) m_log.printf("---------------------------------------------------------------------------\n");

			FESlidingInterface *psi = dynamic_cast<FESlidingInterface*>(&m_CI[i]);
			if (psi)
			{
				m_log.printf("contact interface %d:\n", i+1);
				m_log.printf("\tType                           : sliding with gaps\n");
				m_log.printf("\tPenalty factor                 : %lg\n", psi->m_eps);
				m_log.printf("\tTwo-pass algorithm             : %s\n", (psi->npass==1? "off":"on"));
				m_log.printf("\tAugmented Lagrangian           : %s\n", (psi->m_blaugon? "on" : "off"));
				if (psi->m_blaugon)
					m_log.printf("\tAugmented Lagrangian tolerance : %lg\n", psi->m_atol);
			}

			FETiedInterface *pti = dynamic_cast<FETiedInterface*>(&m_CI[i]);
			if (pti)
			{
				m_log.printf("contact interface %d:\n", i+1);
				m_log.printf("\tType                           : tied\n");
				m_log.printf("\tPenalty factor                 : %lg\n", pti->m_eps);
				m_log.printf("\tAugmented Lagrangian tolerance : %lg\n", pti->m_atol);
			}

			FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(&m_CI[i]);
			if (pri)
			{
				m_log.printf("contact interface %d:\n", i+1);
				m_log.printf("\tType                           : rigid wall\n");
				m_log.printf("\tPenalty factor                 : %lg\n", pri->m_eps);
				m_log.printf("\tAugmented Lagrangian tolerance : %lg\n", pri->m_atol);
				if (dynamic_cast<FEPlane*>(pri->m_mp))
				{
					FEPlane* pp = dynamic_cast<FEPlane*>(pri->m_mp);
					m_log.printf("\tPlane equation:\n");
					double* a = pp->GetEquation();
					m_log.printf("\t\ta = %lg\n", a[0]);
					m_log.printf("\t\tb = %lg\n", a[1]);
					m_log.printf("\t\tc = %lg\n", a[2]);
					m_log.printf("\t\td = %lg\n", a[3]);
				}
				else if (dynamic_cast<FERigidSphere*>(pri->m_mp))
				{
					FERigidSphere* ps = dynamic_cast<FERigidSphere*>(pri->m_mp);
					m_log.printf("\trigid sphere:\n");
					m_log.printf("\t\tx ............... : %lg\n", ps->m_rc.x);
					m_log.printf("\t\ty ............... : %lg\n", ps->m_rc.y);
					m_log.printf("\t\tz ............... : %lg\n", ps->m_rc.z);
					m_log.printf("\t\tR ............... : %lg\n", ps->m_R);
				}
			}
		}
		m_log.printf("\n\n");
	}

	if (m_nrj > 0)
	{
		m_log.printf(" RIGID JOINT DATA\n");
		m_log.printf("===========================================================================\n");
		for (i=0; i<m_nrj; ++i)
		{
			if (i>0) m_log.printf("---------------------------------------------------------------------------\n");

			FERigidJoint& rj = m_RJ[i];
			m_log.printf("rigid joint %d:\n", i+1);
			m_log.printf("\tRigid body A                   : %d\n", m_RB[rj.m_nRBa].m_mat + 1);
			m_log.printf("\tRigid body B                   : %d\n", m_RB[rj.m_nRBb].m_mat + 1);
			m_log.printf("\tJoint                          : (%lg, %lg, %lg)\n", rj.m_q0.x, rj.m_q0.y, rj.m_q0.z);
			m_log.printf("\tPenalty factor                 : %lg\n", rj.m_eps );
			m_log.printf("\tAugmented Lagrangian tolerance : %lg\n", rj.m_atol);
		}
		m_log.printf("\n\n");
	}

	if (m_DE.size())
	{
		m_log.printf(" DISCRETE ELEMENT DATA\n");
		m_log.printf("===========================================================================\n");
		m_log.printf(" Nr of discrete elements : %d\n", m_DE.size());
		for (i=0; i<m_DE.size(); ++i)
		{
			FE_DISCRETE_ELEMENT& de = m_DE[i];
			m_log.printf(" discrete element %d:\n", i+1);
			m_log.printf("\tnodes : %d, %d\n", de.n1+1, de.n2+1);
			m_log.printf("\tE     : %lg\n", de.E);
		}
		m_log.printf("\n\n");
	}

	// reset log mode
	m_log.SetMode(old_mode);

}
