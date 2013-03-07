#include "stdafx.h"
#include "FEBioModel.h"
#include <FECore/FERigid.h>
#include <FECore/FERigidBody.h>
#include "FERigidJoint.h"
#include "FEDiscreteMaterial.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FESlidingInterface.h"
#include "FETiedInterface.h"
#include "FETiedBiphasicInterface.h"
#include "FERigidWallInterface.h"
#include "FEFacet2FacetSliding.h"
#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FEPeriodicBoundary.h"
#include "FESurfaceConstraint.h"
#include "FETransverselyIsotropic.h"
#include "FEPressureLoad.h"
#include "FETractionLoad.h"
#include "FEFluidFlux.h"
#include "FEPoroTraction.h"
#include "FESoluteFlux.h"
#include "FEHeatFlux.h"
#include "FEConvectiveHeatFlux.h"
#include "FEAnalysisStep.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FEElasticTrussDomain.h"
#include "FEHeatSolidDomain.h"
#include "FEDiscreteSpringDomain.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBiphasicSoluteDomain.h"
#include "FEUDGHexDomain.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FEUT4Domain.h"
#include "FEConstBodyForce.h"
#include "FEPointConstraint.h"
#include "FEAugLagLinearConstraint.h"
#include <FECore/FERigidBody.h>
#include "NodeDataRecord.h"
#include "ElementDataRecord.h"
#include "RigidBodyDataRecord.h"
#include "FEBioPlot/LSDYNAPlotFile.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEBioXML/FEBioImport.h"
#include "FEMultiphasic.h"
#include "FEViscoElasticMaterial.h"
#include "FEUncoupledViscoElasticMaterial.h"
#include "FEElasticMultigeneration.h"
#include "FESlidingInterfaceBW.h"
#include "FECore/log.h"
#include "version.h"

//-----------------------------------------------------------------------------
// echo the input data to the log file
extern void echo_input(FEBioModel& fem);

//-----------------------------------------------------------------------------
// Constructor of FEBioModel class.
FEBioModel::FEBioModel()
{
	// --- Direct Solver Data ---
	// set the default linear solver
	// TODO: I don't want to do this here. In fact, I want
	//       to be able to use different solvers even if 
	//       when these directives are not declared.
#ifdef PARDISO
	m_nsolver = PARDISO_SOLVER;
#elif PSLDLT
	m_nsolver = PSLDLT_SOLVER;
#else
	m_nsolver = SKYLINE_SOLVER;
#endif

	// --- I/O-Data ---
	strcpy(m_szplot, "n3plot");
	strcpy(m_szlog , "n3log" );
	strcpy(m_szdump, "n3dump");
	m_sztitle[0] = 0;
	m_debug = false;
	m_becho = true;
	m_plot = 0;
}

//-----------------------------------------------------------------------------
//! Sets the extension of the plot file name.
void FEBioModel::SetPlotFileNameExtension(const char *szext)
{
	char* ch = strrchr(m_szplot, '.');
	if (ch) *ch = 0;
	strcat(m_szplot, szext);
}

//-----------------------------------------------------------------------------
//! Sets the name of the FEBio input file
void FEBioModel::SetInputFilename(const char* szfile)
{ 
	strcpy(m_szfile, szfile); 
	m_szfile_title = strrchr(m_szfile, '/');
	if (m_szfile_title == 0) 
	{
		m_szfile_title = strchr(m_szfile, '\\'); 
		if (m_szfile_title == 0) m_szfile_title = m_szfile; else ++m_szfile_title;
	}
	else ++m_szfile_title;
}

//-----------------------------------------------------------------------------
//! This function resets the FEM data so that a new run can be done.
//! This routine is called from the optimization routine.

bool FEBioModel::Reset()
{
	int i;

	// initialize materials
	FEMaterial* pmat;
	
	for (i=0; i<Materials(); ++i)
	{
		pmat = GetMaterial(i);
		pmat->Init();
	}

	// reset mesh data
	m_mesh.Reset();

	// reset object data
	int nrb = m_Obj.size();
	for (i=0; i<nrb; ++i) m_Obj[i]->Reset();

	// set up rigid joints
	if (!m_NLC.empty())
	{
		int NC = (int) m_NLC.size();
		for (i=0; i<NC; ++i)
		{
			FENLConstraint* plc = m_NLC[i];
			if (dynamic_cast<FERigidJoint*>(plc))
			{
				FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*plc);
				rj.m_F = vec3d(0,0,0);

				FERigidBody& ra = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBa]);
				FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBb]);

				rj.m_qa0 = rj.m_q0 - ra.m_r0;
				rj.m_qb0 = rj.m_q0 - rb.m_r0;
			}
		}
	}

	// set the start time
	m_ftime = 0;
	m_ftime0 = 0;

	// set first time step
	m_pStep = m_Step[0];
	m_nStep = 0;
	m_pStep->m_dt = m_pStep->m_dt0;
	m_pStep->m_ntotref    = 0;		// total nr of stiffness reformations
	m_pStep->m_ntotiter   = 0;		// total nr of non-linear iterations
	m_pStep->m_ntimesteps = 0;		// time steps completed
	m_pStep->m_ntotrhs    = 0;		// total nr of right hand side evaluations

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) m_plot = new LSDYNAPlotFile;

		if (m_plot->Open(*this, m_szplot) == false)
		{
			clog.printf("ERROR : Failed creating PLOT database\n");
			return false;
		}
	}

	// Since it is assumed that for the first timestep
	// there are no loads or initial displacements, the case n=0 is skipped.
	// Therefor we can output those results here.
	// Offcourse we should actually check if this is indeed
	// the case, otherwise we should also solve for t=0
	if (m_pStep->m_nplot != FE_PLOT_NEVER) m_plot->Write(*this);
/*
	// reset the log file
	if (!log.is_valid())
	{
		log.open(m_szlog);

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) log.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Hello();
	}
*/
	// do the callback
	DoCallback();

	// All data is reset successfully
	return true;
}

//-----------------------------------------------------------------------------
//! Find a BC based on its ID. This is needed for restarts.
FEBoundaryCondition* FEBioModel::FindBC(int nid)
{
	int i;
	for (i=0; i<(int) m_DC.size(); ++i) if (m_DC[i]->GetID() == nid) return m_DC[i];

	for (i=0; i<(int) m_FC.size(); ++i) if (m_FC[i]->GetID() == nid) return m_FC[i];

	for (i=0; i<(int) m_SL.size(); ++i) if (m_SL[i]->GetID() == nid) return m_SL[i];

	for (i=0; i<(int) m_RDC.size(); ++i) if (m_RDC[i]->GetID() == nid) return m_RDC[i];

	for (i=0; i<(int) m_RFC.size(); ++i) if (m_RFC[i]->GetID() == nid) return m_RFC[i];

	for (i=0; i<(int) m_RN.size(); ++i) if (m_RN[i]->GetID() == nid) return m_RN[i];

	return 0;
}

//=============================================================================
//    P A R A M E T E R   F U N C T I O N S
//=============================================================================

//-----------------------------------------------------------------------------
//! Helper function returning a pointer to the named variable in a solid mixture

double* FindSolidMixtureParameter(const char* szvar, const int index, FEElasticMixture* pme)
{
	char* ch = strchr((char*)szvar, '.');
	if (ch == 0) return 0;
	*ch = 0;
	const char* szvar2 = ch+1;
	
	int NMAT = pme->Materials();
	for (int i=0; i<NMAT; ++i) 
	{
		FEElasticMaterial* pmi = pme->GetMaterial(i);
		if (strcmp(szvar, pmi->GetName()) == 0)
		{
			// search the nested material parameter list
			FEParameterList& pl = pmi->GetParameterList();
			FEParam* pp = pl.Find(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
	}
	// no match found
	return 0;
}

//-----------------------------------------------------------------------------
//! Helper function returning a pointer to the named variable in a uncoupled solid mixture

double* FindUncoupledSolidMixtureParameter(const char* szvar, const int index, FEUncoupledElasticMixture* pme)
{
	char* ch = strchr((char*)szvar, '.');
	if (ch == 0) return 0;
	*ch = 0;
	const char* szvar2 = ch+1;
	
	FEMaterial* pmat;
	for (int i=0; i<(int) pme->m_pMat.size(); ++i) {
		if (strcmp(szvar, pme->m_pMat[i]->GetName()) == 0)
		{
			// search the nested material parameter list
			pmat = pme->m_pMat[i];
			FEParameterList& pl = pmat->GetParameterList();
			FEParam* pp = pl.Find(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
	}
	// no match found
	return 0;
}

//-----------------------------------------------------------------------------
//! Helper function returning a pointer to the named variable in an elastic material
double* FindElasticMaterialParameter(const char* szvar, int index, FEElasticMaterial* pme)
{
	FEParam* pp = pme->GetParameter(szvar);
	if (pp) return pp->pvalue<double>(index);
	// if material is solid mixture, check individual solid materials
	FEElasticMixture* pmm = dynamic_cast<FEElasticMixture*>(pme);
	if (pmm) return FindSolidMixtureParameter(szvar, index, pmm);
	// if this material is a viscoelastic material, check its elastic solid
	FEViscoElasticMaterial* pmv = dynamic_cast<FEViscoElasticMaterial*>(pme);
	if (pmv)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar3 = ch+1;
				
		if (strcmp(szvar, "elastic") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = pmv->m_pBase;
			FEParam* pp = pme->GetParameter(szvar3);
			if (pp) return pp->pvalue<double>(index);
			// if material is solid mixture, check individual solid materials
			FEElasticMixture* pmm = dynamic_cast<FEElasticMixture*>(pme);
			if (pmm) return FindSolidMixtureParameter(szvar3, index, pmm);
			else return 0;
		}
		else return 0;
	}

	// if this material is an uncoupled viscoelastic material, check its elastic solid
	FEUncoupledViscoElasticMaterial* puv = dynamic_cast<FEUncoupledViscoElasticMaterial*>(pme);
	if (puv)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;
		
		if (strcmp(szvar, "elastic") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = puv->m_pBase;
			FEParam* pp = pme->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			// if material is an uncoupled solid mixture, check individual solid materials
			FEUncoupledElasticMixture* pmm = dynamic_cast<FEUncoupledElasticMixture*>(pme);
			if (pmm) return FindUncoupledSolidMixtureParameter(szvar2, index, pmm);
			else return 0;
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable

//! This function returns a pointer to a named variable. Currently, we only
//! support names of the form:
//!		material_name.parameter_name
//!		material_name.elastic.parameter_name (nested material)
//!		material_name.solid_name.parameter_name (solid mixture)
//!		material_name.solid.parameter_name (biphasic material)
//!		material_name.permeability.parameter_name (biphasic material)
//!		material_name.solid.solid_name.parameter_name (biphasic material with solid mixture)
//! The 'material_name' is a user defined name for a material.
//! The 'parameter_name' is the predefined name of the variable.
//! The keywords 'elastic', 'solid', and 'permeability' must appear as shown.
//! \todo perhaps I should use XPath to refer to material parameters ?

double* FEBioModel::FindParameter(const char* szparam)
{
	int i, nmat;

	char szname[256];
	strcpy(szname, szparam);

	// get the material and parameter name
	char* ch = strchr((char*)szname, '.');
	if (ch == 0) return 0;
	*ch = 0;
	const char* szmat = szname;
	const char* szvar = ch+1;

	// find the material with the same name
	FEMaterial* pmat;

	for (i=0; i<Materials(); ++i)
	{
		pmat = GetMaterial(i);
		nmat = i;

		if (strcmp(szmat, pmat->GetName()) == 0)
		{
			break;
		}

		pmat = 0;
	}

	// make sure we found a material with the same name
	if (pmat == 0) return false;

	// if the variable is a vector, then we require an index
	char* szarg = strchr((char*) szvar, '[');
	int index = 0;
	if (szarg)
	{
		*szarg = 0; szarg++;
		const char* ch = strchr(szarg, ']');
		assert(ch);
		index = atoi(szarg) - 1;	// index is one-based for user
	}

	// find the parameter
	FEParam* pp = pmat->GetParameter(szvar);
	if (pp) return pp->pvalue<double>(index);

	// if material is solid mixture, check individual solid materials
	FEElasticMixture* pme = dynamic_cast<FEElasticMixture*>(pmat);
	if (pme) return FindSolidMixtureParameter(szvar, index, pme);
	FEUncoupledElasticMixture* pmu = dynamic_cast<FEUncoupledElasticMixture*>(pmat);
	if (pmu) return FindUncoupledSolidMixtureParameter(szvar, index, pmu);
	
	// if this material is a viscoelastic material, check its elastic solid
	FEViscoElasticMaterial* pmv = dynamic_cast<FEViscoElasticMaterial*>(pmat);
	if (pmv)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;
		
		if (strcmp(szvar, "elastic") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = pmv->m_pBase;
			FEParam* pp = pme->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			// if material is solid mixture, check individual solid materials
			FEElasticMixture* pmm = dynamic_cast<FEElasticMixture*>(pme);
			if (pmm) return FindSolidMixtureParameter(szvar2, index, pmm);
			else return 0;
		}
	}
	
	// if this material is an uncoupled viscoelastic material, check its elastic solid
	FEUncoupledViscoElasticMaterial* puv = dynamic_cast<FEUncoupledViscoElasticMaterial*>(pmat);
	if (puv)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;
		
		if (strcmp(szvar, "elastic") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = puv->m_pBase;
			FEParam* pp = pme->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			// if material is an uncoupled solid mixture, check individual solid materials
			FEUncoupledElasticMixture* pmm = dynamic_cast<FEUncoupledElasticMixture*>(pme);
			if (pmm) return FindUncoupledSolidMixtureParameter(szvar2, index, pmm);
			else return 0;
		}
	}

	// if this material is a biphasic material, check solid and permeability materials
	FEBiphasic* pmb = dynamic_cast<FEBiphasic*>(pmat);
	if (pmb)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;
		
		if (strcmp(szvar, "solid") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = pmb->m_pSolid;
			return FindElasticMaterialParameter(szvar2, index, pme);
		}
		else if (strcmp(szvar, "permeability") == 0)
		{
			// search the nested material parameter list
			pmat = pmb->m_pPerm;
			FEParam* pp = pmat->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
	}

	// if this material is a biphasic-solute material, check solid, permeability,
	// osmotic_coefficient, and solute materials
	FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*>(pmat);
	if (pbs)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;
		
		if (strcmp(szvar, "solid") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = pbs->m_pSolid;
			return FindElasticMaterialParameter(szvar2, index, pme);
		}
		else if (strcmp(szvar, "permeability") == 0)
		{
			// search the nested material parameter list
			pmat = pbs->m_pPerm;
			FEParam* pp = pmat->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
		else if (strcmp(szvar, "osmotic_coefficient") == 0)
		{
			// search the nested material parameter list
			pmat = pbs->m_pOsmC;
			FEParam* pp = pmat->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
		else if (strcmp(szvar, "solute") == 0)
		{
			char* ch = strchr((char*)szvar2, '.');
			if (ch == 0) return 0;
			*ch = 0;
			const char* szvar3 = ch+1;
			
            FESolute* pms = pbs->m_pSolute;
            if (strcmp(szvar2, "diffusivity") == 0)
            {
                // search the nested material parameter list
                FESoluteDiffusivity* pmd = pms->m_pDiff;
				FEParam* pp = pmd->GetParameter(szvar3);
                if (pp) return pp->pvalue<double>(index);
                else return 0;
            }
            else if (strcmp(szvar2, "solubility") == 0)
            {
                // search the nested material parameter list
                FESoluteSolubility* pmd = pms->m_pSolub;
				FEParam* pp = pmd->GetParameter(szvar3);
                if (pp) return pp->pvalue<double>(index);
                else return 0;
            }
			// no match found
			return 0;
		}
	}
	
	// if this material is a triphasic material, check solid, permeability,
	// osmotic_coefficient, and solute materials
	FETriphasic* ptp = dynamic_cast<FETriphasic*>(pmat);
	if (ptp)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;
		
		if (strcmp(szvar, "solid") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = ptp->m_pSolid;
			return FindElasticMaterialParameter(szvar2, index, pme);
		}
		else if (strcmp(szvar, "permeability") == 0)
		{
			// search the nested material parameter list
			pmat = ptp->m_pPerm;
			FEParam* pp = ptp->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
		else if (strcmp(szvar, "osmotic_coefficient") == 0)
		{
			// search the nested material parameter list
			pmat = ptp->m_pOsmC;
			FEParam* pp = ptp->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
		else if (strcmp(szvar, "solute") == 0)
		{
			char* ch = strchr((char*)szvar2, '.');
			if (ch == 0) return 0;
			*ch = 0;
			const char* szvar3 = ch+1;
			
			for (int i=0; i<(int)ptp->m_pSolute.size(); ++i) {
				if (strcmp(szvar2, ptp->m_pSolute[i]->GetName()) == 0)
				{
					FESolute* pms = ptp->m_pSolute[i];
					char* ch = strchr((char*)szvar3, '.');
					if (ch == 0) return 0;
					*ch = 0;
					const char* szvar4 = ch+1;
					if (strcmp(szvar3, "diffusivity") == 0)
					{
						// search the nested material parameter list
						FESoluteDiffusivity* pmd = pms->m_pDiff;
						FEParam* pp = pmd->GetParameter(szvar4);
						if (pp) return pp->pvalue<double>(index);
						else return 0;
					}
					else if (strcmp(szvar3, "solubility") == 0)
					{
						// search the nested material parameter list
						FESoluteSolubility* pmd = pms->m_pSolub;
						FEParam* pp = pmd->GetParameter(szvar4);
						if (pp) return pp->pvalue<double>(index);
						else return 0;
					}
				}
			}
			// no match found
			return 0;
		}
	}
	
	// if this material is a multiphasic material, check solid, permeability,
	// osmotic_coefficient, and solute materials
	FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
	if (pmp)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;
		
		if (strcmp(szvar, "solid") == 0)
		{
			// search the nested material parameter list
			FEElasticMaterial* pme = pmp->m_pSolid;
			return FindElasticMaterialParameter(szvar2, index, pme);
		}
		else if (strcmp(szvar, "permeability") == 0)
		{
			// search the nested material parameter list
			pmat = pmp->m_pPerm;
			FEParam* pp = pmat->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
		else if (strcmp(szvar, "osmotic_coefficient") == 0)
		{
			// search the nested material parameter list
			pmat = pmp->m_pOsmC;
			FEParam* pp = pmat->GetParameter(szvar2);
			if (pp) return pp->pvalue<double>(index);
			else return 0;
		}
		else if (strcmp(szvar, "solute") == 0)
		{
			char* ch = strchr((char*)szvar2, '.');
			if (ch == 0) return 0;
			*ch = 0;
			const char* szvar3 = ch+1;
			
			for (int i=0; i<(int)pmp->m_pSolute.size(); ++i) {
				if (strcmp(szvar2, pmp->m_pSolute[i]->GetName()) == 0)
				{
					FESolute* pms = pmp->m_pSolute[i];
					char* ch = strchr((char*)szvar3, '.');
					if (ch == 0) return 0;
					*ch = 0;
					const char* szvar4 = ch+1;
					if (strcmp(szvar3, "diffusivity") == 0)
					{
						// search the nested material parameter list
						FESoluteDiffusivity* pmd = pms->m_pDiff;
						FEParam* pp = pmd->GetParameter(szvar4);
						if (pp) return pp->pvalue<double>(index);
						else return 0;
					}
					else if (strcmp(szvar3, "solubility") == 0)
					{
						// search the nested material parameter list
						FESoluteSolubility* pmd = pms->m_pSolub;
						FEParam* pp = pmd->GetParameter(szvar4);
						if (pp) return pp->pvalue<double>(index);
						else return 0;
					}
				}
			}
			// no match found
			return 0;
		}
	}

	// the rigid bodies are dealt with differently
	int nrb = m_Obj.size();
	for (i=0; i<nrb; ++i)
	{
		FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[i]);

		if (rb.m_mat == nmat)
		{
			if (strcmp(szvar, "Fx") == 0) return &rb.m_Fr.x;
			if (strcmp(szvar, "Fy") == 0) return &rb.m_Fr.y;
			if (strcmp(szvar, "Fz") == 0) return &rb.m_Fr.z;
			if (strcmp(szvar, "Mx") == 0) return &rb.m_Mr.x;
			if (strcmp(szvar, "My") == 0) return &rb.m_Mr.y;
			if (strcmp(szvar, "Mz") == 0) return &rb.m_Mr.z;
		}
	}

	// oh, oh, we didn't find it
	return 0;
}

//-----------------------------------------------------------------------------
//! Evaluate a parameter list
void FEBioModel::EvaluateParameterList(FEParameterList &pl)
{
	list<FEParam>::iterator pi = pl.first();
	for (int j=0; j<pl.Parameters(); ++j, ++pi)
	{
		if (pi->m_nlc >= 0)
		{
			double v = GetLoadCurve(pi->m_nlc)->Value();
			switch (pi->m_itype)
			{
			case FE_PARAM_INT   : pi->value<int>() = (int) v; break;
			case FE_PARAM_DOUBLE: pi->value<double>() = pi->m_scl*v; break;
			case FE_PARAM_BOOL  : pi->value<bool>() = (v > 0? true : false); break;
			default: 
				assert(false);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function evaluates material parameter lists. Since some of the materials
//! can have other materials as sub-componenents, we need to set up a recursive
//! call to evaluate the parameter lists of the sub-materials.
void FEBioModel::EvaluateMaterialParameters(FEMaterial* pm)
{
	// evaluate the materials' parameter list
	EvaluateParameterList(pm->GetParameterList());

	// evaluate fiber material properties for trans-iso materials
	FETransverselyIsotropic* pti = dynamic_cast<FETransverselyIsotropic*>(pm);
	if (pti)
	{
		EvaluateMaterialParameters(&pti->m_fib);
	}

	// for elastic and uncoupled elastic mixtures, as well as biphasic
	// and biphasic-solute materials we also need to evaluate
	// the sub-materials
	FEElasticMixture* pem = dynamic_cast<FEElasticMixture*>(pm);
	if (pem)
	{
		int NMAT = pem->Materials();
		for (int i=0; i<NMAT; ++i)
			EvaluateMaterialParameters(pem->GetMaterial(i));
	}
	
	FEUncoupledElasticMixture* pum = dynamic_cast<FEUncoupledElasticMixture*>(pm);
	if (pum)
	{
		for (int i=0; i < (int) pum->m_pMat.size(); ++i)
			EvaluateMaterialParameters(pum->m_pMat[i]);
	}
	
	FEElasticMultigeneration* pmg = dynamic_cast<FEElasticMultigeneration*>(pm);
	if (pmg)
	{
		for (int i=0; i < (int) pmg->m_pMat.size(); ++i)
			EvaluateMaterialParameters(pmg->m_pMat[i]);
	}

	FEBiphasic* pb = dynamic_cast<FEBiphasic*>(pm);
	if (pb)
	{
		EvaluateMaterialParameters(pb->m_pSolid);
		EvaluateMaterialParameters(pb->m_pPerm);
	}

	FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*>(pm);
	if (pbs)
	{
		EvaluateMaterialParameters(pbs->m_pSolid );
		EvaluateMaterialParameters(pbs->m_pPerm  );
		EvaluateMaterialParameters(pbs->m_pOsmC  );
		EvaluateMaterialParameters(pbs->m_pSolute);
	}

	FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pm);
	if (pmp)
	{
		EvaluateMaterialParameters(pmp->m_pSolid);
		EvaluateMaterialParameters(pmp->m_pPerm );
		EvaluateMaterialParameters(pmp->m_pOsmC );
		for (int i=0; i<(int)pmp->m_pSolute.size(); ++i)
			EvaluateMaterialParameters(pmp->m_pSolute[i]);
	}

	FESolute* ps = dynamic_cast<FESolute*>(pm);
	if (ps)
	{
		EvaluateMaterialParameters(ps->m_pDiff );
		EvaluateMaterialParameters(ps->m_pSolub);
		if (ps->m_pSupp) EvaluateMaterialParameters(ps->m_pSupp );
	}
}

//=============================================================================
//    I N P U T
//=============================================================================

//-----------------------------------------------------------------------------
//! This routine reads in an input file and performs some initialization stuff.
//! The rest of the initialization is done in Init

bool FEBioModel::Input(const char* szfile)
{
	// create file reader
	FEFEBioImport fim;

	// Load the file
	if (fim.Load(*this, szfile) == false)
	{
		char szerr[256];
		fim.GetErrorMessage(szerr);
		fprintf(stderr, szerr);

		return false;
	}

	// see if user redefined output filenames
	if (fim.m_szdmp[0]) SetDumpFilename(fim.m_szdmp);
	if (fim.m_szlog[0]) SetLogFilename (fim.m_szlog);
	if (fim.m_szplt[0]) SetPlotFilename(fim.m_szplt);

	// we're done reading
	return true;
}

//=============================================================================
//    O U T P U T
//=============================================================================

//-----------------------------------------------------------------------------
//! Export state to plot file.
void FEBioModel::Write()
{
	m_plot->Write(*this);
}

//-----------------------------------------------------------------------------
//! Write user data to the logfile
void FEBioModel::WriteData()
{
	m_Data.Write();
}

//-----------------------------------------------------------------------------
//! Dump state to archive for restarts
void FEBioModel::DumpData()
{
	DumpFile ar(this);
	if (ar.Create(m_szdump) == false)
	{
		clog.printf("WARNING: Failed creating restart point.\n");
	}
	else 
	{
		Serialize(ar);
		clog.printf("\nRestart point created. Archive name is %s\n", m_szdump);
	}
}

//=============================================================================
//    R E S T A R T
//=============================================================================

//-----------------------------------------------------------------------------
//!  Reads or writes the current state to/from a binary file
//!  This is used to restart the solution from a saved position
//!  or to create a restart point.
//!  A version number is written to file to make sure the same
//!  format is used for reading and writing.
//! \param[in] ar the archive to which the data is serialized
//! \sa DumpFile
bool FEBioModel::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		// --- version number ---
		ar << (int) RSTRTVERSION;
	}
	else
	{
		// --- version ---
		int nversion;
		ar >> nversion;

		// make sure it is the right version
		if (nversion != RSTRTVERSION) return false;
	}

	// --- Load Data ---
	SerializeLoadData(ar);

	// --- Serialize global constants ---
	SerializeConstants(ar);

	// --- Material Data ---
	SerializeMaterials(ar);

	// --- Geometry Data ---
	SerializeGeometry(ar);

	// --- Contact Data ---
	SerializeContactData(ar);

	// --- Boundary Condition Data ---
	SerializeBoundaryData(ar);

	// --- Analysis data ---
	SerializeAnalysisData(ar);

	// --- Save IO Data
	SerializeIOData(ar);

	return true;
}


//-----------------------------------------------------------------------------
//! Serialize load curves
void FEBioModel::SerializeLoadData(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		// save curve data
		ar << LoadCurves();
		for (int i=0; i<LoadCurves(); ++i) GetLoadCurve(i)->Serialize(ar);
	}
	else
	{
		// loadcurve data
		int nlc;
		ar >> nlc;
		m_LC.clear();
		for (int i=0; i<nlc; ++i)
		{
			FELoadCurve* plc = new FELoadCurve();
			plc->Serialize(ar);
			AddLoadCurve(plc);
		}
	}
}

//-----------------------------------------------------------------------------
//! Serialize global constants
void FEBioModel::SerializeConstants(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		int NC = (int) m_Const.size();
		ar << NC;
		if (NC > 0)
		{
			char sz[256] = {0};
			map<string, double>::iterator it;
			for (it = m_Const.begin(); it != m_Const.end(); ++it)
			{
				strcpy(sz, it->first.c_str());
				ar << sz;
				ar << it->second;
			}
		}
	}
	else
	{
		char sz[256] = {0};
		double v;
		int NC;
		ar >> NC;
		m_Const.clear();
		for (int i=0; i<NC; ++i)
		{
			ar >> sz >> v;
			SetGlobalConstant(string(sz), v);
		}
	}
}

//-----------------------------------------------------------------------------
//! Serialize analysis data
void FEBioModel::SerializeAnalysisData(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		// analysis steps
		ar << (int) m_Step.size();
		for (int i=0; i<(int) m_Step.size(); ++i) 
		{
			int ntype = m_Step[i]->GetType();
			ar << ntype;
			m_Step[i]->Serialize(ar);
		}

		ar << m_nStep;
		ar << m_ftime << m_ftime0;
		ar << m_nplane_strain;

		// direct solver data
		ar << m_nsolver;
		ar << m_bwopt;

		// body loads
		ar << (int) m_BL.size();
		for (int i=0; i<(int) m_BL.size(); ++i)
		{
			FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_BL[i]);
			int ntype = -1;
			if (dynamic_cast<FEConstBodyForce*      >(pbf)) ntype = FE_CONST_BODY_FORCE;
			if (dynamic_cast<FENonConstBodyForce*   >(pbf)) ntype = FE_NONCONST_BODY_FORCE;
			if (dynamic_cast<FECentrifugalBodyForce*>(pbf)) ntype = FE_CENTRIFUGAL_BODY_FORCE;
			assert(ntype);
			ar << ntype;
			pbf->Serialize(ar);
		}
	}
	else
	{
		m_Step.clear();
		FEModel* pfem = ar.GetFEModel();

		// analysis steps
		int nsteps, ntype;
		ar >> nsteps;
		for (int i=0; i<nsteps; ++i)
		{
			ar >> ntype;
			FEAnalysisStep* pstep = 0;
			switch (ntype)
			{
			case FE_SOLID         : pstep = new FESolidAnalysis         (*this); break;
			case FE_EXPLICIT_SOLID: pstep = new FEExplicitSolidAnalysis         (*this); break;
			case FE_BIPHASIC      : pstep = new FEBiphasicAnalysis      (*this); break;
			case FE_HEAT          : pstep = new FEHeatTransferAnalysis  (*this); break;
			case FE_POROSOLUTE    : pstep = new FEBiphasicSoluteAnalysis(*this); break;
			case FE_LINEAR_SOLID  : pstep = new FELinearSolidAnalysis   (*this); break;
			case FE_HEAT_SOLID    : pstep = new FEThermoElasticAnalysis (*this); break;
			default:
				assert(false);
			}
			pstep->Serialize(ar);
			m_Step.push_back(pstep);
		}
		ar >> m_nStep;
		ar >> m_ftime >> m_ftime0;
		ar >> m_nplane_strain;

		// direct solver data
		ar >> m_nsolver;
		ar >> m_bwopt;

		// body loads
		int nbl;
		ar >> nbl;
		m_BL.clear();
		for (int i=0; i<nbl; ++i)
		{
			int ntype = -1;
			ar >> ntype;
			FEBodyForce* pbl = 0;
			switch (ntype)
			{
			case FE_CONST_BODY_FORCE      : pbl = new FEConstBodyForce      (pfem); break;
			case FE_NONCONST_BODY_FORCE   : pbl = new FENonConstBodyForce   (pfem); break;
			case FE_CENTRIFUGAL_BODY_FORCE: pbl = new FECentrifugalBodyForce(pfem); break;
			default:
				assert(false);
			}
			assert(pbl);
			pbl->Serialize(ar);
			m_BL.push_back(pbl);
		}

		// set the correct step
		m_pStep = m_Step[m_nStep];
	}
}

//-----------------------------------------------------------------------------
//! serialize material data
void FEBioModel::SerializeMaterials(DumpFile& ar)
{
	FEBioKernel& febio = FEBioKernel::GetInstance();

	if (ar.IsSaving())
	{
		// store the nr of materials
		ar << Materials();

		// store the materials
		for (int i=0; i<Materials(); ++i)
		{
			FEMaterial* pmat = GetMaterial(i);

			// store the type string
			ar << febio.GetTypeStr<FEMaterial>(pmat);

			// store the name
			ar << pmat->GetName();

			// store material parameters
			pmat->Serialize(ar);
		}
	}
	else
	{
		// read the number of materials
		int nmat;
		ar >> nmat;

		// read the material data
		char szmat[256] = {0}, szvar[256] = {0};
		for (int i=0; i<nmat; ++i)
		{
			// read the type string
			ar >> szmat;

			// create a material
			FEMaterial* pmat = febio.Create<FEMaterial>(szmat, this);
			assert(pmat);

			// read the name
			ar >> szmat;
			pmat->SetName(szmat);

			// read all parameters
			pmat->Serialize(ar);

			// call init in case this function initializes other data
			pmat->Init();

			// Add material and parameter list to FEM
			AddMaterial(pmat);
		}
	}
}

//-----------------------------------------------------------------------------
// TODO: serialize nonlinear constraints
void FEBioModel::SerializeGeometry(DumpFile &ar)
{
	// serialize the mesh first 
	SerializeMesh(ar);

	// serialize the other geometry data
	if (ar.IsSaving())
	{
		int i;

		// FE objects
		int nrb = m_Obj.size();
		ar << nrb;
		for (i=0; i<nrb; ++i) m_Obj[i]->Serialize(ar);
	}
	else
	{
		int i;

		// rigid bodies
		int nrb;
		ar >> nrb;
		m_Obj.clear();
		for (i=0; i<nrb; ++i)
		{
			FERigidBody* prb = new FERigidBody(this);
			prb->Serialize(ar);
			m_Obj.push_back(prb);
		}
	}
}

//-----------------------------------------------------------------------------
//! This function is used by the restart feature and reads or writes
//! the mesh data to or from the binary archive
//! \param[in] ar the archive to which the data is serialized
//! \sa FEM::Serialize()
//! \sa DumpFile
//! \todo serialize nodesets

void FEBioModel::SerializeMesh(DumpFile& ar)
{
	FEMesh& m = m_mesh;

	if (ar.IsSaving())
	{
		int i;

		// write nodal data
		int nn = m.Nodes();
		ar << nn;
		for (i=0; i<nn; ++i) ar.write(&m.Node(i), sizeof(FENode), 1);

		// write domain data
		int ND = m.Domains();
		ar << ND;
		for (i=0; i<ND; ++i)
		{
			FEDomain& d = m.Domain(i);
			int ntype = d.Type();
			int ne = d.Elements();
			int nmat = (int) d.GetMaterial()->GetID() - 1; 
			ar << nmat;
			ar << ntype << ne;
			d.Serialize(ar);
		}
	}
	else
	{
		int i;

		// read nodal data
		int nn;
		ar >> nn;
		m.CreateNodes(nn);
		for (i=0; i<nn; ++i) ar.read(&m.Node(i), sizeof(FENode), 1);

		// read domain data
		int ND;
		ar >> ND;
		for (i=0; i<ND; ++i)
		{
			int nmat;
			ar >> nmat;
			FEMaterial* pm = GetMaterial(nmat);
			assert(pm);

			int ntype, ne;
			ar >> ntype >> ne;
			FEDomain* pd = 0;
			switch (ntype)
			{
			case FE_SOLID_DOMAIN          : pd = new FEElasticSolidDomain      (&m, pm); break;
			case FE_SHELL_DOMAIN          : pd = new FEElasticShellDomain      (&m, pm); break;
			case FE_TRUSS_DOMAIN          : pd = new FEElasticTrussDomain      (&m, pm); break;
			case FE_RIGID_SOLID_DOMAIN    : pd = new FERigidSolidDomain        (&m, pm); break;
			case FE_RIGID_SHELL_DOMAIN    : pd = new FERigidShellDomain        (&m, pm); break;
			case FE_UDGHEX_DOMAIN         : pd = new FEUDGHexDomain            (&m, pm); break;
			case FE_HEAT_SOLID_DOMAIN     : pd = new FEHeatSolidDomain         (&m, pm); break;
			case FE_DISCRETE_DOMAIN       : pd = new FEDiscreteSpringDomain    (&m, pm); break;
			case FE_3F_SOLID_DOMAIN       : pd = new FE3FieldElasticSolidDomain(&m, pm); break;
			case FE_BIPHASIC_DOMAIN       : pd = new FEBiphasicSolidDomain          (&m, pm); break;
			case FE_BIPHASIC_SOLUTE_DOMAIN: pd = new FEBiphasicSoluteDomain    (&m, pm); break;
			case FE_UT4_DOMAIN            : pd = new FEUT4Domain               (&m, pm); break;
			default: assert(false);
			}

			assert(pd);
			pd->create(ne);
			pd->Serialize(ar);

			m.AddDomain(pd);
		}

		m.UpdateBox();
	}
}

//-----------------------------------------------------------------------------
//! serialize contact data
void FEBioModel::SerializeContactData(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << ContactInterfaces();
		for (int i=0; i<ContactInterfaces(); ++i)
		{
			ar << m_CI[i]->Type();
			m_CI[i]->Serialize(ar);
		}
	}
	else
	{
		int numci, ntype;
		ar >> numci;
		for (int i=0; i<numci; ++i)
		{
			FEContactInterface* ps;

			// get the interface type
			ar >> ntype;

			// create a new interface
			switch (ntype)
			{
			case FE_CONTACT_SLIDING      : ps = new FESlidingInterface     (this); break;
			case FE_FACET2FACET_SLIDING  : ps = new FEFacet2FacetSliding   (this); break;
			case FE_CONTACT_TIED         : ps = new FETiedInterface        (this); break;
			case FE_CONTACT_RIGIDWALL    : ps = new FERigidWallInterface   (this); break;
			case FE_CONTACT_SLIDING2     : ps = new FESlidingInterface2    (this); break;
			case FE_PERIODIC_BOUNDARY    : ps = new FEPeriodicBoundary     (this); break;
			case FE_SURFACE_CONSTRAINT   : ps = new FESurfaceConstraint    (this); break;
			case FE_CONTACT_SLIDING3     : ps = new FESlidingInterface3    (this); break;
			case FE_CONTACT_TIED_BIPHASIC: ps = new FETiedBiphasicInterface(this); break;
			case FE_CONTACT_SLIDINGBW    : ps = new FESlidingInterfaceBW   (this); break;
			default:
				assert(false);
			}
				
			// serialize interface data from archive
			ps->Serialize(ar);

			// add interface to list
			m_CI.push_back(ps);
		}	
	}
}

//-----------------------------------------------------------------------------
// TODO: do we need to store the m_bActive flag of the boundary conditions?
//
void FEBioModel::SerializeBoundaryData(DumpFile& ar)
{
	int i, n;

	if (ar.IsSaving())
	{
		// displacements
		ar << (int) m_DC.size();
		for (i=0; i<(int) m_DC.size(); ++i) 
		{
			FEPrescribedBC& dc = *m_DC[i];
			ar << dc.GetID() << dc.IsActive();
			ar << dc.bc << dc.lc << dc.node << dc.s;
		}

		// nodal loads
		ar << (int) m_FC.size();
		for (i=0; i<(int) m_FC.size(); ++i)
		{
			FENodalForce& fc = *m_FC[i];
			ar << fc.GetID() << fc.IsActive();
			ar << fc.bc << fc.lc << fc.node << fc.s;
		}

		// surface loads
		ar << (int) m_SL.size();
		for (i=0; i<(int) m_SL.size(); ++i)
		{
			FESurfaceLoad* psl = m_SL[i];

			// get the surface
			FESurface& s = psl->Surface();
			s.Serialize(ar);

			// save the load data
			int ntype = -1;
			if (dynamic_cast<FEPressureLoad*      >(psl)) ntype = FE_PRESSURE_LOAD;
			if (dynamic_cast<FETractionLoad*      >(psl)) ntype = FE_TRACTION_LOAD;
			if (dynamic_cast<FEFluidFlux*         >(psl)) ntype = FE_FLUID_FLUX;
			if (dynamic_cast<FEPoroNormalTraction*>(psl)) ntype = FE_PORO_TRACTION;
			if (dynamic_cast<FESoluteFlux*        >(psl)) ntype = FE_SOLUTE_FLUX;
			if (dynamic_cast<FEHeatFlux*          >(psl)) ntype = FE_HEAT_FLUX;
			if (dynamic_cast<FEConvectiveHeatFlux*>(psl)) ntype = FE_CONV_HEAT_FLUX;
			assert(ntype != -1);
			ar << ntype;
			ar << psl->GetID() << psl->IsActive();
			psl->Serialize(ar);
		}

		// rigid body displacements
		ar << m_RDC.size();
		for (i=0; i<(int) m_RDC.size(); ++i)
		{
			FERigidBodyDisplacement& dc = *m_RDC[i];
			ar << dc.GetID() << dc.IsActive();
			ar << dc.bc << dc.id << dc.lc << dc.sf;
		}

		// rigid body forces
		ar << m_RFC.size();
		for (i=0; i<(int) m_RFC.size(); ++i)
		{
			FERigidBodyForce& fc = *m_RFC[i];
			ar << fc.GetID() << fc.IsActive();
			ar << fc.bc << fc.id << fc.lc << fc.sf;
		}

		// rigid nodes
		ar << m_RN.size();
		for (i=0; i<(int) m_RN.size(); ++i)
		{
			FERigidNode& rn = *m_RN[i];
			ar << rn.GetID() << rn.IsActive();
			ar << rn.nid << rn.rid;
		}

		// linear constraints
		ar << (int) m_LinC.size();
		list<FELinearConstraint>::iterator it = m_LinC.begin();
		for (i=0; i<(int) m_LinC.size(); ++i, ++it) it->Serialize(ar);

		ar << m_LCT;

		// aug lag linear constraints
/*		n = (int) m_LCSet.size();
		ar << n;
		if (m_LCSet.empty() == false)
		{
			for (i=0; i<n; ++i) m_LCSet[i]->Serialize(ar);
		}
*/
		n = m_NLC.size();
		ar << n;
		if (n) 
		{
			for (i=0; i<n; ++i) 
			{
//				ar << m_NLC[i]->Type();
//				m_NLC[i]->Serialize(ar);
			}
		}
	}
	else
	{
		int n;
		int nid; bool bactive;

		// displacements
		ar >> n;
		m_DC.clear();
		for (i=0; i<n; ++i) 
		{
			FEPrescribedBC* pdc = new FEPrescribedBC;
			ar >> nid >> bactive;
			ar >> pdc->bc >> pdc->lc >> pdc->node >> pdc->s >> pdc->br >> pdc->r; // GAA
			pdc->SetID(nid);
			if (bactive) pdc->Activate(); else pdc->Deactivate();
			m_DC.push_back(pdc);
		}
		
		// nodal loads
		ar >> n;
		m_FC.clear();
		for (i=0; i<n; ++i)
		{
			FENodalForce* pfc = new FENodalForce;
			ar >> nid >> bactive;
			ar >> pfc->bc >> pfc->lc >> pfc->node >> pfc->s;
			pfc->SetID(nid);
			if (bactive) pfc->Activate(); else pfc->Deactivate();
			m_FC.push_back(pfc);
		}

		// surface loads
		ar >> n;
		m_SL.clear();
		for (i=0; i<n; ++i)
		{
			// create a new surface
			FESurface* psurf = new FESurface(&m_mesh);
			psurf->Serialize(ar);

			// read load data
			int ntype;
			ar >> ntype;
			FESurfaceLoad* ps = 0;
			switch (ntype)
			{
			case FE_PRESSURE_LOAD : ps = new FEPressureLoad      (psurf); break;
			case FE_TRACTION_LOAD : ps = new FETractionLoad      (psurf); break;
			case FE_FLUID_FLUX    : ps = new FEFluidFlux         (psurf); break;
			case FE_PORO_TRACTION : ps = new FEPoroNormalTraction(psurf); break;
			case FE_SOLUTE_FLUX   : ps = new FESoluteFlux        (psurf); break;
			case FE_HEAT_FLUX     : ps = new FEHeatFlux          (psurf); break;
			case FE_CONV_HEAT_FLUX: ps = new FEConvectiveHeatFlux(psurf); break;
			default:
				assert(false);
			}
			assert(ps);

			ar >> nid >> bactive;
			ps->SetID(nid);
			if (bactive) ps->Activate(); else ps->Deactivate();

			ps->Serialize(ar);
			m_SL.push_back(ps);
		}

		// rigid body displacements
		ar >> n;
		m_RDC.clear();
		for (i=0; i<n; ++i)
		{
			FERigidBodyDisplacement* pdc = new FERigidBodyDisplacement;
			ar >> nid >> bactive;
			ar >> pdc->bc >> pdc->id >> pdc->lc >> pdc->sf;
			pdc->SetID(nid);
			if (bactive) pdc->Activate(); else pdc->Deactivate();
			m_RDC.push_back(pdc);
		}

		// rigid body forces
		ar >> n;
		m_RFC.clear();
		for (i=0; i<n; ++i)
		{
			FERigidBodyForce* pfc = new FERigidBodyForce;
			ar >> nid >> bactive;
			ar >> pfc->bc >> pfc->id >> pfc->lc >> pfc->sf;
			pfc->SetID(nid);
			if (bactive) pfc->Activate(); else pfc->Deactivate();
			m_RFC.push_back(pfc);
		}

		// rigid nodes
		ar >> n;
		m_RN.clear();
		for (i=0; i<n; ++i)
		{
			FERigidNode* prn = new FERigidNode;
			ar >> nid >> bactive;
			ar >> prn->nid >> prn->rid;
			prn->SetID(nid);
			if (bactive) prn->Activate(); else prn->Deactivate();
			m_RN.push_back(prn);
		}

		// linear constraints
		ar >> n;
		FELinearConstraint LC;
		for (i=0; i<n; ++i)
		{
			LC.Serialize(ar);
			m_LinC.push_back(LC);
		}

		ar >> m_LCT;

		// reset the pointer table
		int nlin = m_LinC.size();
		m_LCA.resize(nlin);
		list<FELinearConstraint>::iterator ic = m_LinC.begin();
		for (i=0; i<nlin; ++i, ++ic) m_LCA[i] = &(*ic);

		// aug lag linear constraints
		ar >> n;
//		int ntype;
		m_NLC.clear();
		for (i=0; i<n; ++i)
		{
/*			ar >> ntype;
			FENLConstraint* plc = 0;
			switch (ntype)
			{
			case FE_POINT_CONSTRAINT : plc = new FEPointConstraint    (this); break;
			case FE_LINEAR_CONSTRAINT: plc = new FELinearConstraintSet(this); break;
			default:
				assert(false);
			}
			assert(plc);
			plc->Serialize(ar);
			m_NLC.push_back(plc);
*/		}
	}
}

//-----------------------------------------------------------------------------
//! Serialization of FEBioModel data
void FEBioModel::SerializeIOData(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		// file names
		ar << m_szfile << m_szplot << m_szlog << m_szdump;
		ar << m_sztitle;

		// plot file
		int npltfmt = 0;
		if (dynamic_cast<LSDYNAPlotFile*>(m_plot)) npltfmt = 1;
		else if (dynamic_cast<FEBioPlotFile*>(m_plot)) npltfmt = 2;
		assert(npltfmt != 0);
		ar << npltfmt;

		if (npltfmt == 1)
		{
			LSDYNAPlotFile* plt = dynamic_cast<LSDYNAPlotFile*>(m_plot);

			int* n = plt->m_nfield;
			ar << n[0] << n[1] << n[2] << n[3] << n[4];
		}

		// data records
		SerializeDataStore(ar);
	}
	else
	{
		// file names
		ar >> m_szfile >> m_szplot >> m_szlog >> m_szdump;
		ar >> m_sztitle;

		// don't forget to call store the input file name so
		// that m_szfile_title gets initialized
		SetInputFilename(m_szfile);

		// get the plot file format
		int npltfmt = 0;
		ar >> npltfmt;
		assert(m_plot == 0);

		switch (npltfmt)
		{
		case 1:
			{
				// Open the plot file for appending
				// TODO: We need a better way to create a plotfile
				//		 what if the user created a different output format?
				LSDYNAPlotFile* plt = new LSDYNAPlotFile;
				m_plot = plt;
				if (m_plot->Append(*this, m_szplot) == false)
				{
					printf("FATAL ERROR: Failed reopening plot database %s\n", m_szplot);
					throw "FATAL ERROR";
				}

				// plot file
				int* n = plt->m_nfield;
				ar >> n[0] >> n[1] >> n[2] >> n[3] >> n[4];
			}
			break;
		case 2:
			{
				m_plot = new FEBioPlotFile(*this);
				if (m_plot->Append(*this, m_szplot) == false)
				{
					printf("FATAL ERROR: Failed reopening plot database %s\n", m_szplot);
					throw "FATAL ERROR";
				}
			}
			break;
		};

		// data records
		SerializeDataStore(ar);
	}
}

//-----------------------------------------------------------------------------
void FEBioModel::SerializeDataStore(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		int N = m_Data.Size();
		ar << N;
		for (int i=0; i<N; ++i)
		{
			DataRecord* pd = m_Data.GetDataRecord(i);

			int ntype = -1;
			if (dynamic_cast<NodeDataRecord*>(pd)) ntype = FE_DATA_NODE;
			if (dynamic_cast<ElementDataRecord*>(pd)) ntype = FE_DATA_ELEM;
			if (dynamic_cast<RigidBodyDataRecord*>(pd)) ntype = FE_DATA_RB;
			assert(ntype != -1);
			ar << ntype;
			pd->Serialize(ar);
		}
	}
	else
	{
		int N;
		m_Data.Clear();
		ar >> N;
		for (int i=0; i<N; ++i)
		{
			int ntype;
			ar >> ntype;

			DataRecord* pd = 0;
			switch(ntype)
			{
			case FE_DATA_NODE: pd = new NodeDataRecord(this, 0); break;
			case FE_DATA_ELEM: pd = new ElementDataRecord(this, 0); break;
			case FE_DATA_RB  : pd = new RigidBodyDataRecord(this, 0); break;
			}
			assert(pd);
			pd->Serialize(ar);
			m_Data.AddRecord(pd);
		}
	}
}

//=============================================================================
//    I N I T I A L I Z A T I O N
//=============================================================================

//-----------------------------------------------------------------------------
// Forward declarations
void Hello();

//-----------------------------------------------------------------------------
//! This function performs one-time-initialization stuff. All the different 
//! modules are initialized here as well. This routine also performs some
//! data checks

bool FEBioModel::Init()
{
	int i;

	// Open the logfile
	if (!clog.is_valid()) 
	{
		if (clog.open(m_szlog) == false)
		{
			clog.printbox("FATAL ERROR", "Failed creating log file");
			return false;
		}

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) clog.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Logfile::MODE m = clog.SetMode(Logfile::FILE_ONLY);
		Hello();
		clog.SetMode(m);
	}

	// intitialize time
	m_ftime = 0;
	m_ftime0 = 0;

	// check step data
	// TODO: should I let the Steps take care of this instead?
	for (i=0; i<(int) m_Step.size(); ++i)
	{
		FEAnalysis& step = *m_Step[i];
		if ((step.m_ntime <= 0) && (step.m_final_time <= 0.0)) { clog.printf("Invalid number of time steps for analysis step %d", i+1); return false; }
		if ((step.m_ntime >  0) && (step.m_final_time >  0.0)) { clog.printf("You must either set the number of time steps or the final time but not both.\n"); return false; }
		if (step.m_dt0   <= 0) { clog.printf("Invalid time step size for analysis step %d", i+1); return false; }
		if (step.m_bautostep)
		{
//			if (m_pStep->m_dtmin <= 0) return err("Invalid minimum time step size");
//			if (m_pStep->m_dtmax <= 0) return err("Invalid maximum time step size");
		}
	}

	// evaluate all loadcurves at the initial time
	for (i=0; i<LoadCurves(); ++i) m_LC[i]->Evaluate(0);

	// if the analysis is run in plain-strain mode we fix all the z-dofs of all nodes
	if (m_nplane_strain >= 0)
	{
		int bc = m_nplane_strain;
		for (int i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[bc] = -1;
	}

	// find and remove isolated vertices
	int ni = m_mesh.RemoveIsolatedVertices();
	if (ni != 0) 
	{
		if (ni == 1)
			clog.printbox("WARNING", "%d isolated vertex removed.", ni);
		else
			clog.printbox("WARNING", "%d isolated vertices removed.", ni);
	}

	// create and initialize the rigid body data
	if (CreateRigidBodies() == false) return false;

	// initialize poroelastic/biphasic and solute data
	if (InitPoroSolute() == false) return false;

	// initialize random number generator
	srand((unsigned) time(NULL));

	// initialize mesh data
	// note that this must be done AFTER the elements have been assigned material point data !
	// this is because the mesh data is reset
	// TODO: perhaps I should not reset the mesh data during the initialization
	if (InitMesh() == false) return false;

	// initialize material data
	if (InitMaterials() == false) return false;

	// initialize contact data
	if (InitContact() == false) return false;

	// init some other stuff
	for (i=0; i<(int) m_BL.size(); ++i)
	{
		if (m_BL[i]->Init() == false) return false;
	}

	// initialize nonlinear constraints
	// TODO: This is also initialized in the analysis step. Do I need to do this here?
	for (i=0; i<(int) m_NLC.size(); ++i) m_NLC[i]->Init();

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) m_plot = new LSDYNAPlotFile;

		if (m_plot->Open(*this, m_szplot) == false)
		{
			clog.printf("ERROR : Failed creating PLOT database\n");
			return false;
		}
	}

	// Since it is assumed that for the first timestep
	// there are no loads or initial displacements, the case n=0 is skipped.
	// Therefor we can output those results here.
	// TODO: Offcourse we should actually check if this is indeed
	//       the case, otherwise we should also solve for t=0
	if (m_pStep->m_nplot != FE_PLOT_NEVER) m_plot->Write(*this);

	// do the callback
	DoCallback();

	// Alright, all initialization is done, so let's get busy !
	return true;
}

//-----------------------------------------------------------------------------
//! Initializes contact data
// TODO: Contact interfaces have two initialization functions: Init() and
//   Activate(). The Init member is called here and allocates the required memory
//   for each interface. The Activate() member is called during the step initialization
//   which is called later during the solve phase. However, for global interfaces (i.e.
//   interfaces that are active during the entire simulation), the Activate() member is
//   not called in the solve phase. That is why we have to call it here. Global interfaces
//   can be idenfitifed since they are active during the initialization. 
//
//   I am not entirely a fan of this approach but it does solve the problem that contact
//   interface shoulds only do work (e.g. update projection status) when they are active, but
//   have to allocate memory during the initialization fase.
//
bool FEBioModel::InitContact()
{
	// loop over all contact interfaces
	for (int i=0; i<ContactInterfaces(); ++i)
	{
		// get the contact interface
		FEContactInterface& ci = *m_CI[i];

		// initializes contact interface data
		if (ci.Init() == false) return false;

		// If the contact interface is active
		// we have to call the Activate() member. 
		if (ci.IsActive()) ci.Activate();
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize material data
bool FEBioModel::InitMaterials()
{
	int i;

	// initialize material data
	for (i=0; i<Materials(); ++i)
	{
		// get the material
		FEMaterial* pmat = GetMaterial(i);

		// initialize material data
		try
		{
			pmat->Init();
		}
		catch (MaterialError e)
		{
			clog.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			clog.printf("ERROR: %s\n\n", e.Error());
			return false;
		}
		catch (MaterialRangeError e)
		{
			clog.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			clog.printf("ERROR: parameter \"%s\" out of range ", e.m_szvar);
			if (e.m_bl) clog.printf("["); else clog.printf("(");
			clog.printf("%lg, %lg", e.m_vmin, e.m_vmax);
			if (e.m_br) clog.printf("]"); else clog.printf(")");
			clog.printf("\n\n");
			return false;
		}
		catch (...)
		{
			clog.printf("A fatal error occured during material intialization\n\n");
			return false;
		}
	}

	// initialize discrete materials
	try
	{
		for (i=0; i<(int) m_MAT.size(); ++i)
		{
			FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(m_MAT[i]);
			if (pm)
			{
				if (dynamic_cast<FENonLinearSpring*>(pm))
				{
					FENonLinearSpring* ps = dynamic_cast<FENonLinearSpring*>(pm);
					ps->m_plc = GetLoadCurve(ps->m_nlc);
				}
				pm->Init();
			}
		}
	}
	catch (...)
	{
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function creates the rigid bodies by analyzing the rigid materials
//! and the mesh in the model. 
//!
bool FEBioModel::CreateRigidBodies()
{
	int i, j, n, m, nd;
	// count the number of rigid materials
	int nrm = 0;
	for (i=0; i<Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(i));
		if (pm) nrm++;
	}
	
	// make sure there are rigid materials
	if (nrm == 0) return true;

	// First we need to figure out how many rigid bodies there are.
	// This is not the same as rigid materials, since a rigid body
	// may be composed of different rigid materials (similarly to a deformable
	// body that may contain different materials). Although there can
	// only be one deformable mesh, there can be several rigid bodies.

	// The mrb array will contain an index to the rigid body the material
	// is attached too.
	vector<int> mrb(Materials());
	n = 0;
	for (i=0; i<Materials(); ++i)
	{
		if (dynamic_cast<FERigidMaterial*>(GetMaterial(i))) mrb[i] = n++;
		else mrb[i] = -1;
	}

	// Next, we assign to all nodes a rigid node number
	// This number is preliminary since rigid materials can be merged
	for (nd = 0; nd < m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					el.m_nrigid = el.GetMatID();
					for (j=0; j<el.Nodes(); ++j)
					{
						n = el.m_node[j];
						FENode& node = m_mesh.Node(n);
						node.m_rid = el.GetMatID();
					}
				}
				else el.m_nrigid = -1;
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					el.m_nrigid = el.GetMatID();
					for (j=0; j<el.Nodes(); ++j)
					{
						n = el.m_node[j];
						FENode& node = m_mesh.Node(n);
						node.m_rid = el.GetMatID();
					}		
				}
				else el.m_nrigid = -1;
			}
		}
	}

	// now we can merge rigid materials
	// if a rigid element has two nodes that connect to two different
	// rigid materials we need to merge. 
	bool bdone;
	do
	{
		bdone = true;
		for (nd=0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (el.m_nrigid >= 0)
					{
						m = m_mesh.Node(el.m_node[0]).m_rid;
						for (j=1; j<el.Nodes(); ++j)
						{
							n = m_mesh.Node(el.m_node[j]).m_rid;
							if (mrb[n] != mrb[m])
							{
								if (mrb[n]<mrb[m]) mrb[m] = mrb[n]; else mrb[n] = mrb[m];
								bdone = false;
							}
						}
					}
				}
			}

			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (el.m_nrigid >= 0)
					{
						m = m_mesh.Node(el.m_node[0]).m_rid;
						for (j=1; j<el.Nodes(); ++j)
						{
							n = m_mesh.Node(el.m_node[j]).m_rid;
							if (mrb[n] != mrb[m])
							{
								if (mrb[n]<mrb[m]) mrb[m] = mrb[n]; else mrb[n] = mrb[m];
								bdone = false;
							}
						}
					}
				}
			}
		}
	}
	while (!bdone);

	// since we may have lost a rigid body in the merge process
	// we reindex the RB's.
	int nmat = Materials();
	vector<int> mrc; mrc.assign(nmat, -1);
	for (i=0; i<nmat; ++i) if (mrb[i] >= 0) mrc[mrb[i]] = 0;
	int nrb = 0;
	for (i=0; i<nmat; ++i)
	{
		if (mrc[i] == 0) mrc[i] = nrb++;
	}

	for (i=0; i<nmat; ++i) 
	{
		if (mrb[i] >= 0) mrb[i] = mrc[mrb[i]];
	}

	// set rigid body index for materials
	for (i=0; i<Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(i));
		if (pm)	
		{
			pm->m_nRB = mrb[i];
		}
	}

	// assign rigid body index to rigid elements
	for (nd=0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(el.GetMatID()));
				if (pm)
					el.m_nrigid = pm->m_nRB;
				else
					el.m_nrigid = -1;
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(el.GetMatID()));
				if (pm)
					el.m_nrigid = pm->m_nRB;
				else
					el.m_nrigid = -1;
			}
		}
	}

	// assign rigid body index to nodes
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0) node.m_rid = mrb[ node.m_rid ];
	}

	// Ok, we now know how many rigid bodies there are
	// so let's create them
	m_Obj.clear();
	for (i=0; i<nrb; ++i)
	{
		// create a new rigid body
		FERigidBody* prb = new FERigidBody(this);
		prb->m_nID = i;

		// Since a rigid body may contain several rigid materials
		// we find the first material that this body has and use
		// that materials data to set up the rigid body data
		FERigidMaterial* pm = 0;
		for (j=0; j<Materials(); ++j)
		{
			pm = dynamic_cast<FERigidMaterial*> (GetMaterial(j));

			if (pm && (pm->m_nRB == i))	break;
		}
		assert(j<Materials());
		prb->m_mat = j;

		// initialize center of mass
		if (pm->m_com == 1)
		{
			// grab the com from the material
			prb->m_r0 = prb->m_rt = pm->m_rc;
		}
		else
		{
			// calculate the com
			prb->UpdateCOM();
		}

		// add it to the pile
		m_Obj.push_back(prb);
	}

	// set up rigid joints
	if (!m_NLC.empty())
	{
		FERigidMaterial* pm;
		int NC = m_NLC.size();
		for (i=0; i<NC; ++i)
		{
			FENLConstraint* plc = m_NLC[i];
			if (dynamic_cast<FERigidJoint*>(plc))
			{
				FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*plc);
				rj.m_F = vec3d(0,0,0);

				pm = dynamic_cast<FERigidMaterial*> (GetMaterial(rj.m_nRBa));
				if (pm == 0)
				{
					clog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
					return false;
				}
				rj.m_nRBa = pm->m_nRB;

				pm = dynamic_cast<FERigidMaterial*> (GetMaterial(rj.m_nRBb));
				if (pm == 0)
				{
					clog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
					return false;
				}
				rj.m_nRBb = pm->m_nRB;

				FERigidBody& ra = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBa]);
				FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBb]);

				rj.m_qa0 = rj.m_q0 - ra.m_r0;
				rj.m_qb0 = rj.m_q0 - rb.m_r0;
			}
		}
	}

	// overwrite rigid nodes degrees of freedom
	// We do this so that these dofs do not
	// get equation numbers assigned to them. Later we'll assign
	// the rigid dofs equations numbers to these nodes
	for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_bshell = false;
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				if (el.m_nrigid < 0)
				{
					int n = el.Nodes();
					for (j=0; j<n; ++j) m_mesh.Node(el.m_node[j]).m_bshell = true;
				}
			}
		}
	}

	// The following fixes the degrees of freedom for rigid nodes.
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0)
		{
			node.m_ID[DOF_X] = -1;
			node.m_ID[DOF_Y] = -1;
			node.m_ID[DOF_Z] = -1;
			if (node.m_bshell == false)
			{
				node.m_ID[DOF_U] = -1;
				node.m_ID[DOF_V] = -1;
				node.m_ID[DOF_W] = -1;
			}
		}
	}

	// assign correct rigid body ID's to rigid nodes
	for (i=0; i<(int) m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_RN[i];
		rn.rid = mrb[rn.rid];
	}

	// let's find all rigid surface elements
	// a surface element is rigid when it has no free nodes
	for (int is = 0; is < (int) m_SL.size(); ++is)
	{
		FESurfaceLoad* ps = m_SL[is];
		for (i=0; i<ps->Surface().Elements(); ++i)
		{
			FESurfaceElement& el = ps->Surface().Element(i);
			int N = el.Nodes();
			el.m_nrigid = 0;
			for (j=0; j<N; ++j) 
			{
				FENode& node = m_mesh.Node(el.m_node[j]);
				if (node.m_rid < 0) 
				{
					el.m_nrigid = -1;
					break;
				}
			}
		}
	}

	// the rigid body constraints are still associated with the rigid materials
	// so we now associate them with the rigid bodies
	for (i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_RDC[i];
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(DC.id-1));
		DC.id = pm->m_nRB;
	}
	for (i=0; i<(int) m_RFC.size(); ++i)
	{
		FERigidBodyForce& FC = *m_RFC[i];
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(FC.id-1));
		FC.id = pm->m_nRB;
	}

	// set the rigid body parents
	for (i=0; i<(int) m_Obj.size(); ++i)
	{
		FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[i]);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_MAT[rb.m_mat]);
		assert(pm);
		if (pm->m_pmid > -1)
		{
			FERigidMaterial* ppm = dynamic_cast<FERigidMaterial*>(m_MAT[pm->m_pmid-1]);
			assert(ppm);
			FERigidBody& prb = dynamic_cast<FERigidBody&>(*m_Obj[ppm->m_nRB]);
			rb.m_prb = &prb;

			// we also need to open up all the RB's degree of freedoms
			pm->m_bc[0] = 1;
			pm->m_bc[1] = 1;
			pm->m_bc[2] = 1;
			pm->m_bc[3] = 1;
			pm->m_bc[4] = 1;
			pm->m_bc[5] = 1;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Update the mesh data. This function calculates the initial directors
//! for the shell elements.

// NOTE: This function needs to be called after the rigid bodies have been
// initialized

bool FEBioModel::InitMesh()
{
	int i, j, n, m0, m1, m2, nd;
	int* en;
	vec3d a, b, c;

	int ninverted = 0;

	FEMesh& m = m_mesh;

	// zero initial directors for shell nodes
	for (i=0; i<m.Nodes(); ++i) m.Node(i).m_D0 = vec3d(0,0,0);

	for (nd = 0; nd < m.Domains(); ++nd)
	{
		// check all solid elements to see if they are not initially inverted
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);

				int nint = el.GaussPoints();
				for (int n=0; n<nint; ++n)
				{
					double J0 = pbd->detJ0(el, n);
					if (J0 <= 0)
					{
						clog.printf("**************************** E R R O R ****************************\n");
						clog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.m_nID);
						clog.printf("Jacobian = %lg\n", J0);
						clog.printf("Did you use the right node numbering?\n");
						clog.printf("Nodes:");
						for (int l=0; l<el.Nodes(); ++l)
						{
							clog.printf("%d", el.m_node[l]+1);
							if (l+1 != el.Nodes()) clog.printf(","); else clog.printf("\n");
						}
						clog.printf("*******************************************************************\n\n");
						++ninverted;
					}
				}
			}
		}

		// Calculate the shell directors as the local node normals
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			vec3d r0[4];
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);

				n = el.Nodes();
				en = &el.m_node[0];

				// get the nodes
				for (j=0; j<n; ++j) r0[j] = psd->GetMesh()->Node(en[j]).m_r0;

				for (j=0; j<n; ++j)
				{
					m0 = j;
					m1 = (j+1)%n;
					m2 = (j==0? n-1: j-1);

					a = r0[m0];
					b = r0[m1];
					c = r0[m2];

					m.Node(en[m0]).m_D0 += (b-a)^(c-a);
				}
			}
		}
	}

	// make sure we start with unit directors
	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		node.m_D0.unit();
		node.m_Dt = node.m_D0;
	}

	// Check initially inverted shell elements
	for (nd = 0; nd < m.Domains(); ++nd)
	{
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			// check the connectivity of the shells
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				if (!el.IsRigid())
				{
					int nint = el.GaussPoints();
					for (int n=0; n<nint; ++n)
					{
						double J0 = psd->detJ0(el, n);
						if (J0 <= 0)
						{
							clog.printf("**************************** E R R O R ****************************\n");
							clog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.m_nID);
							clog.printf("Jacobian = %lg\n", J0);
							clog.printf("Did you use the right node numbering?\n");
							clog.printf("Nodes:");
							for (int l=0; l<el.Nodes(); ++l)
							{
								clog.printf("%d", el.m_node[l]+1);
								if (l+1 != el.Nodes()) clog.printf(","); else clog.printf("\n");
							}
							clog.printf("*******************************************************************\n\n");
							++ninverted;
						}
					}
				}
			}
		}
	}

	// report number of inverted elements
	if (ninverted != 0)
	{
		clog.printf("**************************** E R R O R ****************************\n");
		clog.printf(" FEBio found %d initially inverted elements.\n", ninverted);
		clog.printf(" Run will be aborted.\n");
		clog.printf("*******************************************************************\n\n");
		return false;
	}

	// next if a node does not belong to a shell
	// we turn of the rotational degrees of freedom
	// To do this, we first tag all shell nodes
	vector<int> tag(m.Nodes());
	zero(tag);
	for (nd = 0; nd < m.Domains(); ++nd)
	{
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				n = el.Nodes();
				en = &el.m_node[0];
				for (j=0; j<n; ++j) tag[en[j]] = 1;
			}
		}
	}

	// fix rotational degrees of freedom of tagged nodes
	for (i=0; i<m.Nodes(); ++i) 
	{
		FENode& node = m.Node(i);
		if (tag[i] == 0)
		{
			node.m_ID[DOF_U] = -1;
			node.m_ID[DOF_V] = -1;
			node.m_ID[DOF_W] = -1;
		}
	}

	// At this point, the node ID still contains the boundary conditions
	// so we copy that into the m_BC array
	// TODO: perhaps we should put the BC's initially in BC instead of ID.
	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		for (j=0; j<MAX_NDOFS; ++j)	node.m_BC[j] = node.m_ID[j];
	}

	// reset data
	m.Reset();

	// intialize domains
	for (i=0; i<m_mesh.Domains(); ++i) m_mesh.Domain(i).Initialize(*this);

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize solute-poroelastic data.
//! Find all nodes that are not part of a poro-solute domain and fix the 
//! pressure and concentration DOFS. 
//! \todo This function should probably move to the FEAnalysisStep class.
bool FEBioModel::InitPoroSolute()
{
	int i, j, k, nd;

	// make sure this is the poro-solute module
	int nstep = m_pStep->GetType();
	bool bporo = ((nstep == FE_BIPHASIC) || (nstep == FE_POROSOLUTE));
	bool bsolu = (nstep == FE_POROSOLUTE);
	
	if (!bporo)
	{
		// if there is no poroelasticity
		// we set all pressure degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[DOF_P] = -1;
		
		// also remove prescribed pressures
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if (bc == DOF_P) bc = -1;
		}
	}
	
	if (!bsolu)
	{
		// if there is no solute
		// we set all concentration degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (j=0; j<MAX_CDOFS; ++j) {
			k = DOF_C+j;
			for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[k] = -1;
		}
		
		// also remove prescribed concentrations
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if ((bc >= DOF_C) && (bc < MAX_NDOFS)) bc = -1;
		}
	}
	
	if (!bporo)
		// let's go back
		return true;
	
	// see if we are using the symmetric version or not
/*	if (m_pStep->m_bsym_poro == false) 
	{
		m_pStep->m_psolver->m_bsymm = false;
		
		// make sure we are using full-Newton
//		if (m_pStep->m_psolver->m_bfgs.m_maxups != 0)
//		{
//			m_pStep->m_psolver->m_bfgs.m_maxups = 0;
//			clog.printbox("WARNING", "The non-symmetric solver algorithm does not work with BFGS yet.\nThe full-Newton method will be used instead.");
//		}
	}
*/	
	// fix all mixture dofs that are not used
	// that is, that are not part of a biphasic, biphasic-solute, or triphasic or multiphasic element
	// this is done in three steps
	// step 1. mark all mixture nodes
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEDomain& dom = m_mesh.Domain(nd);
		FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(dom.GetMaterial());
		FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(dom.GetMaterial());
		FETriphasic*      btm = dynamic_cast<FETriphasic*     >(dom.GetMaterial());
		FEMultiphasic*    bmm = dynamic_cast<FEMultiphasic*   >(dom.GetMaterial());
		if (bm || bsm || btm || bmm)
		{
			for (i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (j=0; j<N; ++j) 
					if (m_mesh.Node(n[j]).m_ID[DOF_P] == 0) m_mesh.Node(n[j]).m_ID[DOF_P] = 1;
			}
		}
		if (bsm)
		{
			for (i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (j=0; j<N; ++j) {
					int dofc = DOF_C + bsm->m_pSolute->GetSoluteID();
					if (m_mesh.Node(n[j]).m_ID[dofc] == 0) m_mesh.Node(n[j]).m_ID[dofc] = 1;
				}
			}
		}
		else if (btm)
		{
			for (i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (j=0; j<N; ++j) {
					for (k=0; k<2; ++k) {
						int dofc = DOF_C + btm->m_pSolute[k]->GetSoluteID();
						if (m_mesh.Node(n[j]).m_ID[dofc] == 0) m_mesh.Node(n[j]).m_ID[dofc] = 1;
					}
				}
			}
		}
		else if (bmm)
		{
			for (i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				int nsol = bmm->m_pSolute.size();
				for (j=0; j<N; ++j) {
					for (k=0; k<nsol; ++k) {
						int dofc = DOF_C + btm->m_pSolute[k]->GetSoluteID();
						if (m_mesh.Node(n[j]).m_ID[dofc] == 0) m_mesh.Node(n[j]).m_ID[dofc] = 1;
					}
				}
			}
		}	
	}
	
	// step 2. fix mixture dofs of all unmarked nodes
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (j=0; j<N; ++j) {
					if (m_mesh.Node(n[j]).m_ID[DOF_P] != 1) m_mesh.Node(n[j]).m_ID[DOF_P] = -1;
					for (k=0; k<MAX_CDOFS; ++k) {
						int dofc = DOF_C + k;
						if (m_mesh.Node(n[j]).m_ID[dofc] != 1) m_mesh.Node(n[j]).m_ID[dofc] = -1;
					}
				}
			}
		}
		
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (j=0; j<N; ++j) {
					if (m_mesh.Node(n[j]).m_ID[DOF_P] != 1) m_mesh.Node(n[j]).m_ID[DOF_P] = -1;
					for (k=0; k<MAX_CDOFS; ++k) {
						int dofc = DOF_C + k;
						if (m_mesh.Node(n[j]).m_ID[dofc] != 1) m_mesh.Node(n[j]).m_ID[dofc] = -1;
					}
				}
			}
		}
	}
	
	// step 3. free all marked dofs
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (j=0; j<N; ++j) {
					if (m_mesh.Node(n[j]).m_ID[DOF_P] == 1) m_mesh.Node(n[j]).m_ID[DOF_P] = 0;
					for (k=0; k<MAX_CDOFS; ++k) {
						int dofc = DOF_C + k;
						if (m_mesh.Node(n[j]).m_ID[dofc] == 1) m_mesh.Node(n[j]).m_ID[dofc] = 0;
					}
				}
			}
		}
		
		
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (j=0; j<N; ++j) {
					if (m_mesh.Node(n[j]).m_ID[DOF_P] == 1) m_mesh.Node(n[j]).m_ID[DOF_P] = 0;
					for (k=0; k<MAX_CDOFS; ++k) {
						int dofc = DOF_C + k;
						if (m_mesh.Node(n[j]).m_ID[dofc] == 1) m_mesh.Node(n[j]).m_ID[dofc] = 0;
					}
				}
			}
		}
	}
	
	return true;
}

//=============================================================================
//                               S O L V E
//=============================================================================

//-----------------------------------------------------------------------------
//! This is the main solve method. This function loops over all analysis steps
//! and solves each one in turn. 
//! \sa FEAnalysisStep

bool FEBioModel::Solve()
{
	// echo fem data to the logfile
	// we do this here (and not e.g. directly after input)
	// since the data can be changed after input, which is the case,
	// for instance, in the parameter optimization module
	if (m_becho) echo_input(*this);

	// start the total time tracker
	m_TotalTime.start();

	// convergence flag
	bool bconv = true;

	// loop over all analysis steps
	// Note that we don't necessarily from step 0.
	// This is because the user could have restarted
	// the analysis. 
	for (size_t nstep=m_nStep; nstep < m_Step.size(); ++nstep)
	{
		// set the current analysis step
		m_nStep = nstep;
		m_pStep = m_Step[nstep];

		// intitialize step data
		if (m_pStep->Init() == false)
		{
			bconv = false;
			break;
		}

		// solve the analaysis step
		bconv = m_pStep->Solve();

		// break if the step has failed
		if (bconv == false) break;

		// wrap it up
		m_pStep->Finish();
	}

	// close the plot file
	if (m_plot) m_plot->Close();

	// stop total time tracker
	m_TotalTime.stop();

	// get and print elapsed time
	char sztime[64];
	m_TotalTime.time_str(sztime);
	clog.printf("\n Elapsed time : %s\n\n", sztime);

	if (bconv)
	{
		clog.printf("\n N O R M A L   T E R M I N A T I O N\n\n");
	}
	else
	{
		clog.printf("\n E R R O R   T E R M I N A T I O N\n\n");
	}

	// flush the log file
	clog.flush();

	// We're done !
	return bconv;
}
