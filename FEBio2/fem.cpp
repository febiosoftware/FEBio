#include "stdafx.h"
#include "fem.h"
#include <string.h>
#include "FECore/DumpFile.h"
#include "FECore/XMLReader.h"
#include "FEBioLib/FESlidingInterface.h"
#include "FEBioLib/FETiedInterface.h"
#include "FEBioLib/FERigidWallInterface.h"
#include "FEBioLib/FEFacet2FacetSliding.h"
#include "FEBioLib/FESlidingInterface2.h"
#include "FEBioLib/FESlidingInterface3.h"
#include "FEBioLib/FEPeriodicBoundary.h"
#include "FEBioLib/FESurfaceConstraint.h"
#include "FEBioLib/log.h"
#include "plugin.h"
#include "LSDYNAPlotFile.h"
#include "FECore/MathParser.h"
#include "FEBioLib/FEBiphasic.h"
#include "FEBioLib/FEBiphasicSolute.h"
#include "FEBioLib/FETriphasic.h"
#include "FEBioLib/FETransverselyIsotropic.h"
#include "Interrupt.h"
#include "console.h"

// --- Global Constants Data ---
// m_Const and m_SD need a definitions, since static
map<std::string, double> FEModel::m_Const;
vector<FESoluteData*> FEM::m_SD;

//-----------------------------------------------------------------------------
void FEBioProgress::SetProgress(double f)
{
	// get the number of steps
	int nsteps = m_fem.Steps();

	// check debug flag
	bool bdebug = m_fem.GetDebugFlag();

	// obtain a pointer to the console object. We'll use this to
	// set the title of the console window.
	Console* pShell = Console::GetHandle();

	// print progress in title bar
	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s", m_fem.m_nStep+1, nsteps, f, m_fem.GetFileTitle(), (bdebug?"FEBio (debug mode)": "FEBio"));
	else
		pShell->SetTitle("(%.f%%) %s - %s", f, m_fem.GetFileTitle(), (bdebug?"FEBio (debug mode)": "FEBio"));
}

//-----------------------------------------------------------------------------
//! Constructor of the FEM class
//! Initializes default variables
//!
FEM::FEM()
{
	// --- Analysis Data ---
	m_nStep = -1;
	m_nplane_strain = -1;	// don't use plain strain mode

	m_ftime = 0;
	m_ftime0 = 0;

	// add the "zero" loadcurve
	// this is the loadcurve that will be used if a loadcurve is not
	// specified for something that depends on time
	// TODO: I want to get rid of this 
	FELoadCurve* plc = new FELoadCurve();
	plc->Create(2);
	plc->LoadPoint(0).time = 0;
	plc->LoadPoint(0).value = 0;
	plc->LoadPoint(1).time = 1;
	plc->LoadPoint(1).value = 1;
	plc->SetExtendMode(FELoadCurve::EXTRAPOLATE);
	AddLoadCurve(plc);

	m_bInterruptable = true;

	// --- Material Data ---
	// (nothing to initialize yet)

	// --- Load Curve Data ---
	// (nothing to initialize yet)

	// --- Direct Solver Data ---
	// set the default linear solver
#ifdef PARDISO
	m_nsolver = PARDISO_SOLVER;
#elif PSLDLT
	m_nsolver = PSLDLT_SOLVER;
#else
	m_nsolver = SKYLINE_SOLVER;
#endif

	m_bwopt = 0;

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
//! destructor of FEM class.
//! Delete all dynamically allocated data

FEM::~FEM()
{
	size_t i;
	for (i=0; i<m_Step.size(); ++i) delete m_Step[i]; m_Step.clear();
	for (i=0; i<m_CI.size  (); ++i) delete m_CI [i] ; m_CI.clear  ();
	for (i=0; i<m_MAT.size (); ++i) delete m_MAT[i] ; m_MAT.clear ();
	for (i=0; i<m_LC.size  (); ++i) delete m_LC [i] ; m_LC.clear  ();
	for (i=0; i<m_BF.size  (); ++i) delete m_BF [i] ; m_BF.clear  ();
	for (i=0; i<m_DC.size  (); ++i) delete m_DC [i] ; m_DC.clear  ();
	for (i=0; i<m_FC.size  (); ++i) delete m_FC [i] ; m_FC.clear  ();
	for (i=0; i<m_SL.size  (); ++i) delete m_SL [i] ; m_SL.clear  ();
	for (i=0; i<m_RDC.size (); ++i) delete m_RDC[i] ; m_RDC.clear ();
	for (i=0; i<m_RFC.size (); ++i) delete m_RFC[i] ; m_RFC.clear ();
	for (i=0; i<m_RN.size  (); ++i) delete m_RN [i] ; m_RN.clear  ();
	for (i=0; i<m_NLC.size (); ++i) delete m_NLC[i] ; m_NLC.clear ();
	for (i=0; i<m_Obj.size (); ++i) delete m_Obj[i] ; m_Obj.clear ();
//	for (i=0; i<m_PC.size  (); ++i) delete m_PC [i] ; m_PC.clear  ();
//	for (i=0; i<m_LCSet.size(); ++i) delete m_LCSet[i]; m_LCSet.clear();
}

//-----------------------------------------------------------------------------
static FEM* pfem_copy = 0;

//-----------------------------------------------------------------------------
void FEM::PushState()
{
	if (pfem_copy == 0) pfem_copy = new FEM;
	pfem_copy->ShallowCopy(*this);
}

//-----------------------------------------------------------------------------
void FEM::PopState()
{
	assert(pfem_copy);
	ShallowCopy(*pfem_copy);
	delete pfem_copy;
	pfem_copy = 0;
}

//-----------------------------------------------------------------------------
void FEM::Write()
{
	m_plot->Write(*this);
}

//-----------------------------------------------------------------------------
void FEM::WriteData()
{
	m_Data.Write();
}

//-----------------------------------------------------------------------------
// This function is used when pushing the FEM state data. Since we don't need
// to copy all the data, this function only copies the data that needs to be 
// restored for a running restart.

// TODO: Shallow copy nonlinear constraints
void FEM::ShallowCopy(FEM& fem)
{
	int i;

	// copy time data
	m_ftime = fem.m_ftime;

	// copy the mesh
	m_mesh = fem.m_mesh;

	// copy rigid body data
	if (m_Obj.empty())
	{
		for (int i=0; i<(int) fem.m_Obj.size();	++i)
		{
			FERigidBody* prb = new FERigidBody(this);
			m_Obj.push_back(prb);
		}
	}
	assert(m_Obj.size() == fem.m_Obj.size());
	for (i=0; i<(int) m_Obj.size(); ++i) m_Obj[i]->ShallowCopy(fem.m_Obj[i]);

	// copy contact data
	if (m_CI.empty())
	{
		FEContactInterface* pci;
		for (int i=0; i<fem.ContactInterfaces(); ++i)
		{
			switch (fem.m_CI[i]->Type())
			{
			case FE_CONTACT_SLIDING    : pci = new FESlidingInterface  (this); break;
			case FE_FACET2FACET_SLIDING: pci = new FEFacet2FacetSliding(this); break;
			case FE_CONTACT_TIED       : pci = new FETiedInterface     (this); break;
			case FE_CONTACT_RIGIDWALL  : pci = new FERigidWallInterface(this); break;
			case FE_CONTACT_SLIDING2   : pci = new FESlidingInterface2 (this); break;
			case FE_PERIODIC_BOUNDARY  : pci = new FEPeriodicBoundary  (this); break;
			case FE_SURFACE_CONSTRAINT : pci = new FESurfaceConstraint (this); break;
			case FE_CONTACT_SLIDING3   : pci = new FESlidingInterface3 (this); break;
			default:
				assert(false);
			}
			m_CI.push_back(pci);
		}
	}
	assert(ContactInterfaces() == fem.ContactInterfaces());
	for (i=0; i<ContactInterfaces(); ++i) m_CI[i]->ShallowCopy(*fem.m_CI[i]);
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

double* FEM::FindParameter(const char* szparam)
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
	int index = -1;
	if (szarg)
	{
		*szarg = 0; szarg++;
		const char* ch = strchr(szarg, ']');
		assert(ch);
		index = atoi(szarg) - 1;	// index is one-based for user
	}

	// get the material's parameter list
	FEParameterList& pl = pmat->GetParameterList();

	// find the parameter
	FEParam* pp = pl.Find(szvar);

	if (pp) return ReturnParameter(pp, index);

	// if this material is a nested material, we'll need to check the base material
	FENestedMaterial* pmn = dynamic_cast<FENestedMaterial*>(pmat);
	if (pmn)
	{
		char* ch = strchr((char*)szvar, '.');
		if (ch == 0) return 0;
		*ch = 0;
		const char* szvar2 = ch+1;

		if (strcmp(szvar, "elastic") == 0)
		{
			// search the nested material parameter list
			pmat = pmn->m_pBase;
			FEParameterList& pl = pmat->GetParameterList();
			FEParam* pp = pl.Find(szvar2);
			if (pp) return ReturnParameter(pp, index);
			else return 0;
		}
	}

	// if material is solid mixture, check individual solid materials
	FEElasticMixture* pme = dynamic_cast<FEElasticMixture*>(pmat);
	if (pme) return FindSolidMixtureParameter(szvar, index, pme);
	FEUncoupledElasticMixture* pmu = dynamic_cast<FEUncoupledElasticMixture*>(pmat);
	if (pmu) return FindUncoupledSolidMixtureParameter(szvar, index, pmu);
	
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
			FEParameterList& pl = pme->GetParameterList();
			FEParam* pp = pl.Find(szvar2);
			if (pp) return ReturnParameter(pp, index);
			// if material is solid mixture, check individual solid materials
			FEElasticMixture* pmm = dynamic_cast<FEElasticMixture*>(pme);
			if (pmm) return FindSolidMixtureParameter(szvar2, index, pmm);
			else return 0;
		}
		else if (strcmp(szvar, "permeability") == 0)
		{
			// search the nested material parameter list
			pmat = pmb->m_pPerm;
			FEParameterList& pl = pmat->GetParameterList();
			FEParam* pp = pl.Find(szvar2);
			if (pp) return ReturnParameter(pp, index);
			else return 0;
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
//! Return a pointer to the parameter variable

double* FEM::ReturnParameter(FEParam* pp, const int index)
{
	switch (pp->m_itype)
	{
		case FE_PARAM_DOUBLE:
		{
			assert(index<0);
			return &pp->value<double>();
		}
			break;
		case FE_PARAM_DOUBLEV:
		{
			assert((index >= 0) && (index < pp->m_ndim));
			return &(pp->pvalue<double>()[index]);
		}
			break;
		default:
			// other types are not supported yet
			assert(false);
			return 0;
	}
}

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable in a solid mixture

double* FEM::FindSolidMixtureParameter(const char* szvar, const int index, FEElasticMixture* pme)
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
			if (pp) return ReturnParameter(pp, index);
			else return 0;
		}
	}
	// no match found
	return 0;
}

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable in a uncoupled solid mixture

double* FEM::FindUncoupledSolidMixtureParameter(const char* szvar, const int index, FEUncoupledElasticMixture* pme)
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
			if (pp) return ReturnParameter(pp, index);
			else return 0;
		}
	}
	// no match found
	return 0;
}

//-----------------------------------------------------------------------------
//! Evaluate a parameter list
void FEM::EvaluateParameterList(FEParameterList &pl)
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
void FEM::EvaluateMaterialParameters(FEMaterial* pm)
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
		for (int i=0; i < (int) pem->m_pMat.size(); ++i)
			EvaluateMaterialParameters(pem->m_pMat[i]);
	}
	
	FEUncoupledElasticMixture* pum = dynamic_cast<FEUncoupledElasticMixture*>(pm);
	if (pum)
	{
		for (int i=0; i < (int) pum->m_pMat.size(); ++i)
			EvaluateMaterialParameters(pum->m_pMat[i]);
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
		EvaluateMaterialParameters(pbs->m_pSolid);
		EvaluateMaterialParameters(pbs->m_pPerm );
		EvaluateMaterialParameters(pbs->m_pDiff );
		EvaluateMaterialParameters(pbs->m_pSolub);
		EvaluateMaterialParameters(pbs->m_pOsmC );
	}
}

//-----------------------------------------------------------------------------
//! Reads the FEBio configuration file. This file contains some default settings.

bool FEM::Configure(const char *szfile)
{
	// open the configuration file
	XMLReader xml;
	if (xml.Open(szfile) == false)
	{
		fprintf(stderr, "FATAL ERROR: Failed reading FEBio configuration file %s.\n", szfile);
		return false;
	}

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_config", tag) == false) return false;

		if (strcmp(tag.m_szatv[0], "1.0") == 0)
		{
			if (!tag.isleaf())
			{
				// Read version 1.0
				++tag;
				do
				{
					if (tag == "linear_solver")
					{
						const char* szt = tag.AttributeValue("type");
						if      (strcmp(szt, "skyline"           ) == 0) m_nsolver = SKYLINE_SOLVER;
						else if (strcmp(szt, "psldlt"            ) == 0) m_nsolver = PSLDLT_SOLVER;
						else if (strcmp(szt, "superlu"           ) == 0) m_nsolver = SUPERLU_SOLVER;
						else if (strcmp(szt, "superlu_mt"        ) == 0) m_nsolver = SUPERLU_MT_SOLVER;
						else if (strcmp(szt, "pardiso"           ) == 0) m_nsolver = PARDISO_SOLVER;
						else if (strcmp(szt, "rcicg"             ) == 0) m_nsolver = RCICG_SOLVER;
						else if (strcmp(szt, "wsmp"              ) == 0) m_nsolver = WSMP_SOLVER;
					}
					else if (tag == "import")
					{
						const char* szfile = tag.szvalue();
						if (LoadPlugin(szfile) == false) throw XMLReader::InvalidValue(tag);
						printf("Plugin \"%s\" loaded successfully\n", szfile);
					}
					else throw XMLReader::InvalidTag(tag);

					// go to the next tag
					++tag;
				}
				while (!tag.isend());
			}
		}
		else
		{
			clog.printbox("FATAL ERROR", "Invalid version for FEBio configuration file.");
			return false;
		}
	}
	catch (XMLReader::XMLSyntaxError)
	{
		clog.printf("FATAL ERROR: Syntax error (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		clog.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		clog.printf("FATAL ERROR: unrecognized value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		clog.printf("FATAL ERROR: the value for tag \"%s\" is invalid (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::MissingAttribute e)
	{
		clog.printf("FATAL ERROR: Missing attribute \"%s\" of tag \"%s\" (line %d)\n", e.szatt, e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::UnmatchedEndTag e)
	{
		const char* sz = e.tag.m_szroot[e.tag.m_nlevel];
		clog.printf("FATAL ERROR: Unmatched end tag for \"%s\" (line %d)\n", sz, e.tag.m_nstart_line);
		return false;
	}
	catch (...)
	{
		clog.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	xml.Close();

	return true;
}

//-----------------------------------------------------------------------------

FEBoundaryCondition* FEM::FindBC(int nid)
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

//-----------------------------------------------------------------------------
void FEM::SetPlotFileNameExtension(const char *szext)
{
	char* ch = strrchr(m_szplot, '.');
	if (ch) *ch = 0;
	strcat(m_szplot, szext);
}

//-----------------------------------------------------------------------------
void FEM::SetSD(FESoluteData* psd)
{
	m_SD.push_back(psd);
}

//-----------------------------------------------------------------------------
FESoluteData* FEM::FindSD(int nid)
{
	int i;
	for (i=0; i<(int) m_SD.size(); ++i) if (m_SD[i]->m_nID == nid) return m_SD[i];
	
	return 0;
}

//-----------------------------------------------------------------------------
//! Sets the name of the FEBio input file
void FEM::SetInputFilename(const char* szfile)
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
void FEM::CheckInterruption()
{
	if (m_bInterruptable)
	{
		Interruption itr;
		if (itr.m_bsig)
		{
			itr.m_bsig = false;
			itr.interrupt();
		}
	}
}
