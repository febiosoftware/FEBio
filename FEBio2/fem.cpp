#include "stdafx.h"
#include "fem.h"
#include <string.h>
#include "FECore/DumpFile.h"
#include "FECore/XMLReader.h"
#include "FESlidingInterface.h"
#include "FETiedInterface.h"
#include "FERigidWallInterface.h"
#include "FEFacet2FacetSliding.h"
#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FEPeriodicBoundary.h"
#include "FESurfaceConstraint.h"
#include "log.h"
#include "plugin.h"
#include "LSDYNAPlotFile.h"
#include "FEBiphasic.h"

// --- Global Constants Data ---
// m_Const needs a definition, since static
map<string, double> FEM::m_Const;

//-----------------------------------------------------------------------------
//! Constructor of the FEM class
//! Initializes default variables
//!
FEM::FEM()
{
	// --- Analysis Data ---
	m_pStep = 0;
	m_nStep = -1;
	m_bsym_poro = true;			// use symmetric poro implementation
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

	// --- Geometry Data ---
	m_nreq = 0;
	m_nrb = 0;
	m_nrm = 0;

	m_bsymm = true;	// assume symmetric stiffness matrix
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

	m_neq = 0;
	m_npeq = 0;
	m_nceq = 0;
	m_bwopt = 0;

	// --- I/O-Data ---
	strcpy(m_szplot, "n3plot");
	strcpy(m_szlog , "n3log" );
	strcpy(m_szdump, "n3dump");

	m_sztitle[0] = 0;
	m_debug = false;

	m_plot = 0;
}

//-----------------------------------------------------------------------------
//! copy constructor for FEM class.
//! The copy constructor and assignment operator are only used for push/pop'ing
//! to a stack. This is used by the running restart feature. Note that not all 
//! data is copied. We only copy the data that changes during iterations.
//! for push/pop'ing

FEM::FEM(const FEM& fem)
{
	ShallowCopy(const_cast<FEM&>(fem));
}

//-----------------------------------------------------------------------------
//! assignment operator for FEM class.
//! The copy constructor and assignment operator are only used for push/pop'ing
//! to a stack. This is used by the running restart feature. Note that not all 
//! data is copied. We only copy the data that changes during iterations.
//! for push/pop'ing

void FEM::operator =(const FEM& fem)
{
	ShallowCopy(const_cast<FEM&>(fem));
}

//-----------------------------------------------------------------------------
//! destructor of FEM class.
//! Delete all dynamically allocated data

FEM::~FEM()
{
	size_t i;
	for (i=0; i<m_Step.size(); ++i) delete m_Step[i]; m_Step.clear();
	for (i=0; i<m_RJ.size  (); ++i) delete m_RJ [i] ; m_RJ.clear  ();
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

	if (!m_LCSet.empty())
	{
		list<FELinearConstraintSet*>::iterator pi;
		for (pi = m_LCSet.begin(); pi != m_LCSet.end(); ++pi) delete (*pi); 
		m_LCSet.clear();
	}
}

//-----------------------------------------------------------------------------
// This function is used when pushing the FEM to a stack. Since we don't need
// to copy all the data, this function only copies the data that needs to be 
// restored for a running restart.
//
void FEM::ShallowCopy(FEM& fem)
{
	int i;

	// copy time data
	m_ftime = fem.m_ftime;

	// copy the mesh
	m_mesh = fem.m_mesh;

	// copy rigid body data
	m_RB = fem.m_RB;

	// copy rigid joint data
	if (m_RJ.empty())
	{
		for (i=0; i<(int) fem.m_RJ.size(); ++i) m_RJ.push_back(new FERigidJoint(this));
	}
	assert(m_RJ.size() == fem.m_RJ.size());
	for (i=0; i<(int) m_RJ.size(); ++i) m_RJ[i]->ShallowCopy(*fem.m_RJ[i]);

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
// This function adds a callback routine
//
void FEM::AddCallback(FEBIO_CB_FNC pcb, void *pd)
{
	FEBIO_CALLBACK cb;
	cb.m_pcb = pcb;
	cb.m_pd = pd;

	m_pcb.push_back(cb);
}

//-----------------------------------------------------------------------------
// FEM::DoCallback
// Call the callback function if there is one defined
//

void FEM::DoCallback()
{
	list<FEBIO_CALLBACK>::iterator it = m_pcb.begin();
	for (int i=0; i<(int) m_pcb.size(); ++i, ++it)
	{
		// call the callback function
		(it->m_pcb)(this, it->m_pd);
	}
}

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable

//! This function returns a pointer to a named variable. Currently, we only
//! support names of the form material_name.parameter name. The material
//! name is a user defined name for a material and the parameter name is the
//! predefined name of the variable
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
		index = atoi(szarg);
	}

	// get the material's parameter list
	FEParameterList& pl = pmat->GetParameterList();

	// find the parameter
	FEParam* pp = pl.Find(szvar);

	if (pp)
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


	// the rigid bodies are dealt with differently
	for (i=0; i<m_nrb; ++i)
	{
		FERigidBody& rb = m_RB[i];

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
void FEM::EvalParameterList(FEParameterList &pl)
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
void FEM::SetGlobalConstant(const string& s, double v)
{
	m_Const[s] = v;
	return;
}

//-----------------------------------------------------------------------------
double FEM::GetGlobalConstant(const string& s)
{
	return m_Const.find(s)->second;
}

//-----------------------------------------------------------------------------
//! Returns the elastic component of the material, identified by id
FEElasticMaterial* FEM::GetElasticMaterial(int id)
{
	FEMaterial* pm = m_MAT[id];
	while (dynamic_cast<FENestedMaterial*>(pm)) pm = (dynamic_cast<FENestedMaterial*>(pm))->m_pBase;
	while (dynamic_cast<FEBiphasic*>(pm)) pm = (dynamic_cast<FEBiphasic*>(pm))->m_pSolid;
	while (dynamic_cast<FEBiphasicSolute*>(pm)) pm = (dynamic_cast<FEBiphasicSolute*>(pm))->m_pSolid;
	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
	assert(pme);
	return pme;
}

//-----------------------------------------------------------------------------
//! Returns the elastic component of the material pm
FEElasticMaterial* FEM::GetElasticMaterial(FEMaterial* pm)
{
	while (dynamic_cast<FENestedMaterial*>(pm)) pm = (dynamic_cast<FENestedMaterial*>(pm))->m_pBase;
	while (dynamic_cast<FEBiphasic*>(pm)) pm = (dynamic_cast<FEBiphasic*>(pm))->m_pSolid;
	while (dynamic_cast<FEBiphasicSolute*>(pm)) pm = (dynamic_cast<FEBiphasicSolute*>(pm))->m_pSolid;
	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
	assert(pme);
	return pme;
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
