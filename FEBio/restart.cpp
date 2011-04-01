// restart module
#include "stdafx.h"
#include "fem.h"
#include "FERigid.h"
#include "FERestartImport.h"
#include "FESlidingInterface.h"
#include "FETiedInterface.h"
#include "FERigidWallInterface.h"
#include "FEFacet2FacetSliding.h"
#include "FESlidingInterface2.h"
#include "FEPeriodicBoundary.h"
#include "FESurfaceConstraint.h"
#include "FETransverselyIsotropic.h"
#include "log.h"
#include "LSDYNAPlotFile.h"
#include "FEBioPlotFile.h"
#include "version.h"
#include "FEPressureLoad.h"
#include "FETractionLoad.h"
#include "FEFluidFlux.h"
#include "FEPoroTraction.h"
#include "FESoluteFlux.h"
#include "FEHeatFlux.h"

//-----------------------------------------------------------------------------
bool restart(FEM& fem, const char* szfile)
{
	// load restart data
	if (fem.Restart(szfile) == false) return false;

	// continue the analysis
	return fem.Solve();
}

//-----------------------------------------------------------------------------
//!  This routine reads a binary archive that stores a restart point and prepares
//!  the FEM data to be restarted from this point
//!	\param[in] szfile name of the file

bool FEM::Restart(const char* szfile)
{
	// check the extension of the file
	// if the extension is .dmp or not given it is assumed the file
	// is a bindary archive (dump file). Otherwise it is assumed the
	// file is a restart input file.
	const char* ch = strrchr(szfile, '.');
	if ((ch == 0) || (strcmp(ch, ".dmp") == 0) || (strcmp(ch, ".DMP") == 0))
	{
		// the file is binary so just read the dump file and return

		// open the archive
		DumpFile ar(this);
		if (ar.Open(szfile) == false) { fprintf(stderr, "FATAL ERROR: failed opening restart archive\n"); return false; }

		// read the archive
		try
		{
			if (Serialize(ar) == false) { fprintf(stderr, "FATAL ERROR: failed reading restart data from archive %s\n", szfile); return false; }
		}
		catch (...)
		{
			fprintf(stderr, "FATAL ERROR: failed reading restart data from archive %s\n", szfile); 
			return false;
		}
	}
	else
	{
		// the file is assumed to be a xml-text input file
		FERestartImport file;
		if (file.Load(*this, szfile) == false)
		{
			char szerr[256];
			file.GetErrorMessage(szerr);
			fprintf(stderr, "%s", szerr);
			return false;
		}
	}

	// Open the log file for appending
	if (clog.append(m_szlog) == false)
	{
		printf("WARNING: Could not reopen log file. A new log file is created\n");
		clog.open(m_szlog);
		return false;
	}

	// inform the user from where the problem is restarted
	clog.printbox(" - R E S T A R T -", "Restarting from time %lg.\n", m_ftime);

	return true;
}

//-----------------------------------------------------------------------------
//!  Reads or writes the current state to/from a binary file
//!  This is used to restart the solution from a saved position
//!  or to create a restart point.
//!  A version number is written to file to make sure the same
//!  format is used for reading and writing.
//! \param[in] ar the archive to which the data is serialized
//! \sa DumpFile

bool FEM::Serialize(DumpFile& ar)
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

void FEM::SerializeLoadData(DumpFile& ar)
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
void FEM::SerializeConstants(DumpFile& ar)
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

void FEM::SerializeAnalysisData(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		// analysis steps
		ar << (int) m_Step.size();
		for (int i=0; i<(int) m_Step.size(); ++i) m_Step[i]->Serialize(ar);
		ar << m_nStep;
		ar << m_ftime << m_ftime0;

		// direct solver data
		ar << m_nsolver;
		ar << m_neq;
		ar << m_npeq;
		ar << m_nceq;
		ar << m_bwopt;
		ar << m_bsymm;

		// body loads
		ar << (int) m_BF.size();
		for (int i=0; i<(int) m_BF.size(); ++i)
		{
			FEBodyForce* pbf = m_BF[i];
			int ntype = -1;
			if (dynamic_cast<FEConstBodyForce*>(pbf)) ntype = FE_CONST_BODY_FORCE;
			if (dynamic_cast<FENonConstBodyForce*>(pbf)) ntype = FE_NONCONST_BODY_FORCE;
			if (dynamic_cast<FECentrifugalBodyForce*>(pbf)) ntype = FE_CENTRIFUGAL_BODY_FORCE;
			assert(ntype);
			ar << ntype;
			pbf->Serialize(ar);
		}
	}
	else
	{
		m_Step.clear();

		// analysis steps
		int nsteps;
		ar >> nsteps;
		for (int i=0; i<nsteps; ++i)
		{
			FEAnalysis* pstep = new FEAnalysis(*this);
			pstep->Serialize(ar);
			m_Step.push_back(pstep);
		}
		ar >> m_nStep;
		ar >> m_ftime >> m_ftime0;

		// direct solver data
		ar >> m_nsolver;
		ar >> m_neq;
		ar >> m_npeq;
		ar >> m_nceq;
		ar >> m_bwopt;
		ar >> m_bsymm;

		// body loads
		int nbl;
		ar >> nbl;
		m_BF.clear();
		for (int i=0; i<nbl; ++i)
		{
			int ntype = -1;
			ar >> ntype;
			FEBodyForce* pbl = 0;
			switch (ntype)
			{
			case FE_CONST_BODY_FORCE: pbl = new FEConstBodyForce; break;
			case FE_NONCONST_BODY_FORCE: pbl = new FENonConstBodyForce; break;
			case FE_CENTRIFUGAL_BODY_FORCE: pbl = new FECentrifugalBodyForce; break;
			default:
				assert(false);
			}
			assert(pbl);
			pbl->Serialize(ar);
			m_BF.push_back(pbl);
		}

		// set the correct step
		m_pStep = m_Step[m_nStep];
	}
}

//-----------------------------------------------------------------------------

void FEM::SerializeMaterials(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		// store the nr of materials
		ar << Materials();

		// store the materials
		for (int i=0; i<Materials(); ++i)
		{
			FEMaterial* pmat = GetMaterial(i);

			// store the type string
			ar << pmat->GetTypeString();

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
			FEMaterial* pmat = FEMaterialFactory::CreateMaterial(szmat);
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
			AddParameterList(pmat->GetParameterList());
		}

		// we still need to reset the material pointers for the nested materials
		for (int i=0; i<nmat; ++i)
		{
			FENestedMaterial* pmat = dynamic_cast<FENestedMaterial*>(m_MAT[i]);
			if (pmat)
			{
				pmat->m_pBase = dynamic_cast<FESolidMaterial*>(m_MAT[pmat->m_nBaseMat-1]);
				assert(pmat->m_pBase);
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FEM::SerializeGeometry(DumpFile &ar)
{
	// serialize the mesh first 
	m_mesh.Serialize(ar);

	// serialize the other geometry data
	if (ar.IsSaving())
	{
		int i;

		// rigid bodies
		ar << m_nreq << m_nrm << m_nrb;
		for (i=0; i<m_nrb; ++i) m_RB[i].Serialize(ar);

		// rigid joints
		ar << m_nrj;
		for (i=0; i<m_nrj; ++i) m_RJ[i]->Serialize(ar);	
	}
	else
	{
		int i;

		// rigid bodies
		ar >> m_nreq >> m_nrm >> m_nrb;
		if (m_nrb) m_RB.resize(m_nrb); else m_RB.clear();
		for (i=0; i<m_nrb; ++i)
		{
			FERigidBody& rb = m_RB[i];
			rb.Serialize(ar);
			rb.AttachToFEM(this);
		}

		// rigid joints
		ar >> m_nrj;
		for (i=0; i<m_nrj; ++i)
		{
			FERigidJoint* prj = new FERigidJoint(this);
			prj->Serialize(ar);
			m_RJ.push_back(prj);
		}
	}
}

//-----------------------------------------------------------------------------

void FEM::SerializeContactData(DumpFile &ar)
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
			case FE_CONTACT_SLIDING    : ps = new FESlidingInterface(this); break;
			case FE_FACET2FACET_SLIDING: ps = new FEFacet2FacetSliding(this); break;
			case FE_CONTACT_TIED       : ps = new FETiedInterface(this); break;
			case FE_CONTACT_RIGIDWALL  : ps = new FERigidWallInterface(this); break;
			case FE_CONTACT_SLIDING2   : ps = new FESlidingInterface2(this); break;
			case FE_PERIODIC_BOUNDARY  : ps = new FEPeriodicBoundary(this); break;
			case FE_SURFACE_CONSTRAINT : ps = new FESurfaceConstraint(this); break;
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
void FEM::SerializeBoundaryData(DumpFile& ar)
{
	int i, n;

	if (ar.IsSaving())
	{
		// displacements
		ar << (int) m_DC.size();
		for (i=0; i<(int) m_DC.size(); ++i) 
		{
			FENodalDisplacement& dc = *m_DC[i];
			ar << dc.bc << dc.lc << dc.node << dc.s;
		}

		// nodal loads
		ar << (int) m_FC.size();
		for (i=0; i<(int) m_FC.size(); ++i)
		{
			FENodalForce& fc = *m_FC[i];
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
			assert(ntype != -1);
			ar << ntype;
			psl->Serialize(ar);
		}

		// rigid body displacements
		ar << m_RDC.size();
		for (i=0; i<(int) m_RDC.size(); ++i)
		{
			FERigidBodyDisplacement& dc = *m_RDC[i];
			ar << dc.bc << dc.id << dc.lc << dc.sf;
		}

		// rigid body forces
		ar << m_RFC.size();
		for (i=0; i<(int) m_RFC.size(); ++i)
		{
			FERigidBodyForce& fc = *m_RFC[i];
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
		n = (int) m_LCSet.size();
		ar << n;
		if (m_LCSet.empty() == false)
		{
			list<FELinearConstraintSet*>::iterator ic = m_LCSet.begin();
			for (i=0; i<n; ++i, ++ic) (*ic)->Serialize(ar);
		}
	}
	else
	{
		int n;
		// displacements
		ar >> n;
		m_DC.clear();
		for (i=0; i<n; ++i) 
		{
			FENodalDisplacement* pdc = new FENodalDisplacement;
			ar >> pdc->bc >> pdc->lc >> pdc->node >> pdc->s;
			m_DC.push_back(pdc);
		}
		
		// nodal loads
		ar >> n;
		m_FC.clear();
		for (i=0; i<n; ++i)
		{
			FENodalForce* pfc = new FENodalForce;
			ar >> pfc->bc >> pfc->lc >> pfc->node >> pfc->s;
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
			case FE_PRESSURE_LOAD: ps = new FEPressureLoad      (psurf); break;
			case FE_TRACTION_LOAD: ps = new FETractionLoad      (psurf); break;
			case FE_FLUID_FLUX   : ps = new FEFluidFlux         (psurf); break;
			case FE_PORO_TRACTION: ps = new FEPoroNormalTraction(psurf); break;
			case FE_SOLUTE_FLUX  : ps = new FESoluteFlux        (psurf); break;
			case FE_HEAT_FLUX    : ps = new FEHeatFlux          (psurf); break;
			default:
				assert(false);
			}
			assert(ps);
			ps->Serialize(ar);
			m_SL.push_back(ps);
		}

		// rigid body displacements
		ar >> n;
		m_RDC.clear();
		for (i=0; i<n; ++i)
		{
			FERigidBodyDisplacement* pdc = new FERigidBodyDisplacement;
			ar >> pdc->bc >> pdc->id >> pdc->lc >> pdc->sf;
			m_RDC.push_back(pdc);
		}

		// rigid body forces
		ar >> n;
		m_RFC.clear();
		for (i=0; i<n; ++i)
		{
			FERigidBodyForce* pfc = new FERigidBodyForce;
			ar >> pfc->bc >> pfc->id >> pfc->lc >> pfc->sf;
			m_RFC.push_back(pfc);
		}

		// rigid nodes
		ar >> n;
		m_RN.clear();
		for (i=0; i<n; ++i)
		{
			int nid; bool bactive;
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
		for (i=0; i<n; ++i)
		{
			FELinearConstraintSet* plc = new FELinearConstraintSet(this);
			plc->Serialize(ar);
			m_LCSet.push_back(plc);
		}
	}
}

//-----------------------------------------------------------------------------
//! Serialization of FEM data
void FEM::SerializeIOData(DumpFile &ar)
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
		m_Data.Serialize(ar);
	}
	else
	{
		// file names
		ar >> m_szfile >> m_szplot >> m_szlog >> m_szdump;
		ar >> m_sztitle;

		// don't forget to call SetInputFilename so
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
				m_plot = new FEBioPlotFile;
				if (m_plot->Append(*this, m_szplot) == false)
				{
					printf("FATAL ERROR: Failed reopening plot database %s\n", m_szplot);
					throw "FATAL ERROR";
				}
			}
			break;
		};

		// data records
		m_Data.Serialize(ar);
	}
}
