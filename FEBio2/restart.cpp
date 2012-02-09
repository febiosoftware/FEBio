// restart module
#include "stdafx.h"
#include "fem.h"
#include "FEBioLib/FERigid.h"
#include "FERestartImport.h"
#include "FEBioLib/FESlidingInterface.h"
#include "FEBioLib/FETiedInterface.h"
#include "FEBioLib/FERigidWallInterface.h"
#include "FEBioLib/FEFacet2FacetSliding.h"
#include "FEBioLib/FESlidingInterface2.h"
#include "FEBioLib/FESlidingInterface3.h"
#include "FEBioLib/FEPeriodicBoundary.h"
#include "FEBioLib/FESurfaceConstraint.h"
#include "FEBioLib/FETransverselyIsotropic.h"
#include "FEBioLib/FEPressureLoad.h"
#include "FEBioLib/FETractionLoad.h"
#include "FEBioLib/FEFluidFlux.h"
#include "FEBioLib/FEPoroTraction.h"
#include "FEBioLib/FESoluteFlux.h"
#include "FEBioLib/FEHeatFlux.h"
#include "FEBioLib/FEAnalysisStep.h"
#include "FEBioLib/FEElasticSolidDomain.h"
#include "FEBioLib/FEElasticShellDomain.h"
#include "FEBioLib/FEElasticTrussDomain.h"
#include "FEBioLib/FEHeatSolidDomain.h"
#include "FEBioLib/FEDiscreteSpringDomain.h"
#include "FEBioLib/FEBiphasicSolidDomain.h"
#include "FEBioLib/FEBiphasicSoluteDomain.h"
#include "FEBioLib/FEUDGHexDomain.h"
#include "FEBioLib/FERigidSolidDomain.h"
#include "FEBioLib/FERigidShellDomain.h"
#include "FEBioLib/FE3FieldElasticSolidDomain.h"
#include "FEBioLib/FEUT4Domain.h"
#include "FEBioLib/FEConstBodyForce.h"
#include "FEBioLib/FEPointConstraint.h"
#include "FEBioLib/FEAugLagLinearConstraint.h"
#include "FEBioLib/log.h"
#include "LSDYNAPlotFile.h"
#include "FEBioPlotFile.h"
#include "version.h"

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
		ar << (int) m_BF.size();
		for (int i=0; i<(int) m_BF.size(); ++i)
		{
			FEBodyForce* pbf = m_BF[i];
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
		FEM* pfem = dynamic_cast<FEM*>(ar.GetFEModel());

		// analysis steps
		int nsteps, ntype;
		ar >> nsteps;
		for (int i=0; i<nsteps; ++i)
		{
			ar >> ntype;
			FEAnalysisStep* pstep = 0;
			switch (ntype)
			{
			case FE_SOLID       : pstep = new FESolidAnalysis         (*this); break;
			case FE_BIPHASIC    : pstep = new FEBiphasicAnalysis      (*this); break;
			case FE_HEAT        : pstep = new FEHeatTransferAnalysis  (*this); break;
			case FE_POROSOLUTE  : pstep = new FEBiphasicSoluteAnalysis(*this); break;
			case FE_LINEAR_SOLID: pstep = new FELinearSolidAnalysis   (*this); break;
			case FE_HEAT_SOLID  : pstep = new FEThermoElasticAnalysis (*this); break;
			case FE_TRIPHASIC   : pstep = new FETriphasicAnalysis     (*this); break;
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
		m_BF.clear();
		for (int i=0; i<nbl; ++i)
		{
			int ntype = -1;
			ar >> ntype;
			FEBodyForce* pbl = 0;
			switch (ntype)
			{
			case FE_CONST_BODY_FORCE: pbl = new FEConstBodyForce(pfem); break;
			case FE_NONCONST_BODY_FORCE: pbl = new FENonConstBodyForce(pfem); break;
			case FE_CENTRIFUGAL_BODY_FORCE: pbl = new FECentrifugalBodyForce(pfem); break;
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
// TODO: serialize nonlinear constraints
void FEM::SerializeGeometry(DumpFile &ar)
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

void FEM::SerializeMesh(DumpFile& ar)
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
			case FE_CONTACT_SLIDING3   : ps = new FESlidingInterface3(this); break;
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
				ar << m_NLC[i]->Type();
				m_NLC[i]->Serialize(ar);
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
		int ntype;
		m_NLC.clear();
		for (i=0; i<n; ++i)
		{
			ar >> ntype;
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
