#include "stdafx.h"
#include "FEBioModel.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioXML/FERestartImport.h"
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"
#include "FECore/ObjectDataRecord.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FECoreKernel.h"
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
	m_debug = false;
	m_becho = true;
	m_plot = 0;
}

//-----------------------------------------------------------------------------
FEBioModel::~FEBioModel()
{
	// close the plot file
	if (m_plot) { delete m_plot; m_plot = 0; }
}

//-----------------------------------------------------------------------------
Timer& FEBioModel::GetTotalTimer()
{
	return m_TotalTime;
}

//=============================================================================
//
//		FEBioModel: I-O Functions
//
//=============================================================================

//-----------------------------------------------------------------------------
//! Return the data store
DataStore& FEBioModel::GetDataStore()
{
	return m_Data;
}

//-----------------------------------------------------------------------------
//! Add a data record to the data store
void FEBioModel::AddDataRecord(DataRecord* pd)
{
	m_Data.AddRecord(pd); 
}

//-----------------------------------------------------------------------------
//! Get the plot file
PlotFile* FEBioModel::GetPlotFile()
{
	return m_plot;
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
//! Set the name of the log file
void FEBioModel::SetLogFilename(const char* szfile) 
{ 
	strcpy(m_szlog , szfile); 
}

//-----------------------------------------------------------------------------
//! Set the name of the plot file
void FEBioModel::SetPlotFilename(const char* szfile) 
{ 
	strcpy(m_szplot, szfile); 
}

//-----------------------------------------------------------------------------
//! Set the name of the restart archive (i.e. the dump file)
void FEBioModel::SetDumpFilename (const char* szfile) 
{ 
	strcpy(m_szdump, szfile); 
}

//-----------------------------------------------------------------------------
//! Return the name of the input file
const char* FEBioModel::GetInputFileName()
{ 
	return m_szfile; 
}

//-----------------------------------------------------------------------------
//! Return the name of the log file
const char* FEBioModel::GetLogfileName()
{ 
	return m_szlog;  
}

//-----------------------------------------------------------------------------
//! Return the name of the plot file
const char* FEBioModel::GetPlotFileName()
{ 
	return m_szplot; 
}

//-----------------------------------------------------------------------------
//! get the file title (i.e. name of input file without the path)
//! \todo Do I actually need to store this?
const char* FEBioModel::GetFileTitle()
{ 
	return m_szfile_title; 
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

	// set the plot file
	FEBioPlotFile* pplt = new FEBioPlotFile(*this);
	m_plot = pplt;

	// define the plot file variables
	FEMesh& mesh = GetMesh();
	int NP = (int) fim.m_plot.size();
	for (int i=0; i<NP; ++i)
	{
		FEFEBioImport::FEPlotVariable& var = fim.m_plot[i];

		vector<int> item = var.m_item;
		if (item.empty() == false)
		{
			// TODO: currently, this is only supported for domain variables, where
			//       the list is a list of materials
			vector<int> lmat = var.m_item;

			// convert the material list to a domain list
			mesh.DomainListFromMaterial(lmat, item);
		}

		// add the plot output variable
		if (pplt->AddVariable(var.m_szvar, item) == false) 
		{
			felog.printf("FATAL ERROR: Output variable \"%s\" is not defined\n", var.m_szvar);
			return false;
		}
	}

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
	if (m_plot) m_plot->Write(*this);
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
		felog.printf("WARNING: Failed creating restart point.\n");
	}
	else 
	{
		Serialize(ar);
		felog.printf("\nRestart point created. Archive name is %s\n", m_szdump);
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

	// --- Global Data ---
	SerializeGlobals(ar);

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
//! Serialize global data
void FEBioModel::SerializeGlobals(DumpFile& ar)
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
		int nGD = GlobalDataItems();
		ar << nGD;
		for (int i=0; i<nGD; i++) 
		{
			FEGlobalData* pgd = GetGlobalData(i);
			ar << pgd->GetTypeStr();
			pgd->Serialize(ar);
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
		int nGD;
		ar >> nGD;
		if (nGD) 
		{
			char sztype[256];
			for (int i=0; i<nGD; ++i)
			{
				ar >> sztype;
				FEGlobalData* pgd = fecore_new<FEGlobalData>(FEGLOBALDATA_ID, sztype, this);
				pgd->Serialize(ar);
				FEModel::AddGlobalData(pgd);
			}
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
			ar << m_Step[i]->GetTypeStr();
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
			FEBodyLoad* pbl = m_BL[i];
			ar << pbl->GetTypeStr();
			pbl->Serialize(ar);
		}
	}
	else
	{
		m_Step.clear();
		FEModel* pfem = ar.GetFEModel();

		char sztype[256] = {0};

		// analysis steps
		int nsteps;
		ar >> nsteps;
		for (int i=0; i<nsteps; ++i)
		{
			ar >> sztype;
			FEAnalysis* pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, sztype, this); assert(pstep);
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
		char szbl[256] = {0};
		for (int i=0; i<nbl; ++i)
		{
			ar >> szbl;
			FEBodyLoad* pbl = fecore_new<FEBodyLoad>(FEBODYLOAD_ID, szbl, this);
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
	FECoreKernel& febio = FECoreKernel::GetInstance();

	if (ar.IsSaving())
	{
		// store the nr of materials
		ar << Materials();

		// store the materials
		for (int i=0; i<Materials(); ++i)
		{
			FEMaterial* pmat = GetMaterial(i);

			// store the type string
			ar << pmat->GetTypeStr();

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
			FEMaterial* pmat = fecore_new<FEMaterial>(FEMATERIAL_ID, szmat, this);
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
//! \todo Serialize nonlinear constraints
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
		FECoreKernel& febio = FECoreKernel::GetInstance();

		// read nodal data
		int nn;
		ar >> nn;
		m.CreateNodes(nn);
		for (int i=0; i<nn; ++i) ar.read(&m.Node(i), sizeof(FENode), 1);

		// read domain data
		int ND;
		ar >> ND;
		for (int i=0; i<ND; ++i)
		{
			int nmat;
			ar >> nmat;
			FEMaterial* pm = GetMaterial(nmat);
			assert(pm);

			int ntype, ne;
			ar >> ntype >> ne;
			FEDomain* pd = febio.CreateDomain(ntype, &m, pm);
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
	FECoreKernel& febio = FECoreKernel::GetInstance();

	if (ar.IsSaving())
	{
		ar << SurfacePairInteractions();
		for (int i=0; i<SurfacePairInteractions(); ++i)
		{
			FESurfacePairInteraction* pci = SurfacePairInteraction(i);

			// store the type string
			ar << pci->GetTypeStr();

			pci->Serialize(ar);
		}
	}
	else
	{
		int numci;
		ar >> numci;

		char szci[256] = {0};
		for (int i=0; i<numci; ++i)
		{
			// get the interface type
			ar >> szci;

			// create a new interface
			FESurfacePairInteraction* pci = fecore_new<FESurfacePairInteraction>(FESURFACEPAIRINTERACTION_ID, szci, this);

			// serialize interface data from archive
			pci->Serialize(ar);

			// add interface to list
			AddSurfacePairInteraction(pci);
		}	
	}
}

//-----------------------------------------------------------------------------
//! \todo Do we need to store the m_bActive flag of the boundary conditions?
void FEBioModel::SerializeBoundaryData(DumpFile& ar)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();

	if (ar.IsSaving())
	{
		// displacements
		ar << (int) m_DC.size();
		for (int i=0; i<(int) m_DC.size(); ++i) 
		{
			FEPrescribedBC& dc = *m_DC[i];
			ar << dc.GetID() << dc.IsActive();
			ar << dc.bc << dc.lc << dc.node << dc.s << dc.br << dc.r; //DSR
		}

		// nodal loads
		ar << (int) m_FC.size();
		for (int i=0; i<(int) m_FC.size(); ++i)
		{
			FENodalForce& fc = *m_FC[i];
			ar << fc.GetID() << fc.IsActive();
			ar << fc.bc << fc.lc << fc.node << fc.s;
		}

		// surface loads
		ar << (int) m_SL.size();
		for (int i=0; i<(int) m_SL.size(); ++i)
		{
			FESurfaceLoad* psl = m_SL[i];

			// get the surface
			FESurface& s = psl->Surface();
			s.Serialize(ar);

			// save the load data
			ar << psl->GetTypeStr();
			ar << psl->GetID() << psl->IsActive();
			psl->Serialize(ar);
		}

		// fixed rigid body dofs
		ar << m_RBC.size();
		for (int i=0; i<(int) m_RBC.size(); ++i)
		{
			FERigidBodyFixedBC& bc = *m_RBC[i];
			ar << bc.GetID() << bc.IsActive();
			ar << bc.bc << bc.id;
		}

		// rigid body displacements
		ar << m_RDC.size();
		for (int i=0; i<(int) m_RDC.size(); ++i)
		{
			FERigidBodyDisplacement& dc = *m_RDC[i];
			ar << dc.GetID() << dc.IsActive();
			ar << dc.bc << dc.id << dc.lc << dc.sf;
		}

		// rigid body forces
		ar << m_RFC.size();
		for (int i=0; i<(int) m_RFC.size(); ++i)
		{
			FERigidBodyForce& fc = *m_RFC[i];
			ar << fc.GetID() << fc.IsActive();
			ar << fc.bc << fc.id << fc.lc << fc.sf;
		}

		// rigid nodes
		ar << m_RN.size();
		for (int i=0; i<(int) m_RN.size(); ++i)
		{
			FERigidNode& rn = *m_RN[i];
			ar << rn.GetID() << rn.IsActive();
			ar << rn.nid << rn.rid;
		}

		// linear constraints
		ar << (int) m_LinC.size();
		list<FELinearConstraint>::iterator it = m_LinC.begin();
		for (int i=0; i<(int) m_LinC.size(); ++i, ++it) it->Serialize(ar);

		ar << m_LCT;

		// aug lag linear constraints
/*		n = (int) m_LCSet.size();
		ar << n;
		if (m_LCSet.empty() == false)
		{
			for (i=0; i<n; ++i) m_LCSet[i]->Serialize(ar);
		}
*/
		int n = m_NLC.size();
		ar << n;
		if (n) 
		{
			for (int i=0; i<n; ++i) 
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
		for (int i=0; i<n; ++i) 
		{
			FEPrescribedBC* pdc = new FEPrescribedBC(this);
			ar >> nid >> bactive;
			ar >> pdc->bc >> pdc->lc >> pdc->node >> pdc->s >> pdc->br >> pdc->r; // GAA
			pdc->SetID(nid);
			if (bactive) pdc->Activate(); else pdc->Deactivate();
			m_DC.push_back(pdc);
		}
		
		// nodal loads
		ar >> n;
		m_FC.clear();
		for (int i=0; i<n; ++i)
		{
			FENodalForce* pfc = new FENodalForce(this);
			ar >> nid >> bactive;
			ar >> pfc->bc >> pfc->lc >> pfc->node >> pfc->s;
			pfc->SetID(nid);
			if (bactive) pfc->Activate(); else pfc->Deactivate();
			m_FC.push_back(pfc);
		}

		// surface loads
		ar >> n;
		m_SL.clear();
		for (int i=0; i<n; ++i)
		{
			// create a new surface
			FESurface* psurf = new FESurface(&m_mesh);
			psurf->Serialize(ar);

			// read load data
			char sztype[256] = {0};
			ar >> sztype;
			FESurfaceLoad* ps = fecore_new<FESurfaceLoad>(FESURFACELOAD_ID, sztype, this);
			assert(ps);
			ps->SetSurface(psurf);

			ar >> nid >> bactive;
			ps->SetID(nid);
			if (bactive) ps->Activate(); else ps->Deactivate();

			ps->Serialize(ar);
			m_SL.push_back(ps);
		}

		// fixed rigid body dofs
		ar >> n;
		m_RBC.clear();
		for (int i=0; i<n; ++i)
		{
			FERigidBodyFixedBC* pbc = new FERigidBodyFixedBC(this);
			ar >> nid >> bactive;
			ar >> pbc->bc >> pbc->id;
			pbc->SetID(nid);
			if (bactive) pbc->Activate(); else pbc->Deactivate();
			m_RBC.push_back(pbc);
		}

		// rigid body displacements
		ar >> n;
		m_RDC.clear();
		for (int i=0; i<n; ++i)
		{
			FERigidBodyDisplacement* pdc = new FERigidBodyDisplacement(this);
			ar >> nid >> bactive;
			ar >> pdc->bc >> pdc->id >> pdc->lc >> pdc->sf;
			pdc->SetID(nid);
			if (bactive) pdc->Activate(); else pdc->Deactivate();
			m_RDC.push_back(pdc);
		}

		// rigid body forces
		ar >> n;
		m_RFC.clear();
		for (int i=0; i<n; ++i)
		{
			FERigidBodyForce* pfc = new FERigidBodyForce(this);
			ar >> nid >> bactive;
			ar >> pfc->bc >> pfc->id >> pfc->lc >> pfc->sf;
			pfc->SetID(nid);
			if (bactive) pfc->Activate(); else pfc->Deactivate();
			m_RFC.push_back(pfc);
		}

		// rigid nodes
		ar >> n;
		m_RN.clear();
		for (int i=0; i<n; ++i)
		{
			FERigidNode* prn = new FERigidNode(this);
			ar >> nid >> bactive;
			ar >> prn->nid >> prn->rid;
			prn->SetID(nid);
			if (bactive) prn->Activate(); else prn->Deactivate();
			m_RN.push_back(prn);
		}

		// linear constraints
		ar >> n;
		FELinearConstraint LC;
		for (int i=0; i<n; ++i)
		{
			LC.Serialize(ar);
			m_LinC.push_back(LC);
		}

		ar >> m_LCT;

		// reset the pointer table
		int nlin = m_LinC.size();
		m_LCA.resize(nlin);
		list<FELinearConstraint>::iterator ic = m_LinC.begin();
		for (int i=0; i<nlin; ++i, ++ic) m_LCA[i] = &(*ic);

		// aug lag linear constraints
		ar >> n;
//		int ntype;
		m_NLC.clear();
		for (int i=0; i<n; ++i)
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
		if (dynamic_cast<FEBioPlotFile*>(m_plot)) npltfmt = 2;
		assert(npltfmt != 0);
		ar << npltfmt;

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
			if (dynamic_cast<ObjectDataRecord*>(pd)) ntype = FE_DATA_RB;
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
			case FE_DATA_RB  : pd = new ObjectDataRecord(this, 0); break;
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
	// Open the logfile
	if (!felog.is_valid()) 
	{
		if (felog.open(m_szlog) == false)
		{
			felog.printbox("FATAL ERROR", "Failed creating log file");
			return false;
		}

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) felog.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Logfile::MODE m = felog.SetMode(Logfile::FILE_ONLY);
		Hello();
		felog.SetMode(m);
	}

	// initialize model data
	FEModel::Init();

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) 
		{
			m_plot = new FEBioPlotFile(*this);
		}

		if (m_plot->Open(*this, m_szplot) == false)
		{
			felog.printf("ERROR : Failed creating PLOT database\n");
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
	DoCallback(CB_MAJOR_ITERS);

	// Alright, all initialization is done, so let's get busy !
	return true;
}

//-----------------------------------------------------------------------------
//! This function resets the FEM data so that a new run can be done.
//! This routine is called from the optimization routine.

bool FEBioModel::Reset()
{
	// Reset model data
	FEModel::Reset();

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) 
		{
			m_plot = new FEBioPlotFile(*this);
		}

		if (m_plot->Open(*this, m_szplot) == false)
		{
			felog.printf("ERROR : Failed creating PLOT database\n");
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
	DoCallback(CB_MAJOR_ITERS);

	// All data is reset successfully
	return true;
}

//=============================================================================
//                               S O L V E
//=============================================================================

//-----------------------------------------------------------------------------
//! This is the main solve method. This function loops over all analysis steps
//! and solves each one in turn. 
//! \sa FEAnalysis

bool FEBioModel::Solve()
{
	// echo fem data to the logfile
	// we do this here (and not e.g. directly after input)
	// since the data can be changed after input, which is the case,
	// for instance, in the parameter optimization module
	if (m_becho) echo_input(*this);

	// start the total time tracker
	m_TotalTime.start();

	// solve the FE model
	bool bconv = FEModel::Solve();

	// stop total time tracker
	m_TotalTime.stop();

	// get and print elapsed time
	char sztime[64];
	m_TotalTime.time_str(sztime);
	felog.printf("\n Elapsed time : %s\n\n", sztime);

	if (bconv)
	{
		felog.printf("\n N O R M A L   T E R M I N A T I O N\n\n");
	}
	else
	{
		felog.printf("\n E R R O R   T E R M I N A T I O N\n\n");
	}

	// flush the log file
	felog.flush();

	// We're done !
	return bconv;
}

//-----------------------------------------------------------------------------
//!  This routine reads a binary archive that stores a restart point and prepares
//!  the FEM data to be restarted from this point
//!	\param[in] szfile name of the file

bool FEBioModel::Restart(const char* szfile)
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

		// see if user redefined restart file name
		if (file.m_szdmp[0]) SetDumpFilename(file.m_szdmp);
	}

	// Open the log file for appending
	if (felog.append(m_szlog) == false)
	{
		printf("WARNING: Could not reopen log file. A new log file is created\n");
		felog.open(m_szlog);
		return false;
	}

	// inform the user from where the problem is restarted
	felog.printbox(" - R E S T A R T -", "Restarting from time %lg.\n", m_ftime);

	return true;
}
