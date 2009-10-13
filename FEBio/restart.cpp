// restart module
#include "stdafx.h"
#include "fem.h"
#include "FERigid.h"
#include "FERestartImport.h"
#include "FEFacet2FacetSliding.h"
#include "FESlidingInterface2.h"
#include "log.h"

//-----------------------------------------------------------------------------
//!  This routine reads a binary archive that stores a restart point and prepares
//!  the FEM data to be restarted from this point
//!	\param[in] szfile name of the file

bool FEM::Restart(const char* szfile)
{
	// Open the restart file
	FERestartImport file;
	if (file.Load(*this, szfile) == false)
	{
		char szerr[256];
		file.GetErrorMessage(szerr);
		fprintf(stderr, "%s", szerr);
		return false;
	}

	// get the logfile
	Logfile& log = GetLogfile();

	// Open the log file for appending
	if (log.append(m_szlog) == false)
	{
		printf("WARNING: Could not reopen log file. A new log file is created\n");
		log.open(m_szlog);
		return false;
	}

	// Open the plot file for appending
	if (m_szplot)
	{
		if (m_plot.Append(*this, m_szplot) == false)
		{
			printf("FATAL ERROR: Failed reopening plot database %s\n", m_szplot);
			return false;
		}
	}

	// inform the user from where the problem is restarted
	log.printbox(" - R E S T A R T -", "Restarting from time %lg.\n", m_ftime);

	return true;
}

//-----------------------------------------------------------------------------
//!  Reads or writes the current state to/from a binary file
//!  This is used to restart the solution from a saved position
//!  or to create a restart point.
//!  A version number is written to file to make sure the same
//!  format is used for reading and writing.
//! \param[in] ar the archive to which the data is serialized
//! \sa Archive

bool FEM::Serialize(Archive& ar)
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

void FEM::SerializeLoadData(Archive& ar)
{
	if (ar.IsSaving())
	{
		// load curve data
		ar << LoadCurves();
		for (int i=0; i<LoadCurves(); ++i) GetLoadCurve(i)->Serialize(ar);
	}
	else
	{
		// loadcurve data
		int nlc;
		ar >> nlc;
		for (int i=0; i<nlc; ++i)
		{
			FELoadCurve* plc = new FELoadCurve;
			plc->Serialize(ar);
			AddLoadCurve(plc);
		}
	}
}

//-----------------------------------------------------------------------------

void FEM::SerializeAnalysisData(Archive &ar)
{
	if (ar.IsSaving())
	{
		// analysis steps
		int i;
		ar << (int) m_Step.size();
		for (i=0; i<(int) m_Step.size(); ++i) m_Step[i].Serialize(ar);
		ar << m_nStep;
		ar << m_ftime;
		ar << m_nhex8;

		// body forces
		ar.write(m_BF  ,sizeof(FE_BODYFORCE), 1);
		ar.write(m_BF+1,sizeof(FE_BODYFORCE), 1);
		ar.write(m_BF+2,sizeof(FE_BODYFORCE), 1);

		ar << m_acc;

		// direct solver data
		ar << m_nsolver;
		ar << m_neq;
		ar << m_bwopt;
	}
	else
	{
		m_Step.clear();
		// analysis steps
		int nsteps, i;
		ar >> nsteps;
		for (i=0; i<nsteps; ++i)
		{
			FEAnalysis* pstep = new FEAnalysis(*this);
			pstep->Serialize(ar);
			m_Step.add(pstep);
		}
		ar >> m_nStep;
		ar >> m_ftime;
		ar >> m_nhex8;

		// body forces
		ar.read(m_BF  ,sizeof(FE_BODYFORCE), 1);
		ar.read(m_BF+1,sizeof(FE_BODYFORCE), 1);
		ar.read(m_BF+2,sizeof(FE_BODYFORCE), 1);

		ar >> m_acc;

		// direct solver data
		ar >> m_nsolver;
		ar >> m_neq;
		ar >> m_bwopt;
	}
}

//-----------------------------------------------------------------------------

void FEM::SerializeMaterials(Archive& ar)
{
	if (ar.IsSaving())
	{
		ar << Materials();
		for (int i=0; i<Materials(); ++i)
		{
			FEMaterial* pmat = GetMaterial(i);

			// store the type string
			ar << pmat->GetTypeString();

			// store the name
			ar << pmat->GetName();

			// store all parameters
			auto_ptr<FEParameterList> pl(pmat->GetParameterList());
			int n = pl->Parameters();
			ar << n;
			list<FEParam>::iterator it = pl->first();
			for (int j=0; j<n; ++j, ++it)
			{
				// store the value
				switch (it->m_itype)
				{
				case FE_PARAM_INT    : ar << it->value<int   >(); break;
				case FE_PARAM_BOOL   : ar << it->value<bool  >(); break;
				case FE_PARAM_DOUBLE : ar << it->value<double>(); break;
				case FE_PARAM_DOUBLEV: { for (int k=0; k<it->m_ndim; ++k) ar << it->pvalue<double>()[k]; } break;
				case FE_PARAM_INTV   : { for (int k=0; k<it->m_ndim; ++k) ar << it->pvalue<int   >()[k]; } break;
				default:
					assert(false);
				}

				// store parameter data
				ar << it->m_nlc;
			}

			// not all parameters can be serialized through the parameter lists
			// so we have to save those parameters the hard way

			if (dynamic_cast<FETransverselyIsotropic*>(pmat))
			{
				FETransverselyIsotropic* pm = dynamic_cast<FETransverselyIsotropic*>(pmat);
				ar << pm->lcna;
				ar << pm->m_ascl;
				ar << pm->ca0;
				ar << pm->beta;
				ar << pm->l0;
				ar << pm->refl;
			}

			// TODO: do we really need to store this data?
			if (dynamic_cast<FERigid*>(pmat))
			{
				FERigid* pm = dynamic_cast<FERigid*>(pmat);
				ar.write(pm->m_bc, sizeof(int), 6);
				ar.write(pm->m_fc, sizeof(int), 6);
				ar.write(pm->m_fs, sizeof(double), 6);
			}
		}
	}
	else
	{
		int nmat;
		char szmat[256] = {0}, szvar[256] = {0};
		ar >> nmat;
		for (int i=0; i<nmat; ++i)
		{
			// read the type string
			ar >> szmat;

			// create a material
			FEMaterial* pmat = FEMaterialFactory::CreateMaterial(szmat);
			assert(pmat);
			AddMaterial(pmat);

			// read the name
			ar >> szmat;
			pmat->SetName(szmat);

			// read all parameters
			FEParameterList* pl = pmat->GetParameterList();
			AddParameterList(pl);
			int n = 0;
			ar >> n;
			assert(n == pl->Parameters());
			list<FEParam>::iterator it = pl->first();
			for (int j=0; j<n; ++j, ++it)
			{
				// read the value
				switch (it->m_itype)
				{
				case FE_PARAM_INT    : ar >> it->value<int   >(); break;
				case FE_PARAM_BOOL   : ar >> it->value<bool  >(); break;
				case FE_PARAM_DOUBLE : ar >> it->value<double>(); break;
				case FE_PARAM_DOUBLEV: { for (int k=0; k<it->m_ndim; ++k) ar >> it->pvalue<double>()[k]; } break;
				case FE_PARAM_INTV   : { for (int k=0; k<it->m_ndim; ++k) ar >> it->pvalue<int   >()[k]; } break;
				default:
					assert(false);
				}

				// read parameter data
				ar >> it->m_nlc;
			}

			// not all parameters can be serialized through the parameter lists
			// so we have to save those parameters the hard way

			if (dynamic_cast<FETransverselyIsotropic*>(pmat))
			{
				FETransverselyIsotropic* pm = dynamic_cast<FETransverselyIsotropic*>(pmat);
				ar >> pm->lcna;
				ar >> pm->m_ascl;
				ar >> pm->ca0;
				ar >> pm->beta;
				ar >> pm->l0;
				ar >> pm->refl;

				if (pm->lcna >= 0) pm->m_plc = GetLoadCurve(pm->lcna);
			}

			if (dynamic_cast<FERigid*>(pmat))
			{
				FERigid* pm = dynamic_cast<FERigid*>(pmat);
				ar.read(pm->m_bc, sizeof(int), 6);
				ar.read(pm->m_fc, sizeof(int), 6);
				ar.read(pm->m_fs, sizeof(double), 6);
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FEM::SerializeGeometry(Archive &ar)
{
	// serialize the mesh first 
	m_mesh.Serialize(ar);

	// serialize the other geometry data
	if (ar.IsSaving())
	{
		int i, j, n;

		// write solid element state data
		for (i=0; i<m_mesh.SolidElements(); ++i)
		{
			FESolidElement& el = m_mesh.SolidElement(i);
			for (j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}

		// write shell element state data
		for (i=0; i<m_mesh.ShellElements(); ++i)
		{
			FEShellElement& el = m_mesh.ShellElement(i);
			for (j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}

		// surface elements
		n = m_psurf->Elements();
		ar << n;
		for (i=0; i<n; ++i)
		{
			FESurfaceElement& el = m_psurf->Element(i);
			ar << el.Type();
			ar << el.GetMatID();
			ar << el.m_nID;
			ar << el.m_nrigid;
			ar << el.m_node;
			ar << el.m_lnode;
		}

		// rigid bodies
		ar << m_nreq << m_nrm << m_nrb;
		for (i=0; i<m_nrb; ++i) m_RB[i].Serialize(ar);

		// rigid joints
		ar << m_nrj;
		for (i=0; i<m_nrj; ++i) m_RJ[i].Serialize(ar);	

		// discrete elements
		ar << m_DE.size();
		for (i=0; i<m_DE.size(); ++i)
		{
			FE_DISCRETE_ELEMENT& de = m_DE[i];
			ar << de.n1 << de.n2;
			ar << de.E;
		}
	}
	else
	{
		int i, j, n, m, mat;

		// read solid element state data
		for (i=0; i<m_mesh.SolidElements(); ++i)
		{
			FESolidElement& el = m_mesh.SolidElement(i);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}

		// read shell element state data
		for (i=0; i<m_mesh.ShellElements(); ++i)
		{
			FEShellElement& el = m_mesh.ShellElement(i);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}

		// surface elements
		ar >> n;
		if (n) 
		{
			m_psurf->Create(n);

			for (i=0; i<n; ++i)
			{
				FESurfaceElement& el = m_psurf->Element(i);
				ar >> m;
				el.SetType(m);

				ar >> mat; el.SetMatID(mat);
				ar >> el.m_nID;
				ar >> el.m_nrigid;
				ar >> el.m_node;
				ar >> el.m_lnode;
			}

			// initialize surface data
			m_psurf->Init();
		}

		// rigid bodies
		ar >> m_nreq >> m_nrm >> m_nrb;
		if (m_nrb) m_RB.create(m_nrb);
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
			m_RJ.add(prj);
		}

		// discrete elements
		int nde;
		ar >> nde;
		if (nde > 0) m_DE.setsize(nde);
		for (i=0; i<nde; ++i)
		{
			FE_DISCRETE_ELEMENT& de = m_DE[i];
			ar >> de.n1 >> de.n2;
			ar >> de.E;
		}
	}
}

//-----------------------------------------------------------------------------

void FEM::SerializeContactData(Archive &ar)
{
	if (ar.IsSaving())
	{
		ar << m_bcontact;
		ar << ContactInterfaces();
		for (int i=0; i<ContactInterfaces(); ++i)
		{
			ar << m_CI[i].Type();
			m_CI[i].Serialize(ar);
		}
	}
	else
	{
		int numci, ntype;
		ar >> m_bcontact;
		ar >> numci;
		for (int i=0; i<numci; ++i)
		{
			FEContactInterface* ps;

			// get the interface type
			ar >> ntype;

			// create a new interface
			switch (ntype)
			{
			case FE_CONTACT_SLIDING  : ps = new FESlidingInterface(this); break;
			case FE_FACET2FACET_SLIDING: ps = new FEFacet2FacetSliding(this); break;
			case FE_CONTACT_TIED     : ps = new FETiedInterface(this); break;
			case FE_CONTACT_RIGIDWALL: ps = new FERigidWallInterface(this); break;
			case FE_CONTACT_SLIDING2 : ps = new FESlidingInterface2(this); break;
			default:
				assert(false);
			}
				
			// serialize interface data from archive
			ps->Serialize(ar);

			// add interface to list
			m_CI.add(ps);
		}	
	}
}

//-----------------------------------------------------------------------------
// TODO: do we need to store the m_bActive flag of the boundary conditions?
//
void FEM::SerializeBoundaryData(Archive& ar)
{
	int i;

	if (ar.IsSaving())
	{
		// displacements
		ar << m_DC.size();
		for (i=0; i<m_DC.size(); ++i) 
		{
			FENodalDisplacement& dc = m_DC[i];
			ar << dc.bc << dc.lc << dc.node << dc.s;
		}

		// nodal loads
		ar << m_FC.size();
		for (i=0; i<m_FC.size(); ++i)
		{
			FENodalForce& fc = m_FC[i];
			ar << fc.bc << fc.lc << fc.node << fc.s;
		}

		// pressure forces
		ar << m_PC.size();
		for (i=0; i<m_PC.size(); ++i)
		{
			FEPressureLoad& pc = m_PC[i];
			ar << pc.blinear << pc.face << pc.lc;
			ar << pc.s[0] << pc.s[1] << pc.s[2] << pc.s[3];
			ar << pc.bc;
		}

		// rigid body displacements
		ar << m_RDC.size();
		for (i=0; i<m_RDC.size(); ++i)
		{
			FERigidBodyDisplacement& dc = m_RDC[i];
			ar << dc.bc << dc.id << dc.lc << dc.sf;
		}

		// rigid body forces
		ar << m_RFC.size();
		for (i=0; i<m_RFC.size(); ++i)
		{
			FERigidBodyForce& fc = m_RFC[i];
			ar << fc.bc << fc.id << fc.lc << fc.sf;
		}

		// rigid nodes
		ar << m_RN.size();
		for (i=0; i<m_RN.size(); ++i)
		{
			FERigidNode& rn = m_RN[i];
			ar << rn.nid << rn.rid;
		}

		// linear constraints
		ar << (int) m_LinC.size();
		list<FELinearConstraint>::iterator it = m_LinC.begin();
		for (i=0; i<(int) m_LinC.size(); ++i, ++it) it->Serialize(ar);

		ar << m_LCT;

		// aug lag linear constraints
		ar << (int) m_LCSet.size();
		list<FELinearConstraintSet*>::iterator ic = m_LCSet.begin();
		for (i=0; i< (int) m_LCSet.size(); ++i, ++ic) (*ic)->Serialize(ar);

		// aug lag traction constraints
		ar << (int) m_RCSet.size();
		list<FETractionConstraintSet*>::iterator tc = m_RCSet.begin();
		for (i=0; i< (int) m_RCSet.size(); ++i, ++tc) (*tc)->Serialize(ar);
	}
	else
	{
		int n;
		// displacements
		ar >> n;
		if (n) m_DC.create(n);
		for (i=0; i<n; ++i) 
		{
			FENodalDisplacement& dc = m_DC[i];
			ar >> dc.bc >> dc.lc >> dc.node >> dc.s;
		}
		
		// nodal loads
		ar >> n;
		if (n) m_FC.create(n);
		for (i=0; i<n; ++i)
		{
			FENodalForce& fc = m_FC[i];
			ar >> fc.bc >> fc.lc >> fc.node >> fc.s;
		}

		// pressure forces
		ar >> n;
		if (n) m_PC.create(n);
		for (i=0; i<n; ++i)
		{
			FEPressureLoad& pc = m_PC[i];
			ar >> pc.blinear >> pc.face >> pc.lc;
			ar >> pc.s[0] >> pc.s[1] >> pc.s[2] >> pc.s[3];
			ar >> pc.bc;
		}

		// rigid body displacements
		ar >> n;
		if (n) m_RDC.create(n);
		for (i=0; i<n; ++i)
		{
			FERigidBodyDisplacement& dc = m_RDC[i];
			ar >> dc.bc >> dc.id >> dc.lc >> dc.sf;
		}

		// rigid body forces
		ar >> n;
		if (n) m_RFC.create(n);
		for (i=0; i<n; ++i)
		{
			FERigidBodyForce& fc = m_RFC[i];
			ar >> fc.bc >> fc.id >> fc.lc >> fc.sf;
		}

		// rigid nodes
		ar >> n;
		if (n) m_RN.create(n);
		for (i=0; i<n; ++i)
		{
			FERigidNode& rn = m_RN[i];
			ar >> rn.nid >> rn.rid;
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
		m_LCA.create(nlin);
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

		// aug lag linear constraints
		ar >> n;
		for (i=0; i<n; ++i)
		{
			FETractionConstraintSet* prc = new FETractionConstraintSet(this);
			prc->Serialize(ar);
			m_RCSet.push_back(prc);
		}
	}
}

//-----------------------------------------------------------------------------
// TODO: serialize data records

void FEM::SerializeIOData(Archive &ar)
{
	if (ar.IsSaving())
	{
		// file names
		ar << m_szfile << m_szplot << m_szlog << m_szdump;
		ar << m_sztitle;

		// plot file
		int* n = m_plot.m_nfield;
		ar << n[0] << n[1] << n[2] << n[3] << n[4];
	}
	else
	{
		// file names
		ar >> m_szfile >> m_szplot >> m_szlog >> m_szdump;
		ar >> m_sztitle;

		// don't forget to call SetInputFilename so
		// that m_szfile_title gets initialized
		SetInputFilename(m_szfile);

		// plot file
		int* n = m_plot.m_nfield;
		ar >> n[0] >> n[1] >> n[2] >> n[3] >> n[4];
	}
}
