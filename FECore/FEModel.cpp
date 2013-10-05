#include "stdafx.h"
#include "FEModel.h"
#include <string>
using namespace std;

//-----------------------------------------------------------------------------
// Global Constants Data 
// m_Const and m_SD need a definitions, since static
map<std::string, double> FEModel::m_Const;
vector<FESoluteData*> FEModel::m_SD;
vector<FESBMData*> FEModel::m_SBM;

//-----------------------------------------------------------------------------
FEModel::FEModel(void)
{
	// --- Analysis Data ---
	m_pStep = 0;
	m_nStep = -1;
	m_nplane_strain = -1;	// don't use plain strain mode
	m_ftime = 0;
	m_ftime0 = 0;
	m_bwopt = 0;
}

//-----------------------------------------------------------------------------
//! Delete all dynamically allocated data
FEModel::~FEModel(void)
{
	size_t i;
	for (i=0; i<m_Step.size(); ++i) delete m_Step[i]; m_Step.clear();
	for (i=0; i<m_CI.size  (); ++i) delete m_CI [i] ; m_CI.clear  ();
	for (i=0; i<m_MAT.size (); ++i) delete m_MAT[i] ; m_MAT.clear ();
	for (i=0; i<m_LC.size  (); ++i) delete m_LC [i] ; m_LC.clear  ();
	for (i=0; i<m_BL.size  (); ++i) delete m_BL [i] ; m_BL.clear  ();
	for (i=0; i<m_DC.size  (); ++i) delete m_DC [i] ; m_DC.clear  ();
	for (i=0; i<m_FC.size  (); ++i) delete m_FC [i] ; m_FC.clear  ();
	for (i=0; i<m_SL.size  (); ++i) delete m_SL [i] ; m_SL.clear  ();
	for (i=0; i<m_RDC.size (); ++i) delete m_RDC[i] ; m_RDC.clear ();
	for (i=0; i<m_RFC.size (); ++i) delete m_RFC[i] ; m_RFC.clear ();
	for (i=0; i<m_RN.size  (); ++i) delete m_RN [i] ; m_RN.clear  ();
	for (i=0; i<m_NLC.size (); ++i) delete m_NLC[i] ; m_NLC.clear ();
	for (i=0; i<m_Obj.size (); ++i) delete m_Obj[i] ; m_Obj.clear ();
}

//-----------------------------------------------------------------------------
void FEModel::ClearBCs()
{
	for (size_t i=0; i<m_DC.size  (); ++i) delete m_DC[i];
	m_DC.clear();
}

//=============================================================================
//    P A R A M E T E R   F U N C T I O N S
//=============================================================================

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

double* FEModel::FindParameter(const char* szparam)
{
	char szname[256];
	strcpy(szname, szparam);

	// get the material and parameter name
	char* ch = strchr((char*)szname, '.');
	if (ch == 0) return 0;
	*ch = 0;
	const char* szmat = szname;
	const char* szvar = ch+1;

	// find the material with the same name
	FEMaterial* pmat = 0;
	int nmat = -1;
	for (int i=0; i<Materials(); ++i)
	{
		pmat = GetMaterial(i);
		nmat = i;
		if (strcmp(szmat, pmat->GetName()) == 0) break;
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
		index = atoi(szarg);
	}

	// find the material parameter
	ParamString sz(szvar);
	FEParam* pp = pmat->GetParameter(sz);
	if (pp) return pp->pvalue<double>(index);

	// the rigid bodies are dealt with differently
	int nrb = m_Obj.size();
	for (int i=0; i<nrb; ++i)
	{
		FEObject& ob = *m_Obj[i];
		if (ob.GetMaterialID() == nmat)
		{
			FEParam* pp = ob.GetParameter(sz);
			if (pp) return pp->pvalue<double>(index);
		}
	}

	// oh, oh, we didn't find it
	return 0;
}

//-----------------------------------------------------------------------------
//! Evaluate a parameter list
void FEModel::EvaluateParameterList(FEParameterList &pl)
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
void FEModel::EvaluateMaterialParameters(FEMaterial* pm)
{
	// evaluate the materials' parameter list
	EvaluateParameterList(pm->GetParameterList());

	// evaluate the material properties
	int N = pm->Properties();
	for (int i=0; i<N; ++i)
	{
		FEMaterial* pmi = pm->GetProperty(i);
		if (pmi) EvaluateMaterialParameters(pmi);
	}
}

//-----------------------------------------------------------------------------
// This function adds a callback routine
//
void FEModel::AddCallback(FEBIO_CB_FNC pcb, unsigned int nwhen, void *pd)
{
	FEBIO_CALLBACK cb;
	cb.m_pcb = pcb;
	cb.m_pd = pd;
	cb.m_nwhen = nwhen;

	m_pcb.push_back(cb);
}

//-----------------------------------------------------------------------------
// Call the callback function if there is one defined
//
void FEModel::DoCallback(unsigned int nevent)
{
	list<FEBIO_CALLBACK>::iterator it = m_pcb.begin();
	for (int i=0; i<(int) m_pcb.size(); ++i, ++it)
	{
		// call the callback function
		if (it->m_nwhen & nevent) (it->m_pcb)(this, it->m_pd);
	}
}

//-----------------------------------------------------------------------------
void FEModel::SetGlobalConstant(const string& s, double v)
{
	m_Const[s] = v;
	return;
}

//-----------------------------------------------------------------------------
double FEModel::GetGlobalConstant(const string& s)
{
	return (m_Const.count(s) ? m_Const.find(s)->second : 0);
}

//-----------------------------------------------------------------------------
void FEModel::SetSD(FESoluteData* psd)
{
	m_SD.push_back(psd);
}

//-----------------------------------------------------------------------------
FESoluteData* FEModel::FindSD(int nid)
{
	int i;
	for (i=0; i<(int) m_SD.size(); ++i) if (m_SD[i]->m_nID == nid) return m_SD[i];
	
	return 0;
}

//-----------------------------------------------------------------------------
void FEModel::SetSBM(FESBMData* psd)
{
	m_SBM.push_back(psd);
}

//-----------------------------------------------------------------------------
FESBMData* FEModel::FindSBM(int nid)
{
	int i;
	for (i=0; i<(int) m_SBM.size(); ++i) if (m_SBM[i]->m_nID == nid) return m_SBM[i];
	
	return 0;
}

//-----------------------------------------------------------------------------
//! Set the title of the model
void FEModel::SetTitle(const char* sz)
{ 
	strcpy(m_sztitle, sz); 
}

//-----------------------------------------------------------------------------
//! Return the title of the model
const char* FEModel::GetTitle()
{ 
	return m_sztitle; 
}

//-----------------------------------------------------------------------------
//! Find a BC based on its ID. This is needed for restarts.
FEBoundaryCondition* FEModel::FindBC(int nid)
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
