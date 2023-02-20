/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FECoreKernel.h"
#include "LinearSolver.h"
#include "Timer.h"
#include "FEModule.h"
#include <stdarg.h>
using namespace std;

//-----------------------------------------------------------------------------
FECoreKernel* FECoreKernel::m_pKernel = 0;

//-----------------------------------------------------------------------------
FECoreKernel& FECoreKernel::GetInstance()
{
	if (m_pKernel == 0) m_pKernel = new FECoreKernel;
	return *m_pKernel;
}

//-----------------------------------------------------------------------------
// This function is used by plugins to make sure that the plugin and the executable
// are using the same kernel
void FECoreKernel::SetInstance(FECoreKernel* pkernel)
{
	m_pKernel = pkernel;
}

//-----------------------------------------------------------------------------
const char* FECoreKernel::SuperClassString(unsigned int sid)
{
	FECoreKernel& fecore = GetInstance();
	if (fecore.m_sidMap.find(sid) != fecore.m_sidMap.end())
		return fecore.m_sidMap[sid];
	else
		return nullptr;
}

//-----------------------------------------------------------------------------
std::map<unsigned int, const char*>	FECoreKernel::GetSuperClassMap()
{
	FECoreKernel& fecore = GetInstance();
	return fecore.m_sidMap;
}

//-----------------------------------------------------------------------------

#define ADD_SUPER_CLASS(a) m_sidMap[a] = #a
FECoreKernel::FECoreKernel()
{
	m_activeModule = -1;
	m_alloc_id = 0;
	m_next_alloc_id = 1;
	m_nspec = -1;
	m_default_solver = nullptr;
	m_blockEvents = true;

	// build the super class ID table
	ADD_SUPER_CLASS(FEINVALID_ID);
	ADD_SUPER_CLASS(FEOBJECT_ID);
	ADD_SUPER_CLASS(FETASK_ID);
	ADD_SUPER_CLASS(FESOLVER_ID);
	ADD_SUPER_CLASS(FEMATERIAL_ID);
	ADD_SUPER_CLASS(FEMATERIALPROP_ID);
	ADD_SUPER_CLASS(FEDISCRETEMATERIAL_ID);
	ADD_SUPER_CLASS(FELOAD_ID);
	ADD_SUPER_CLASS(FENLCONSTRAINT_ID);
	ADD_SUPER_CLASS(FEPLOTDATA_ID);
	ADD_SUPER_CLASS(FEANALYSIS_ID);
	ADD_SUPER_CLASS(FESURFACEINTERFACE_ID);
	ADD_SUPER_CLASS(FELOGNODEDATA_ID);
	ADD_SUPER_CLASS(FELOGFACEDATA_ID);
	ADD_SUPER_CLASS(FELOGELEMDATA_ID);
	ADD_SUPER_CLASS(FELOGOBJECTDATA_ID);
	ADD_SUPER_CLASS(FELOGDOMAINDATA_ID);
	ADD_SUPER_CLASS(FELOGNLCONSTRAINTDATA_ID);
	ADD_SUPER_CLASS(FELOGSURFACEDATA_ID);
	ADD_SUPER_CLASS(FELOGMODELDATA_ID);
	ADD_SUPER_CLASS(FEBC_ID);
	ADD_SUPER_CLASS(FEGLOBALDATA_ID);
	ADD_SUPER_CLASS(FECALLBACK_ID);
	ADD_SUPER_CLASS(FESOLIDDOMAIN_ID);
	ADD_SUPER_CLASS(FESHELLDOMAIN_ID);
	ADD_SUPER_CLASS(FEBEAMDOMAIN_ID);
	ADD_SUPER_CLASS(FEDISCRETEDOMAIN_ID);
	ADD_SUPER_CLASS(FEDOMAIN2D_ID);
	ADD_SUPER_CLASS(FESURFACE_ID);
	ADD_SUPER_CLASS(FEIC_ID);
	ADD_SUPER_CLASS(FEMESHDATAGENERATOR_ID);
	ADD_SUPER_CLASS(FELOADCONTROLLER_ID);
	ADD_SUPER_CLASS(FEMODEL_ID);
	ADD_SUPER_CLASS(FESCALARVALUATOR_ID);
	ADD_SUPER_CLASS(FEVEC3DVALUATOR_ID);
	ADD_SUPER_CLASS(FEMAT3DVALUATOR_ID);
	ADD_SUPER_CLASS(FEMAT3DSVALUATOR_ID);
	ADD_SUPER_CLASS(FEFUNCTION1D_ID);
	ADD_SUPER_CLASS(FELINEARSOLVER_ID);
	ADD_SUPER_CLASS(FEMESHADAPTOR_ID);
	ADD_SUPER_CLASS(FEMESHADAPTORCRITERION_ID);
	ADD_SUPER_CLASS(FENEWTONSTRATEGY_ID);
	ADD_SUPER_CLASS(FETIMECONTROLLER_ID);
	ADD_SUPER_CLASS(FEEIGENSOLVER_ID);
	ADD_SUPER_CLASS(FEDATARECORD_ID);
	ADD_SUPER_CLASS(FECLASS_ID);
}

//-----------------------------------------------------------------------------
// Generate a allocator ID
int FECoreKernel::GenerateAllocatorID()
{
	return m_next_alloc_id++;
}

//-----------------------------------------------------------------------------
FECoreFactory* FECoreKernel::SetDefaultSolverType(const char* sztype)
{
	FECoreFactory* fac = FindFactoryClass(FELINEARSOLVER_ID, sztype);
	if (fac) m_default_solver_type = sztype;
	return fac;
}

//-----------------------------------------------------------------------------
void FECoreKernel::SetDefaultSolver(FEClassDescriptor* linsolve)
{
	delete m_default_solver;
	m_default_solver = linsolve;

	if (linsolve)
	{
		m_default_solver_type = linsolve->ClassType();
	}
	else
	{
		m_default_solver_type.clear();
	}
}

//-----------------------------------------------------------------------------
//! get the linear solver type
const char* FECoreKernel::GetLinearSolverType() const
{
	return m_default_solver_type.c_str();
}

//-----------------------------------------------------------------------------
LinearSolver* FECoreKernel::CreateDefaultLinearSolver(FEModel* fem)
{
	if (m_default_solver == nullptr)
	{
		const char* sztype = m_default_solver_type.c_str();
		FECoreFactory* fac = FindFactoryClass(FELINEARSOLVER_ID, sztype);
		return (LinearSolver*)CreateInstance(fac, fem);
	}
	else
	{
		return (LinearSolver*)Create(FELINEARSOLVER_ID, fem, *m_default_solver);
	}
}

//-----------------------------------------------------------------------------
void FECoreKernel::RegisterFactory(FECoreFactory* ptf)
{
	unsigned int activeID = 0;
	if (m_activeModule != -1)
	{
		FEModule& activeModule = *m_modules[m_activeModule];
		activeID = activeModule.GetModuleID();
	}

	// see if the name already exists
	for (int i=0; i<m_Fac.size(); ++i)
	{
		FECoreFactory* pfi = m_Fac[i];

		if ((pfi->GetSuperClassID() == ptf->GetSuperClassID()) && 
			(strcmp(pfi->GetTypeStr(), ptf->GetTypeStr()) == 0))
		{
			// A feature with the same is already registered. 
			// We need to check the module to see if this would create an ambiguity
			unsigned int modId = pfi->GetModuleID();

			// If the same feature is defined in the active module,
			// then this feature will replace the existing one. 
			if ((modId == activeID) && (pfi->GetSpecID() == ptf->GetSpecID()))
			{
				fprintf(stderr, "WARNING: \"%s\" feature is redefined\n", ptf->GetTypeStr());
				m_Fac[i] = ptf;
				return;
			}
		}
	}

	// it doesn't so add it
	ptf->SetModuleID(activeID);
	ptf->SetAllocatorID(m_alloc_id);
	m_Fac.push_back(ptf);
}

//-----------------------------------------------------------------------------
bool FECoreKernel::UnregisterFactory(FECoreFactory* ptf)
{
	for (vector<FECoreFactory*>::iterator it = m_Fac.begin(); it != m_Fac.end(); ++it)
	{
		FECoreFactory* pfi = *it;
		if (pfi == ptf)
		{
			m_Fac.erase(it);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! unregister factories from allocator
void FECoreKernel::UnregisterFactories(int alloc_id)
{
	for (vector<FECoreFactory*>::iterator it = m_Fac.begin(); it != m_Fac.end();)
	{
		FECoreFactory* pfi = *it;
		if (pfi->GetAllocatorID() == alloc_id)
		{
			it = m_Fac.erase(it);
		}
		else ++it;
	}
}

//-----------------------------------------------------------------------------
//! set the current allocator ID
void FECoreKernel::SetAllocatorID(int alloc_id)
{
	m_alloc_id = alloc_id;
}

//-----------------------------------------------------------------------------
//! Create an object. An object is created by specifying the super-class id
//! and the type-string. 
FECoreBase* FECoreKernel::Create(int superClassID, const char* sztype, FEModel* pfem)
{
	FECoreFactory* fac = FindFactoryClass(superClassID, sztype);
	if (fac == nullptr) return nullptr;
	return CreateInstance(fac, pfem);
}

//-----------------------------------------------------------------------------
//! Create an object. An object is created by specifying the super-class id
//! and the type-string. 
FECoreBase* FECoreKernel::Create(const char* szbase, const char* sztype, FEModel* pfem)
{
	if (szbase == nullptr) return nullptr;
	if (sztype == nullptr) return nullptr;

	unsigned int activeID = 0;
	vector<int> moduleDepends;
	if (m_activeModule != -1)
	{
		FEModule& activeModule = *m_modules[m_activeModule];
		activeID = activeModule.GetModuleID();
		moduleDepends = activeModule.GetDependencies();
	}

	// first check active module
	if ((activeID > 0) || (activeID == 0))
	{
		std::vector<FECoreFactory*>::iterator pf;
		for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
		{
			FECoreFactory* pfac = *pf;
			if (strcmp(pfac->GetBaseClassName(), szbase) == 0) {

				// see if we can match module first
				unsigned int mid = pfac->GetModuleID();
				if ((mid == activeID) || (mid == 0))
				{
					// see if the type name matches
					if ((strcmp(pfac->GetTypeStr(), sztype) == 0))
					{
						// check the spec (TODO: What is this for?)
						int nspec = pfac->GetSpecID();
						if ((nspec == -1) || (m_nspec <= nspec))
						{
							return CreateInstance(pfac, pfem);
						}
					}
				}
			}
		}
	}

	// check dependencies in order in which they are defined
	std::vector<FECoreFactory*>::iterator pf;
	for (int i = 0; i < moduleDepends.size(); ++i)
	{
		unsigned modId = moduleDepends[i];
		for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
		{
			FECoreFactory* pfac = *pf;
			if (strcmp(pfac->GetBaseClassName(), szbase) == 0) {

				// see if we can match module first
				unsigned int mid = pfac->GetModuleID();
				if ((mid == 0) || (mid == modId))
				{
					// see if the type name matches
					if ((strcmp(pfac->GetTypeStr(), sztype) == 0))
					{
						// check the spec (TODO: What is this for?)
						int nspec = pfac->GetSpecID();
						if ((nspec == -1) || (m_nspec <= nspec))
						{
							return CreateInstance(pfac, pfem);
						}
					}
				}
			}
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! Create a specific class
FECoreBase* FECoreKernel::CreateClass(const char* szclassName, FEModel* fem)
{
	std::vector<FECoreFactory*>::iterator pf;
	for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		const char* szfacName = pfac->GetClassName();
		if (szfacName && (strcmp(szfacName, szclassName) == 0))
		{
			return CreateInstance(pfac, fem);
		}
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
//! Create a class from a class descriptor
FECoreBase* FECoreKernel::Create(int superClassID, FEModel* pfem, const FEClassDescriptor& cd)
{
	const FEClassDescriptor::ClassVariable* root = cd.Root();
	FECoreBase* pc = (FECoreBase*)Create(superClassID, root->m_type.c_str(), pfem);
	if (pc == nullptr) return nullptr;
	pc->SetParameters(cd);
	return pc;
}

//-----------------------------------------------------------------------------
bool FECoreKernel::IsModuleActive(int moduleID)
{
	if (moduleID <= 0) return true;
	if (m_activeModule < 0) return false;

	FEModule* mod = GetActiveModule();
	if (mod == nullptr) return false;

	if (mod->GetModuleID() == moduleID) return true;

	std::vector<int> deps = mod->GetDependencies();
	for (int i = 0; i < deps.size(); ++i)
	{
		if (deps[i] == moduleID) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
FECoreBase* FECoreKernel::CreateInstance(const FECoreFactory* fac, FEModel* fem)
{
	FECoreBase* pc = fac->CreateInstance(fem);
	if ((m_blockEvents == false) && pc && (m_createHandlers.empty() == false))
	{
		for (int i = 0; i < m_createHandlers.size(); ++i)
		{
			FECreateHandler* ph = m_createHandlers[i];
			if (ph && (IsModuleActive(ph->GetModuleID())))
			{
				ph->handle(pc);
			}
		}
	}
	return pc;
}

//-----------------------------------------------------------------------------
int FECoreKernel::Count(SUPER_CLASS_ID sid)
{
	int N = 0;
	std::vector<FECoreFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == sid) N++;
	}
	return N;
}

//-----------------------------------------------------------------------------
void FECoreKernel::List(SUPER_CLASS_ID sid)
{
  std::vector<FECoreFactory*>::iterator pf;
  for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
    {
      FECoreFactory* pfac = *pf;
      if (pfac->GetSuperClassID() == sid) fprintf(stdout, "%s\n", pfac->GetTypeStr());
    }
}

//-----------------------------------------------------------------------------
int FECoreKernel::FactoryClasses()
{
	return (int) m_Fac.size();
}

//-----------------------------------------------------------------------------
const FECoreFactory* FECoreKernel::GetFactoryClass(int i)
{
	if ((i < 0) || (i >= m_Fac.size())) return nullptr;
	else return m_Fac[i];
}

//-----------------------------------------------------------------------------
//! return a factory class
const FECoreFactory* FECoreKernel::GetFactoryClass(int classID, int i)
{
	int n = 0;
	for (int j = 0; j < m_Fac.size(); ++j)
	{
		FECoreFactory* fac = m_Fac[j];
		if (fac->GetSuperClassID() == classID)
		{
			if (i == n) return fac;
			n++;
		}
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
//! return a factory class
int FECoreKernel::GetFactoryIndex(int superClassId, const char* sztype)
{
	FEModule* mod = GetActiveModule();
	if (mod == nullptr) return -1;
	for (int j = 0; j < m_Fac.size(); ++j)
	{
		FECoreFactory* fac = m_Fac[j];
		int modId = fac->GetModuleID();

		// check the super-class first
		if (fac->GetSuperClassID() == superClassId)
		{
			// check the string name 
			if (strcmp(sztype, fac->GetTypeStr()) == 0)
			{
				// make sure it's part of the active module
				if ((modId == 0) || (mod->HasDependent(modId))) return j;
			}
		}
	}
	return -1;
}


//-----------------------------------------------------------------------------
FECoreFactory* FECoreKernel::FindFactoryClass(int superID, const char* sztype)
{
	if (sztype == nullptr) return nullptr;

	unsigned int activeID = 0;
	vector<int> moduleDepends;
	if (m_activeModule != -1)
	{
		FEModule& activeModule = *m_modules[m_activeModule];
		activeID = activeModule.GetModuleID();
		moduleDepends = activeModule.GetDependencies();
	}

	// first check active module
	if ((activeID > 0) || (activeID == 0))
	{
		std::vector<FECoreFactory*>::iterator pf;
		for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
		{
			FECoreFactory* pfac = *pf;
			if (pfac->GetSuperClassID() == superID) {

				// see if we can match module first
				unsigned int mid = pfac->GetModuleID();
				if ((mid == activeID) || (mid == 0))
				{
					// see if the type name matches
					if ((strcmp(pfac->GetTypeStr(), sztype) == 0))
					{
						// check the spec (TODO: What is this for?)
						int nspec = pfac->GetSpecID();
						if ((nspec == -1) || (m_nspec <= nspec))
						{
							return pfac;
						}
					}
				}
			}
		}
	}

	// check dependencies in order in which they are defined
	std::vector<FECoreFactory*>::iterator pf;
	for (int i = 0; i < moduleDepends.size(); ++i)
	{
		unsigned modId = moduleDepends[i];
		for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
		{
			FECoreFactory* pfac = *pf;
			if (pfac->GetSuperClassID() == superID) {

				// see if we can match module first
				unsigned int mid = pfac->GetModuleID();
				if ((mid == 0) || (mid == modId))
				{
					// see if the type name matches
					if ((strcmp(pfac->GetTypeStr(), sztype) == 0))
					{
						// check the spec (TODO: What is this for?)
						int nspec = pfac->GetSpecID();
						if ((nspec == -1) || (m_nspec <= nspec))
						{
							return pfac;
						}
					}
				}
			}
		}
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
//! set the active module
bool FECoreKernel::SetActiveModule(const char* szmod)
{
	// See if user want to deactivate modules
	if (szmod == 0)
	{
		m_activeModule = -1;
		return true;
	}

	// see if the module exists or not
	for (size_t i=0; i<m_modules.size(); ++i) 
	{
		FEModule& mi = *m_modules[i];
		if (strcmp(mi.GetName(), szmod) == 0)
		{
			m_activeModule = (int) i;
			return true;
		}
	}

	// couldn't find it
	m_activeModule = -1;
	return false;
}

//-----------------------------------------------------------------------------
bool FECoreKernel::SetActiveModule(int moduleId)
{
	if (moduleId < 0) return false;
	if (GetActiveModuleID() == moduleId) return true;

	// see if the module exists or not
	for (size_t i = 0; i < m_modules.size(); ++i)
	{
		FEModule& mi = *m_modules[i];
		if (mi.GetModuleID() == moduleId)
		{
			m_activeModule = (int)i;
			return true;
		}
	}

	// couldn't find it
	m_activeModule = -1;
	return false;
}

//-----------------------------------------------------------------------------
// return the active module's ID
int FECoreKernel::GetActiveModuleID()
{
	if (m_activeModule == -1) return -1;
	return m_modules[m_activeModule]->GetModuleID();
}

//-----------------------------------------------------------------------------
FEModule* FECoreKernel::GetActiveModule()
{
	if (m_activeModule == -1) return nullptr;
	return m_modules[m_activeModule];
}

//-----------------------------------------------------------------------------
//! count modules
int FECoreKernel::Modules() const
{
	return (int)m_modules.size();
}

//-----------------------------------------------------------------------------
//! create a module
bool FECoreKernel::CreateModule(const char* szmod, const char* description)
{
	FEModule* mod = new FEModule();
	return CreateModule(mod, szmod, description);
}

//-----------------------------------------------------------------------------
//! create a module
bool FECoreKernel::CreateModule(FEModule* pmodule, const char* szmod, const char* description)
{
	assert(pmodule);
	if (pmodule == nullptr) return false;

	m_activeModule = -1;
	if (szmod == 0) return false;

	// see if this module already exist
	if (SetActiveModule(szmod) == false)
	{
		// The module does not exist, so let's add it.
		unsigned int newID = (unsigned int)m_modules.size() + 1;

		pmodule->SetName(szmod);
		pmodule->SetID(newID);
		pmodule->SetDescription(description);
		m_modules.push_back(pmodule);

		// make this the active module
		m_activeModule = (int)m_modules.size() - 1;
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Get a module
const char* FECoreKernel::GetModuleName(int i) const
{
	if ((i<0) || (i >= m_modules.size())) return nullptr;
	return m_modules[i]->GetName();
}

const char* FECoreKernel::GetModuleDescription(int i) const
{
	if ((i < 0) || (i >= m_modules.size())) return nullptr;
	return m_modules[i]->GetDescription();
}

int FECoreKernel::GetModuleStatus(int i) const
{
	if ((i < 0) || (i >= m_modules.size())) return -1;
	return m_modules[i]->GetStatus();
}

//! Get a module
const char* FECoreKernel::GetModuleNameFromId(int id) const
{
	for (size_t n = 0; n < m_modules.size(); ++n)
	{
		const FEModule& mod = *m_modules[n];
		if (mod.GetModuleID() == id) return mod.GetName();
	}
	return 0;
}

//! Get a module's dependencies
vector<int> FECoreKernel::GetModuleDependencies(int i) const
{
	vector<int> md;
	if ((i >= 0) && (i < m_modules.size()))
	{
		md = m_modules[i]->GetDependencies();
	}
	return md;
}


//-----------------------------------------------------------------------------
//! set the spec ID. Features with a matching spec ID will be preferred
//! set spec ID to -1 to stop caring
void FECoreKernel::SetSpecID(int nspec)
{
	m_nspec = nspec;
}

//-----------------------------------------------------------------------------
//! set a dependency on a module
bool FECoreKernel::AddModuleDependency(const char* szmodule)
{
	if (m_activeModule == -1) return false;
	FEModule& activeModule = *m_modules[m_activeModule];

	if (szmodule == 0)
	{
		// clear dependencies
		activeModule.ClearDependencies();
		return true;
	}

	// find the module
	for (size_t i = 0; i<m_modules.size(); ++i)
	{
		FEModule& mi = *m_modules[i];
		if (strcmp(mi.GetName(), szmodule) == 0)
		{
			// add the module to the active module's dependency list
			activeModule.AddDependency(mi);

			return true;
		}
	}

	// oh, oh, couldn't find it
	return false;
}

//-----------------------------------------------------------------------------
//! Register a new domain class
void FECoreKernel::RegisterDomain(FEDomainFactory* pf, bool pushFront)
{
	if (pushFront)
		m_Dom.insert(m_Dom.begin(), pf);
	else
		m_Dom.push_back(pf); 
}

//-----------------------------------------------------------------------------
FEDomain* FECoreKernel::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	for (int i=0; i<(int)m_Dom.size(); ++i)
	{
		FEDomain* pdom = m_Dom[i]->CreateDomain(spec, pm, pmat);
		if (pdom != 0) return pdom;
	}
	return 0;
}

//-----------------------------------------------------------------------------
FEDomain* FECoreKernel::CreateDomainExplicit(int superClass, const char* sztype, FEModel* fem)
{
	FEDomain* domain = (FEDomain*)Create(superClass, sztype, fem);
	return domain;
}

//-----------------------------------------------------------------------------
void FECoreKernel::OnCreateEvent(FECreateHandler* pf)
{
	pf->SetModuleID(GetActiveModuleID());
	m_createHandlers.push_back(pf);
}

//-----------------------------------------------------------------------------
void FECoreKernel::BlockEvents(bool b)
{
	m_blockEvents = b;
}
