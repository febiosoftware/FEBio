#pragma once
#include "FEBioFactory.h"
#include <vector>

class FEModel;
class FEBioTask;
class FEBodyForce;
class FEMaterial;

//-----------------------------------------------------------------------------
typedef FEBioFactory_T<FEBioTask>		FEBioTaskFactory;
typedef FEBioFactory_T<FEBodyForce>		FEBodyForceFactory;
typedef FEBioFactory_T<FEMaterial>		FEMaterialFactory;

//-----------------------------------------------------------------------------
// This is the FEBio kernel class which manages the interactions between the 
// different executable modules. In particular, it manages the factory classes
// which are responsible for the life of different classes
class FEBioKernel
{
	struct TASK_DESCR
	{
		const char*			sztag;
		FEBioTaskFactory*	pfac;
	};

	struct BODY_FORCE_DESCR
	{
		const char*			sztag;
		FEBodyForceFactory*	pfac;
	};

	struct MATERIAL_DESCR
	{
		const char* sztag;
		FEMaterialFactory* pfac;
	};

public:
	static FEBioKernel& GetInstance();

public:
	void RegisterTask(FEBioTaskFactory* ptf, const char* sztag);
	FEBioTask* CreateTask(const char* sztag, FEModel* pfem);

	void RegisterBodyForce(FEBodyForceFactory* ptf, const char* sztag);
	FEBodyForce* CreateBodyForce(const char* sztag, FEModel* pfem);

	void RegisterMaterial(FEMaterialFactory* pmf, const char* stag);
	FEMaterial* CreateMaterial(const char* sztag, FEModel* pfem);

protected:
	std::vector<TASK_DESCR>				m_Task;
	std::vector<BODY_FORCE_DESCR>		m_BF;
	std::vector<MATERIAL_DESCR>			m_Mat;

private:
	FEBioKernel(){}
	FEBioKernel(const FEBioKernel&){}
	static FEBioKernel* m_pKernel;
};
