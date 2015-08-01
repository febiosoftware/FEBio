#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
//! forward declaration of the FEModel class.
//! All classes inherited from FEModelComponent should take the model as a parameter
//! to the constructor.
class FEModel;

//-----------------------------------------------------------------------------
//! This class serves as a base class for many of the FECore classes. It defines
//! activation and deactivation functions which is used in multi-step analyses.
//! A model component is basically anything that affects the state of a model.
//! For instance, boundary conditions, loads, contact definitions, etc.
//! This class also generates a unique class ID (not be confused with the super class ID)
//! which is used for instance by the analysis steps during serialization.
class FEModelComponent : public FECoreBase
{
public:
	//! constructor
	FEModelComponent(SUPER_CLASS_ID, FEModel* pfem);

	//! destructor
	virtual ~FEModelComponent();

	//-----------------------------------------------------------------------------------
	//! This function is called during initialization, prior to solving the model.
	//! Classes should use this to allocate data structures. Although this component
	//! may not be used for a while (e.g. if it is not active during the first time step)
	//! classes should still attempt to allocate and initialize all data as well as
	//! perform any error checking. Use the function Activate to initialize any additional
	//! data that depends on the model state.
	virtual bool Init();

	//-----------------------------------------------------------------------------------
	//! This function checks if the component is active in the current step. 
	bool IsActive();

	//-----------------------------------------------------------------------------------
	//! Activate the component.
	//! This function is called during the step initialization, right before the step is solved.
	//! This function can be used to initialize any data that could depend on the model state. 
	//! Data allocation and initialization of data that does not depend on the model state should
	//! be done in Init().
	virtual void Activate();

	//! Deactivate the component
	virtual void Deactivate();

public:
	//! return the FE model
	FEModel* GetFEModel() const;

	//! return the name of the component
	const char* GetName();

	//! Set the name of the component
	void SetName(const char* sz);

	//! Get the class ID
	int GetClassID() const;

	//! set the class ID
	//! Class ID's are set automatically. There is no need to call this function
	//! This function is currently only used during serialization
	void SetClassID(int n);

	//! Get the ID
	int GetID() const;

public:
	//! serialization
	void Serialize(DumpFile& ar);

protected:
	//! used by derived class to set the ID
	void SetID(int n);

private:
	int			m_nID;
	int			m_nClassID;	//!< the class ID
	bool		m_bactive;	//!< flag indicating whether the component is active
	FEModel*	m_pfem;		//!< model that this component belongs too
	char*		m_szname;	//!< name of component
};
