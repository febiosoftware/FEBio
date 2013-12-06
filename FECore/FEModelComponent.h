#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
//! forward declaration of the FEModel class.
//! All classes inherited from FEModelComponent should take the model as a parameter
//! to the constructor.
class FEModel;

//-----------------------------------------------------------------------------
//! This class serves as a base class for many of the FECore classes. It defines
//! activation and deactivation functions which is used in multi-step analyses
class FEModelComponent : public FECoreBase
{
public:
	FEModelComponent(SUPER_CLASS_ID, FEModel* pfem);
	virtual ~FEModelComponent();

	//! return the FE model
	FEModel* GetFEModel();

	//! Is this component active
	bool IsActive();

	//! Activate the component
	virtual void Activate();

	//! Deactivate the component
	virtual void Deactivate();

private:
	bool		m_bactive;	//!< flag indicating whether the component is active
	FEModel*	m_pfem;		//!< model that this component belongs too
};
