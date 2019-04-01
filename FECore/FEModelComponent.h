/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
//! forward declaration of the FEModel class.
//! All classes inherited from FEModelComponent should take the model as a parameter
//! to the constructor.
class FEModel;

//-----------------------------------------------------------------------------
//! This class serves as a base class for many of the FECore classes. It defines
//! activation and deactivation functions which are used in multi-step analyses to determine which
//! components are active during an analysis.
//! A model component is basically anything that affects the state of a model.
//! For instance, boundary conditions, loads, contact definitions, etc.
//! This class also generates a unique class ID (not be confused with the super class ID)
//! which is used, for instance, by the analysis steps during serialization.
class FECORE_API FEModelComponent : public FECoreBase
{
public:
	//! constructor
	FEModelComponent(FEModel* fem);

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
	bool IsActive() const;

	//-----------------------------------------------------------------------------------
	//! Activate the component.
	//! This function is called during the step initialization, right before the step is solved.
	//! This function can be used to initialize any data that could depend on the model state. 
	//! Data allocation and initialization of data that does not depend on the model state should
	//! be done in Init().
	virtual void Activate();

	//! Deactivate the component
	virtual void Deactivate();

	//-----------------------------------------------------------------------------------
	//! Update the component
	//! This is called whenever the model is updated.
	virtual void Update();

public:

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
	void Serialize(DumpStream& ar);

protected:
	//! used by derived class to set the ID
	void SetID(int n);

private:
	int			m_nID;
	int			m_nClassID;	//!< the class ID
	bool		m_bactive;	//!< flag indicating whether the component is active
};
