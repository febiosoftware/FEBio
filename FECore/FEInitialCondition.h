#pragma once
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
//! Base class for defining initial conditions.
//! Initial conditions can be used to set the initial state of the model in an analysis. 
class FEInitialCondition : public FEModelComponent
{
public:
	FEInitialCondition(FEModel* pfem);
};

//-----------------------------------------------------------------------------
class FEInitialVelocity : public FEInitialCondition
{
	struct ITEM
	{
		int		nid;	//!< node ID
		vec3d	v0;		//!< initial velocity
	};

public:
	FEInitialVelocity(FEModel* pfem) : FEInitialCondition(pfem){}

	void Serialize(DumpFile& ar);

	void Activate();

	void Add(int nid, vec3d v) { ITEM it = {nid, v}; m_item.push_back(it); }

public:
	vector<ITEM>	m_item;
};

//-----------------------------------------------------------------------------
class FEInitialPressure : public FEInitialCondition
{
	struct ITEM
	{
		int		nid;	//!< node ID
		double	p0;		//!< initial pressure
	};

public:
	FEInitialPressure(FEModel* pfem) : FEInitialCondition(pfem){}

	void Serialize(DumpFile& ar);

	void Activate();

	void Add(int nid, double p) { ITEM it = {nid, p}; m_item.push_back(it); }

public:
	vector<ITEM>	m_item;
};

//-----------------------------------------------------------------------------
class FEInitialConcentration : public FEInitialCondition
{
	struct ITEM
	{
		int		nid;	//!< node ID
		double	c0;		//!< initial concentration
	};

public:
	FEInitialConcentration(FEModel* pfem) : FEInitialCondition(pfem){}

	void Serialize(DumpFile& ar);

	void Activate();

	void Add(int nid, double c) { ITEM it = {nid, c}; m_item.push_back(it); }

	void SetSoluteID(int nsol) { m_nsol = nsol; }

public:
	vector<ITEM>	m_item;
	int				m_nsol;
};

//-----------------------------------------------------------------------------
class FEInitialTemperature : public FEInitialCondition
{
	struct ITEM
	{
		int		nid;	//!< node ID
		double	T0;		//!< initial temperature
	};

public:
	FEInitialTemperature(FEModel* pfem) : FEInitialCondition(pfem){}

	void Serialize(DumpFile& ar);

	void Activate();

	void Add(int nid, double u) { ITEM it = {nid, u}; m_item.push_back(it); }

public:
	vector<ITEM>	m_item;
};

//-----------------------------------------------------------------------------
class FEInitialDilatation : public FEInitialCondition
{
    struct ITEM
    {
        int		nid;	//!< node ID
        double	e0;		//!< initial dilatation
    };
    
public:
    FEInitialDilatation(FEModel* pfem) : FEInitialCondition(pfem){}
    
    void Serialize(DumpFile& ar);
    
    void Activate();
    
    void Add(int nid, double e) { ITEM it = {nid, e}; m_item.push_back(it); }
    
public:
    vector<ITEM>	m_item;
};
