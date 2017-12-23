#pragma once
#include "FEBiphasic.h"
#include <map>

typedef std::map<int,double> idmap;     //!< map integer id with double value
typedef std::map<int,double>::iterator itridmap;

//-----------------------------------------------------------------------------
// This class implements a material that has a solvent supply following
// Starling's equation

class FESolventSupplyStarling :	public FESolventSupply
{
public:
	//! constructor
	FESolventSupplyStarling(FEModel* pfem);
	
	//! Solute supply
	double Supply(FEMaterialPoint& pt) override;
	
	//! Tangent of supply with respect to strain
	mat3ds Tangent_Supply_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of supply with respect to pressure
	double Tangent_Supply_Pressure(FEMaterialPoint& mp) override;
	
	//! Tangent of supply with respect to concentration
	double Tangent_Supply_Concentration(FEMaterialPoint& mp, const int isol);
	
    //! set parameter attribute for indexed solute parameters
	bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval) override;
    
	//! set value of indexed parameters
	void SetIndexedParameter(idmap& p, int id, double val) { p.insert(std::pair<int, double>(id, val)); }
    
public:
	double		m_kp;				//!< coefficient of pressure drop
	double		m_pv;				//!< prescribed (e.g., vascular) pressure
	vector<double>		m_qc;       //!< coefficients of concentration drops
	vector<double>		m_cv;       //!< prescribed (e.g., vascular) concentrations
    double  m_qctmp;                //!< helper variable for reading in m_qc
    idmap	m_qcinp;                //!< m_qc for each solute (input)
    double  m_cvtmp;                //!< helper variable for reading in m_cv
    idmap	m_cvinp;                //!< m_cv for each solute (input)
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
