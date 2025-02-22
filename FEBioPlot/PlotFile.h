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



#pragma once
#include "FECore/FEMesh.h"
#include "FECore/FEPlotData.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class implements the facilities to write to a plot database. 
//!
class PlotFile
{
public:
	// size of name variables
	enum { STR_SIZE = 64 };

	// Dictionary entry
	class DICTIONARY_ITEM
	{
	public:
		DICTIONARY_ITEM();
		DICTIONARY_ITEM(const DICTIONARY_ITEM& item);

	public:
		FEPlotData* m_psave;
		unsigned int	m_ntype;	// data type
		unsigned int	m_nfmt;		// storage format
		unsigned int	m_arraySize;	// size of arrays (only used by arrays)
		std::vector<string>	m_arrayNames;	// names of array components (optional)
		char			m_szname[STR_SIZE];
		char			m_szunit[STR_SIZE];
	};

	class Dictionary
	{
	public:
		bool AddVariable(FEModel* pfem, const char* szname, std::vector<int>& item, const char* szdom = "");

		int GlobalVariables() { return (int)m_Glob.size(); }
		int NodalVariables() { return (int)m_Node.size(); }
		int DomainVariables() { return (int)m_Elem.size(); }
		int SurfaceVariables() { return (int)m_Face.size(); }
		int EdgeVariables() { return (int)m_Edge.size(); }

		void Defaults(FEModel& fem);

		void Clear();

	public:
		list<DICTIONARY_ITEM>& GlobalVariableList() { return m_Glob; }
		list<DICTIONARY_ITEM>& MaterialVariableList() { return m_Mat; }
		list<DICTIONARY_ITEM>& NodalVariableList() { return m_Node; }
		list<DICTIONARY_ITEM>& DomainVariableList() { return m_Elem; }
		list<DICTIONARY_ITEM>& SurfaceVariableList() { return m_Face; }
		list<DICTIONARY_ITEM>& EdgeVariableList() { return m_Edge; }

	protected:
		bool AddGlobalVariable(FEPlotData* ps, const char* szname);
		bool AddMaterialVariable(FEPlotData* ps, const char* szname);
		bool AddNodalVariable(FEPlotData* ps, const char* szname, std::vector<int>& item);
		bool AddDomainVariable(FEPlotData* ps, const char* szname, std::vector<int>& item);
		bool AddSurfaceVariable(FEPlotData* ps, const char* szname, std::vector<int>& item);
		bool AddEdgeVariable(FEPlotData* ps, const char* szname, std::vector<int>& item);

	protected:
		list<DICTIONARY_ITEM>	m_Glob;		// Global variables
		list<DICTIONARY_ITEM>	m_Mat;		// Material variables
		list<DICTIONARY_ITEM>	m_Node;		// Node variables
		list<DICTIONARY_ITEM>	m_Elem;		// Domain variables
		list<DICTIONARY_ITEM>	m_Face;		// Surface variables
		list<DICTIONARY_ITEM>	m_Edge;		// Edge variables

		friend class PlotFile;
	};

public:
	//! constructor
	PlotFile(FEModel* fem);

	//! descructor
	virtual ~PlotFile();

	//! close the plot database
	virtual void Close();

	//! Open the plot database
	virtual bool Open(const char* szfile) = 0;

	//! Open for appending
	virtual bool Append(const char* szfile) = 0;

	//! Write current FE state to plot database
	virtual bool Write(float ftime, int flag = 0) = 0;

	//! see if the plot file is valid
	virtual bool IsValid() const = 0;

	virtual void Serialize(DumpStream& ar) {}

public:
	Dictionary& GetDictionary() { return m_dic; }

protected:
	FEModel* GetFEModel() { return m_pfem; }

	// build the dictionary
	void BuildDictionary();

	//! Add a variable to the dictionary
	bool AddVariable(FEPlotData* ps, const char* szname);
	bool AddVariable(const char* sz);
	bool AddVariable(const char* sz, std::vector<int>& item, const char* szdom = "");

private:
	Dictionary	m_dic;	//!< dictionary
	FEModel*	m_pfem;	//!< pointer to FE model
};
