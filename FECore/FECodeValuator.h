/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FEScalarValuator.h>
#include "fecore_api.h"

class FESurface;

class FECORE_API FECodeValuator : public FEScalarValuator
{
	class Imp;

public:
	FECodeValuator(FEModel* fem);
	~FECodeValuator();
	
	bool Init() override;

	int AddGlobalDouble(const std::string& name);

	std::string GetScriptName() const { return m_scriptName; }
	void SetScriptName(const std::string& name) { m_scriptName = name; }

	void CompileDerivative(const std::string& varName);

public:
	double operator()(const FEMaterialPoint& pt) override;
	FEScalarValuator* copy() override;
	bool isConst() override { return false; }
	double* constValue() override { return nullptr; }
	void Serialize(DumpStream& ar) override;

	double run(const FEMaterialPoint& pt, std::vector < std::pair<int, double>>& globals) const;

private:
	bool CompileScript();

private:
	std::string m_scriptName; // user specifies the script name

private:
	Imp& m;
	DECLARE_FECORE_CLASS()
};

// helper function to see if a script compiles. This is used in FEBio studio.
FECORE_API bool ValidateScript(const std::string& script, std::string& err);
