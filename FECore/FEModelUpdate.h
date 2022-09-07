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
#include "fecore_api.h"
#include "FECoreBase.h"
#include <functional>

// This interface class is used in the DoModelUpdate create handlers. 
class FECORE_API FEModelUpdate
{
public:
	FEModelUpdate();
	virtual ~FEModelUpdate();

	static FEModelUpdate* Instance();

public:
	virtual void AddPlotVariable(const char* szplt) = 0;

	static FEModelUpdate* m_pThis;
};

// template class the applies a model update when an FECore class is created. 
template <class T> class DoModelUpdate : public FECreateHandler
{
public:
	DoModelUpdate(std::function<void(FEModelUpdate&)> f) : m_f(f) {}

	void handle(FECoreBase* pc) override
	{
		FEModelUpdate* fem = FEModelUpdate::Instance();
		if (fem && dynamic_cast<T*>(pc)) m_f(*fem);
	}

private:
	std::function<void(FEModelUpdate&)>	m_f;
};

// special function that allocates a model update class that adds a plot variable
// when an FECore class of a particular type is allocated.
template <class T> FECreateHandler* AddPlotVariableWhenCreating(const char* szplt)
{
	return new DoModelUpdate<T>([=](FEModelUpdate& fem) { fem.AddPlotVariable(szplt); });
}

template <class T> FECreateHandler* UpdateModelWhenCreating(std::function<void(FEModelUpdate&)> f)
{
	return new DoModelUpdate<T>(f);
}

template <class T> class CreateHandlerFunction : public FECreateHandler
{
public:
	CreateHandlerFunction(std::function<void(T*)> f) : m_f(f) {}

	void handle(FECoreBase* pc) override
	{
		T* pT = dynamic_cast<T*>(pc);
		if (pT) m_f(pT);
	}

private:
	std::function<void(T* pc)> m_f;
};

template <class T> FECreateHandler* CallWhenCreating(std::function<void(T*)> f)
{
	return new CreateHandlerFunction<T>(f);
}
