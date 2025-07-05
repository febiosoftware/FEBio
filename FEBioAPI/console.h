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
#include <FEBioLib/LogStream.h>
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
//! The Console class manages the shell window. This class is implemented as
//! a singleton, i.e. there can only be one console class in the entire
//! application. Users obtain a pointer to the Console by calling the GetHandle
//! function. 
class Console  
{
public:
	//! return the pointer to the one and only console object
	static Console* GetHandle();

public:
	Console() { m_bActive = true; }
	~Console();

	void CleanUp();

	//! set the title of the console
	void SetTitle(const char* sz, ...);

	void Activate() { m_bActive = true; } 
	void Deactivate() { m_bActive = false; }

	void GetCommand(int& nargs, char** argv);

	//! waits for user input (similar to system("pause"))
	void Wait();

	void Draw(unsigned char* img, int nx, int ny);

	void Write(const char* sz, unsigned short att);

	void SetProgress(double pct);

	const std::vector<std::string>& GetHistory() const;

protected:
	bool	m_bActive;
	std::vector<std::string>	m_history;	//!< command history

protected:
	static	Console*			m_pShell;	//!< pointer to the one and only console class
};

//-----------------------------------------------------------------------------
// we use the console to log output 
class ConsoleStream : public LogStream
{
public:
	void print(const char* sz);
};
