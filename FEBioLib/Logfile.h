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
#include <stdio.h>
#include <string>
#include "LogFileStream.h"

//-----------------------------------------------------------------------------
//! Class that is used for logging purposes

//! This class can output to different 
//! files at the same time.
//! At this time it outputs data to the screen (stdout) and to an external text file.
//! Note that this class is implemented as a singleton, in other words, only one
//! instance can be created.

class FEBIOLIB_API Logfile
{
public:
	enum MODE { LOG_NEVER = 0, LOG_FILE = 1, LOG_SCREEN, LOG_FILE_AND_SCREEN };

public:
	//! constructor
	Logfile();

	//! destructor
	~Logfile();

	//! open a new logfile
	bool open(const char* szfile);

	//! append to existing file
	bool append(const  char* szfile);

	//! formatted printing
	void printf(const char* sz, ...);

	//! print a nice box
	void printbox(const char* sztitle, const char* sz, ...);

	//! set the loggin mode
	MODE SetMode(MODE mode);

	//! get the loggin mode
	MODE GetMode();

	//! flush the logfile
	void flush();

	//! close the logfile
	void close();

	//! return the file name
	const std::string& FileName() { return m_fp->GetFileName(); }

	//! returns if the logfile is ready to be written to
	bool is_valid() { return (m_fp != 0); }

	// set the log stream
	void SetLogStream(LogStream* ps) { m_ps = ps; }

	// set the file log stream
	void SetFileStream(LogFileStream* fp) { m_fp = fp; }

	// return the file handle
	operator FILE* () { return (m_fp ? m_fp->GetFileHandle() : 0); }

private:
	//! copy constructor is private so that you cannot create it directly
	Logfile(const Logfile& log){}

protected:
	LogFileStream*	m_fp;	//!< the actual log file

	LogStream*	m_ps;	//!< This stream is used to output to the screen

	MODE	m_mode;	//!< mode of log file
};
