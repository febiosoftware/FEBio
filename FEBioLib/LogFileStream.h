#pragma once
#include "LogStream.h"
#include "stdio.h"
#include <string>

//-----------------------------------------------------------------------------
// A stream that outputs to a file
class FEBIOLIB_API LogFileStream : public LogStream
{
public:
	// constructor
	LogFileStream();

	// destructor
	~LogFileStream();

	// open the file
	bool open(const char* szfile);

	// open for appending
	bool append(const char* szfile);

	// close the file stream
	void close();

	// get the file handle
	FILE* GetFileHandle() { return m_fp; }

	// get the file name
	const std::string& GetFileName() const { return m_fileName; }

public:
	// print text to the file
	void print(const char* sz);

	// flush the stream
	void flush();

private:
	FILE*		m_fp;
	std::string	m_fileName;
};
