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
