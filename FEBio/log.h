#pragma once

#include "Logfile.h"

//! helper function to obtain a reference to the logfile
inline Logfile& GetLogfile() { return *Logfile::GetInstance(); }
