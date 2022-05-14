/*This file is part of the FEBio Studio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio-Studio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
// units are defined as strings that use the following
// characters to represent the physical quantities that
// define the SI base units.
// (Note that the symbols for time and temperature differ
// from the conventional SI dimension symbols)
//
// L = length
// M = mass
// t = time
// T = temperature
// l = electric current
// n = substance
//
// In addition, the following symbols for derived units are also defined:
//
// F = force
// P = pressure
// E = energy
// W = power
// d = angles in degrees
// r = angles in radians
//
// or use one of these predefined constants
#define UNIT_NONE	""
#define UNIT_LENGTH "L"
#define UNIT_MASS   "M"
#define UNIT_TIME   "t"
#define UNIT_TEMPERATURE "T"
#define UNIT_CURRENT "I"
#define UNIT_SUBSTANCE	"n"
#define UNIT_FORCE "F"
#define UNIT_MOMENT "F.L"
#define UNIT_PRESSURE "P"
#define UNIT_ENERGY	"E"
#define UNIT_POWER "W"
#define UNIT_VOLTAGE "V"
#define UNIT_CONCENTRATION "c"

#define UNIT_DEGREE "d"
#define UNIT_RADIAN	"r"

#define UNIT_AREA   "L^2"
#define UNIT_VOLUME "L^3"
#define UNIT_VELOCITY "L/t"
#define UNIT_ACCELERATION "L/t^2"
#define UNIT_ANGULAR_VELOCITY "r/t"
#define UNIT_DENSITY "M/L^3"
#define UNIT_PERMEABILITY	"L^4/F.t"
#define UNIT_DIFFUSIVITY  "L^2/t"
#define UNIT_MOLAR_MASS "M/n"
#define UNIT_MOLAR_VOLUME "L^3/n"
#define UNIT_GAS_CONSTANT "E/n.T"
#define UNIT_FARADAY_CONSTANT "I.t/n"
#define UNIT_VISCOSITY "P.t"
#define UNIT_FILTRATION "L^2/F.t"
#define UNIT_STIFFNESS "F/L"
#define UNIT_CURRENT_CONDUCTIVITY "I/V.L^2"
