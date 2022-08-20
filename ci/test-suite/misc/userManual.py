import os, glob, shutil, re, platform

if platform.system() == "Windows":
    from windowsLocals import *
elif platform.system() == "Linux":
    from linuxLocals import *
else:
    from macOSLocals import *
    
FEBioDocFolder = FEBIODIR + "Documentation/"
lxyFile = FEBioDocFolder + "FEBio_User_Manual.lyx"

ManualsFolder = "/home/mherron/Projects/Manuals/"
FEBioUserFolder = ManualsFolder + "FEBioUser/"

FBSDocFolder = FBSDIR + "Documentation/"
WebDefinesPath = FBSDIR + "FEBioStudio/WebDefines.h"

MANUAL_PREFIX = "FEBio_um_3-5"

# Create a dictionary relating part numbers to part titles by traversing the 
# table of contents at the top of the page.
numNames = {}

order = []

with open(FEBioUserFolder + MANUAL_PREFIX + "-toc.html") as toc:
    tocLines = toc.readlines()
    
    for line in tocLines:
        if 'href="' + MANUAL_PREFIX + '-' in line:
            if "Bibliography.html" in line:
                continue
            
            number = line.split('href="' + MANUAL_PREFIX + '-')[1].split(".html#toc")[0]
            if "-" in number:
                number = number.split('-')[1]
            
            
            name = line.split('target="contents">')[1][:-5]
            if ":" in name:
                name = name[name.find(':') + 2:]
                
            if "</span>" in name:
                name = name.split("</span>")[1]
                
            name = name.replace(" ", "_").replace("-", "_").replace(",", "").replace(":", "").replace("’", "").replace("‘", "").replace("&#8216;", "").replace("&#8217;", "")
            
            numNames[number] = name
            order.append(number)
            

# ~ for item in partNames:
    # ~ print(item + ": " + partNames[item])
    
def num2Link(number):
    sectionType = ""
    
    dots = number.count(".")
    
    if dots == 0:
        try:
            int(number)
            sectionType = "-Chapter"
        except:
             sectionType =  "-Appendix"
    elif dots == 1:
        sectionType = "-Section"
    elif dots == 2:
        sectionType = "-Subsection"
        
    return MANUAL_PREFIX + sectionType + "-" + number + ".html" 
        

linkNames = {}
order2 = []
for number in order:
    name = numNames[number].upper() + "_HTML"
    
    tempNum = number
    while name in order2:
        if "." in tempNum:
            tempNum = tempNum[:tempNum.rfind(".")]
            name = numNames[tempNum].upper() + "_" + name
        else:
            print("We have a problem")
            sys.exit()
        
    link = num2Link(number)
    
    linkNames[name] = link
    order2.append(name)

with open(WebDefinesPath, "w") as f:
    f.write('''/*This file is part of the FEBio Studio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio-Studio.txt for details.

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

#define MANUAL_PATH "https://help.febio.org/FebioUser/"

''')
    
    for name in order2:
        f.write("#define " + name + ' "' + linkNames[name] + '"\n')
