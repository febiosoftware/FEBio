import sys
import windowsGoldStandards as win, macOSGoldStandards as mac, linuxGoldStandards as lin

if '-c' in sys.argv:
    BLUE = '\033[96m'
    GREEN = '\033[92m'
    RED = '\033[91m'
    WHITE = '\033[37m'
else:
    BLUE = GREEN = RED = WHITE = ''

smallDiff = set()
bigDiff = set()

def getNumStr(pdiff, key):
    diffString = ""
    
    if abs(pdiff) > 10:
        diffString += RED
        bigDiff.add(key)
    elif abs(pdiff) > 1:
        diffString += BLUE
        smallDiff.add(key)
    elif pdiff != 0:
        diffString += GREEN
        
    # Only use sci notation when necessary 
    # Makes it much easier to read
    if pdiff == 0:
        diffString += "0"
    elif round(pdiff,2) != 0:
        diffString += str(round(pdiff,2)) + WHITE
    else:
        diffString += "{:.2e}".format(pdiff) + WHITE
    
    return diffString
    
def printDiffs(OS, osStd, key):
    length = len(win.stdResults[key])
    
    diffString = "\t" + OS + ":\t["
    for index in range(1, length - 1):
        winVal = win.stdResults[key][index]
        osVal = osStd[key][index]
        pdiff = (osVal-winVal)/winVal*100
        
        diffString += getNumStr(pdiff, key) + ", "
    
    # The last element in the list is a string representation of a float
    # Ensure that neither is an empty string
    if win.stdResults[key][length-1] and mac.stdResults[key][length-1]:
        winVal = float(win.stdResults[key][length-1])
        osVal = float(osStd[key][length-1])
        pdiff = (osVal-winVal)/winVal*100
        
        diffString += getNumStr(pdiff, key) + "]"
    else:
        diffString += "~]"
        
    print(diffString)
    

keys = list(lin.stdResults.keys())
keys.sort()

problemCats = {}

macCount = 0
linCount = 0
totalCount = 0
for key in keys:
    macDiff = win.stdResults[key] != mac.stdResults[key]
    linDiff = win.stdResults[key] != lin.stdResults[key]
    
    if macDiff or linDiff:
        print(key + ":")
        print("\tWin:\t" + str(win.stdResults[key][1:]))
        totalCount += 1
        
        try:
            problemCats[key[0:2]] = problemCats[key[0:2]] + 1
        except:
            problemCats[key[0:2]] = 1
        
    if macDiff:
        macCount += 1
        
        printDiffs("Mac", mac.stdResults, key)

    if linDiff:
        linCount += 1
        
        printDiffs("Lin", lin.stdResults, key)

print()        
print("Problem Categories:")
for key in problemCats.keys():
    print("\t" + key + ": " + str(problemCats[key]))

# Remove problems in bigDiff from smallDiff so they aren't 
# reported twice
smallDiff = smallDiff - bigDiff

print()
smallList = list(smallDiff)
smallList.sort()
print("Problems with 10 > %diff > 1:")
for item in smallList:
        print("\t" + item)
print("Total: " + str(len(smallList)))

print()
bigList = list(bigDiff)
bigList.sort()
print("Problems with %diff > 10:")
for item in bigList:
        print("\t" + item)
print("Total: " + str(len(bigList)))
    
print()
print("MacDiff:\t" + str(macCount))
print("LinDiff:\t" + str(linCount))
print("TotalDiff:\t" + str(totalCount))
