#!/usr/bin/env python3
# import pdb;

import os, subprocess, sys, glob, datetime, time, shutil, atexit, platform, smtplib, ssl, email, re, fnmatch, copy
from email.message import EmailMessage
from universal import *
from textwrap import TextWrapper

if platform.system() == "Windows":
    from windowsLocals import *
    from windowsGoldStandards import *
elif platform.system() == "Linux":
    from linuxLocals import *
    from linuxGoldStandards import *
else:
    from macOSLocals import *
    from macOSGoldStandards import *

BLUE = '\033[96m'
RED = '\033[91m'
GREEN = '\033[92m'
WHITE = '\033[37m'

def runTests(plugins, numCores=NUMCORES, exp=None):
    
    # Stop subprocesses on death    
    def cleanup():
        for proc in running.values():
            proc.kill()
    atexit.register(cleanup)
    
    # Define the test problems list.
    testFiles = glob.glob(VERIFYDIR + "*.feb")
    testFiles.sort()
    
    # Remove exempt problems from list
    for name in exemptTests:
        try:
            testFiles.remove(VERIFYDIR + name + ".feb")
        except:
            pass
    
    for name in localExemptTests:
        try:
            testFiles.remove(VERIFYDIR + name + ".feb")
        except:
            pass
    
    # Manually move some long problems to the top of the list to help 
    # the suite run faster
    for name in longTests:
        try:
            testFiles.remove(VERIFYDIR + name + ".feb")
            testFiles.insert(0, VERIFYDIR + name + ".feb")
        except:
            pass

    err = []
    ioErr = []
    running = {}
    results = {}
    current = 0
    # Loop over the tests and run them
    
    tic = time.time()
    while True:
        while current < len(testFiles) and len(running) < numCores:
            fileName = testFiles[current]
            if OSNAME=="WINDOWS":
                baseName = fileName.split('\\')[-1].split('.')[0]
            else:
                baseName = fileName.split('/')[-1].split('.')[0]

            current += 1
            
            # skip optimization files
            if baseName[:2] == "oi":
                continue
            
            # Skip tests not matching optional exp argument
            if exp:
                cont = True
                for ex in exp:
                    if re.search(ex, baseName):
                        cont = False
                        break
                if cont:
                    continue
            
            # Define the log and plt files
            logname = WORKINGDIR + baseName + '.log'
            pltname = WORKINGDIR + baseName + '.xplt'
        
            opt = False
        
            # Test for parameter optimization problems
            if baseName[:2] == "op":
                opt = True
                optFile = VERIFYDIR + "oi" + baseName[2:] + ".feb"
                command = [FEBIOUPLOADPATH, '-i', optFile, '-s', fileName, '-o', logname, '-p', pltname]
            # Test for plugin problems
            elif baseName in pluginTests.keys():
                if not plugins:
                    continue
                pluginPath = pluginPaths[pluginTests[baseName]]
                command = [FEBIOUPLOADPATH, '-import', pluginPath, '-i', fileName, '-o', logname, '-p', pltname]
            else:
                command = [FEBIOUPLOADPATH, '-i', fileName, '-o', logname, '-p', pltname]
                
            # Create a variable that will store the results of the test
                
            #Optimization input file
            # 0: Normal/Error termination status
            # 1: Final objective value
            # 2: Total iterations
            # 3: Plot file size
            if opt: results[baseName] = ["", 0.0, 0, 0]

            # Normal input file
            # 0: Normal/Error termination status
            # 1: Number of time steps
            # 2: Number of iterations
            # 3: Number of RHS evaluations
            # 4: Number of reformations
            # 5: Plot file size
            # 6: Data field value
            else: results[baseName] = ["", 0, 0, 0, 0, 0, ""]
            
            # Run the test problem
            print(' '.join(command))
            
            cmd_env = os.environ.copy()
            cmd_env["OMP_NUM_THREADS"] = "1"
            running[baseName] = subprocess.Popen(command, env=cmd_env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        if len(running) == 0:
            break
                
        finished = []
        for baseName in running:
            if running[baseName].poll() != None:
                finished.append(baseName)

        for baseName in finished:
            # Define the log and plt files
            logname = WORKINGDIR + baseName + '.log'
            pltname = WORKINGDIR + baseName + '.xplt'
            
            # Check if test is an optimization problem
            opt = "op" in baseName
            
            # If this file is in the dataField list
            try:
                df_time = dataField[baseName]
                df_tline = "Time = " + df_time + "\n"
                df_flg = 1
                found = 0
                data1 = 0
                line_num = 0
            # Else
            except:
                df_flg = 0
            
            # check the return value
            if running[baseName].returncode == 0:
                results[baseName][0] = 'Normal'
            else:
                err.append(baseName)
                results[baseName][0] = 'Error'
                
            # search the log file for the convergence info
            try:
                flog = open(logname, 'r')
                
                if opt:
                    for line in flog:
                        if line.find("Final objective value") !=-1: results[baseName][1] = line.split()[5]
                        if line.find("Total iterations"     ) !=-1: results[baseName][2] = int(line[43:])
                else:
                    for line in flog:
                        if  line.find("Number of time steps completed"        ) != -1: results[baseName][1] = int(line[55:])
                        if  line.find("Total number of equilibrium iterations") != -1: results[baseName][2] = int(line[55:])
                        if  line.find("Total number of right hand evaluations") != -1: results[baseName][3] = int(line[55:])
                        if  line.find("Total number of stiffness reformations") != -1: results[baseName][4] = int(line[55:])
                        if  line.find("Time in linear solver") != -1:
                            slv_hr  = int(line[24:25])
                            slv_min = int(line[26:28])
                            slv_sec = int(line[29:31])
                            new_slv_time = slv_hr*3600 + slv_min*60 + slv_sec
                            #print "New solve time", new_slv_time
                        if  line.find("elapsed time") != -1:
                            el_hr  = int(line[37:38])
                            el_min = int(line[39:41])
                            el_sec = int(line[42:44])
                            new_el_time = el_hr*3600 + el_min*60 + el_sec
                        if df_flg:
                            if line.find("Data Record #1") !=-1: data1 = 1
                            if data1: line_num += 1
                            if line.find(df_tline) !=-1 and data1: found = 1
                            if line_num == 6:
                                if found: results[baseName][6] = line.rstrip("\n").split(" ")[1]
                                found = 0
                                line_num = 0
                                data1 = 0

                # get the size of the plotfile and delete it
                results[baseName][5-2*opt] = int(os.path.getsize(pltname))
            except IOError:
                pass
                
                
            # Print convergence info
            # determine which, if any, convergence criteria failed
            fail = [False] * len(results[baseName])
            
            try:
                index = 0
                for index in range(1, len(results[baseName])):
                    if results[baseName][index] != stdResults[baseName][index]:
                        fail[index] = True
            except:
                pass
                  
            # if the test error terminated
            if results[baseName][0] != "Normal":
                print(RED + baseName + ": " + str(results[baseName]) + WHITE)
            # if the test failed convergence criteria
            elif True in fail:                
                resultStrings = str(results[baseName]).replace("[", "").replace("]", "").split(",")
                
                for index in range(1, len(fail)):
                    if fail[index]:
                        resultStrings[index] += "(" + str(stdResults[baseName][index]) + ")"
                        
                msg = '['
                for rStr in resultStrings:
                    msg += rStr + ","
                    
                msg = msg[:-1] +']'

                print(BLUE + baseName + ": " + msg + WHITE)
            # if the test passed
            else:
                print(GREEN + baseName + ": " + str(results[baseName]) + WHITE)
            
            del running[baseName]
            
        time.sleep(0.01)
        
    toc = time.time()
    
    # Calculate time it took to run suite
    hours, remainder = divmod(int(toc-tic), 3600)
    minutes, seconds = divmod(remainder, 60)
    
    elapsedTime = '%s:%s:%s' % (hours, minutes, seconds)
    
    print("Test Suite ran in " + elapsedTime)

    # Find which tests failed the convergence criteria, and which cirteria they failed
    failed = {}
    for key in results:
        try:
            fail = [False] * len(results[key])
            
            if results[key][0] != "Normal":
                # fail[0] = True
                # failed[key] = fail
                continue
            
            index = 0
            for index in range(1, len(results[key])):
                if results[key][index] != stdResults[key][index]:
                    fail[index] = True
                    
            if True in fail:
                failed[key] = fail
        except:
            pass
    
    success = True
    subject = ""
    message = ""
    
    if len(err) > 0:
        success = False
        
        subject += "Error terminations or crashes"
        
        message += "The following tests error terminated, or crashed:\n"
        for name in err:
            message += "\t" + name + "\n"
        message += "\n"
    
    if len(failed) > 0:
        success = False
        
        if subject != "":
            subject += ", "
        
        subject += "Bad Convergence Criteria"
        
        message += "The following tests failed due incorrect convergence criteria:\n\n"
        
        for name in failed:
            if name not in stdResults.keys():
                continue
            
            failedMsg = "    " + name + ": "
            
            resultStrings = str(results[name]).replace("[", "").replace("]", "").split(",")
            
            for index in range(1, len(failed[name])):
                if failed[name][index]:
                    resultStrings[index] += "(" + str(stdResults[name][index]) + ")"
            
            if "op" in name:
                failedMsg += "{:^15}|{:^15} | ".format(*resultStrings[1:-1])
                failedMsg += resultStrings[-1]
            else:
                failedMsg += "{:^15}|{:^15}|{:^15}|{:^15}|{:^20} | ".format(*resultStrings[1:-1])
                failedMsg += resultStrings[-1]
                
            
            message += failedMsg + "\n"
        
        message += "\n"
        
    if success:
        subject = "All tests passed"
        message = "All tests in the suite passed.\n"
        
    message += "Test Suite ran in " + elapsedTime + "\n\n"
        
    # Remove crashed tests from results so that they don't get added as new tests
    for name in err:
        del results[name]
        
    # Add any new tests to the std file
    if len(results) > len(stdResults):
        
        if subject != "":
            subject += ", "
        
        subject += "New Tests Added"
        
        message += "The following tests were newly added to the test suite:"
        
        for key in results:
            if key not in stdResults.keys() and results[key][0] == "Normal":
                message += "\t" + key + "\n"
                
                stdResults[key] = results[key]
    
        with open(AUTOMATIONDIR + GOLDSTANDARDS, "w") as ldata:
            ldata.write("stdResults = " + str(stdResults).replace("], ", "],\n        ") + "\n\n")
            
        # Commit changes and push them
        subprocess.run(["git", "-C", AUTOMATIONDIR, "stage", AUTOMATIONDIR + GOLDSTANDARDS], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run(["git", "-C", AUTOMATIONDIR, "commit", "-m", "Updated " + GOLDSTANDARDS + " with new tests."], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run(["git", "-C", AUTOMATIONDIR, "pull", "--no-edit"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run(["git", "-C", AUTOMATIONDIR, "push"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
    
    with open(LOGDIR + str(datetime.date.today()) + ".txt", "w") as log:
        log.write(subject)
        log.write(message)
        log.write(str(results).replace("], ", "],\n"))
        
    return success, subject, message


def build(name, baseDir, buildDir, path):
    print("Building " + name)
    os.chdir(buildDir)
    
    try:
        oldTime = os.path.getmtime(path)
    except:
        oldTime = 0

    #cmakeout = subprocess.run(["cmake", baseDir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #if cmakeout.returncode != 0:
    #    return False, cmakeout
    #
    ## The build and clean commands for windows have to include the name of the solution. 
    #if platform.system() == "Windows":
    #    buildCommand = BUILDCOMMAND.copy()
    #    buildCommand.insert(1, name + ".sln")
    #    
    #    cleanCommand = CLEANCOMMAND.copy()
    #    cleanCommand.insert(1, name + ".sln")
    #else:
    #    buildCommand = BUILDCOMMAND
    #    cleanCommand = CLEANCOMMAND
    #
    #buildout = subprocess.run(buildCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    ## If the first build fails, clean it then try again
    ## I don't clean every time because it makes it annoying to test
    #if buildout.returncode != 0:
    #    buildout = subprocess.run(cleanCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #    buildout = subprocess.run(buildCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #if buildout.returncode == 0:
    #    print("Build successful")
    #else:
    #    print("Build Failed")
    #    return False, buildout
    #
    newTime = os.path.getmtime(path)

    return True, newTime > oldTime
    
def macPostBuild():
    print("Running postbuild script")
    os.system(FBSBUILDDIR + "postbuild.sh >/dev/null 2>&1")
    
def update(name, path, branch = None):
    print("Updating " + name)
    # os.chdir(path)

    # if branch:
    #     gitout = subprocess.run(["git", "checkout", branch], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #     if gitout.returncode == 0:
    #         print("Switched to " + branch + " branch.")
    #     else:
    #         return False, gitout

    # gitout = subprocess.run(["git", "pull"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # if gitout.returncode == 0:
    #     print("Update successful")
    # else:
    #     print("Update Failed")
    #     return False, gitout
    #     
    return (True,)

def upload(name, path, remoteDir):
    success = False
    attempts = 0
    while attempts < 5:
        print("Uploading " + name + " to server. Atempt: " + str(attempts + 1))
        uploadOut = subprocess.run(["scp", path, "root@repo.febio.org:" + remoteDir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if uploadOut.returncode == 0:
            success = True
            break
            
        attempts += 1
    
    if not success:
        print("Failed to upload " + name)
    
    return success
    
def makeRelease():
    success = False
    attempts = 0
    while attempts < 5:
        print("Running makeDevRelease script on server. Attempt: " + str(attempts + 1))
        sshOut = subprocess.run(["ssh", "root@repo.febio.org", REMOTESCRIPT], 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
        if sshOut.returncode == 0:
            success = True
            break
        else:
            print("Failed to upload FEBio Studio")
            
        attempts += 1
        
    return success
    
def sendEmail(success, subject, message, gitInfo = False):
    msg = EmailMessage()    

    if success:
        msg['Subject'] = OSNAME + " PASSED - " + subject
    else:
        msg['Subject'] = OSNAME + " FAILED - " + subject
    
    if gitInfo:
        gitout = subprocess.run(["git", "-C", FEBIODIR, "log", '--after="1 day ago"'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        gitMessage = "The following commits to the FEBio repository occured in the last 24 hours:\n\n"
        gitMessage += "    " + str(gitout.stdout, "Utf-8").replace("\n", "\n    ") + "\n"
        
        message = gitMessage + message
    
    msg.set_content(message)
    
    print("Subject: " + subject)
    print("Message: " + message)
    
    port = 465  # For SSL
    sender_email = "febioreports@gmail.com"
    password = "reportsadmin2021{}"
    receiver_email = "febio-test@sci.utah.edu"
    
    msg['From'] = sender_email
    msg['To'] = receiver_email

    # Create a secure SSL context
    ssl._create_default_https_context = ssl._create_unverified_context
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
        server.login(sender_email, password)

        server.send_message(msg)
        
def acceptChanges(logFile, exp = None):
    if not os.path.isfile(logFile):
        print(logFile + " does not exist.")
        
    oldResults = copy.deepcopy(stdResults)
    
    updatedTests = set()

    with open(logFile) as log:
        for line in log:
            if not line.startswith('    '):
                continue
            
            splitLine = line.split(':');
            
            test = splitLine[0].strip();         
            results = splitLine[1].split('|')
            
            # Skip tests not matching optional exp argument
            if exp:
                cont = True
                for ex in exp:
                    if re.search(ex, test):
                        cont = False
                        break
                if cont:
                    continue
                    
            for index in range(len(results)):
                # If an incorrect result is recorded, the correct result will
                # be shown in parenthesis next to it
                if '(' in results[index]:
                    # If we find an incorrect result, add the test to the
                    # list of updated tests
                    updatedTests.add(test)
                    
                    result = results[index].split('(')[0].strip()
                    
                    if "'" in result:
                        stdResults[test][index + 1] = result.replace("'", '')
                    else:
                        stdResults[test][index + 1] = int(result)
                    
    with open(AUTOMATIONDIR + GOLDSTANDARDS, "w") as ldata:
        ldata.write("stdResults = " + str(stdResults).replace("], ", "],\n        ") + "\n\n")
        
    if len(updatedTests) == 0:
        print("No results were updated.")
    else:
        print("The following results were updated:")
        
        sortedList = list(updatedTests)
        sortedList.sort()
        for test in sortedList:
            print(test + ":")
            print("\tOld: " + str(oldResults[test]))
            print("\tNew: " + str(stdResults[test]))
            
        
        # Commit changes and push them
        print("Committing new standard results.")
        subprocess.run(["git", "-C", AUTOMATIONDIR, "stage", AUTOMATIONDIR + GOLDSTANDARDS], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run(["git", "-C", AUTOMATIONDIR, "commit", "-m", "Updated " + GOLDSTANDARDS + " with new results."], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Pulling and auto-merging")
        subprocess.run(["git", "-C", AUTOMATIONDIR, "pull", "--no-edit"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Pushing new commit")
        subprocess.run(["git", "-C", AUTOMATIONDIR, "push"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def showHelp():
    wrapper = TextWrapper()
    wrapper.width = os.get_terminal_size()[0]
    wrapper.initial_indent = wrapper.subsequent_indent = "              "
    
    print("tools.py accepts the following flags:\n")
    
    print("-h            Show this help text.")
    
    print("")
    print("-a (logfile)  Accept Changes in logfile")
    print(wrapper.fill("Accepts all changes to convergence criteria as reported in (logfile). "
        "If (logfile) is not provided, defaults to newest file in LOGDIR. "
        "Updates the gold standard results for the current OS and then pushes those changes "
        "to the repository."))

    print("")
    print("-r            Run test suite.")
    print(wrapper.fill("Runs the problems in the test stute. Does not update or build anything. "
        "Does not send email. It will update the gold standards and push them to the repository "
        "if new tests are found."))
    print("")
    print(wrapper.fill("Also accepts optional argument to specify the number of simultaneous tests to run:"))
    print("")
    print("    -c [num]  Runs suite with [num] simultaneous problems.")
    
    print("")
    print("-e [exp] ...  Selects problems matching [exp]")
    print(wrapper.fill("Limits the other functions provided by tools.py to problems that match any "
        "number of regular expressions. For example:"))
    print(wrapper.fill("To limit runs to only biphasic and contact problems: "))
    print(wrapper.fill("    tools.py -r -e bi co"))
    print(wrapper.fill("To only accept changes from a few specific probelms: "))
    print(wrapper.fill("    tools.py -a [logfile] -e dm17 co01 fl36"))

def searchFiles(searchString):
    print("Searching files in ", VERIFYDIR)
    files = os.listdir(VERIFYDIR)
    matches = 0
    for filename in files:
        if fnmatch.fnmatch(filename, "*.feb"):
            file = open(VERIFYDIR + "/" + filename)
            linenr = 0
            for line in file:
                linenr += 1
                if re.search(searchString, line):
                    print(filename + " (line " + str(linenr) + "):" + line[:-1])
                    matches += 1
    if matches==0:
        print("no matches found")


if __name__ == "__main__":
    
    exp = []

    if '-e' in sys.argv:
        index = sys.argv.index('-e')
        
        for i in range(index+1, len(sys.argv)):
            if sys.argv[i][0] == '-' and len(sys.argv[i]) == 2:
                break
                
            exp.append(sys.argv[i])
        
        if not exp:
            print("Pass in regular expression after '-e'")
    
    
    if '-a' in sys.argv:
        index = sys.argv.index('-a')
        
        if len(sys.argv) > index + 1:
            acceptChanges(sys.argv[index+1], exp)
        else:
            logs = glob.glob(LOGDIR + "*.txt")
            latestFile = max(logs, key=os.path.getctime)
            
            print("Accepting changes from latest log file: " + latestFile + "\n")
            acceptChanges(latestFile, exp)
    
    elif '-r' in sys.argv:
        plugins = True
        numCores = NUMCORES
        
        if '-c' in sys.argv:
            index = sys.argv.index('-c')
            
            coreFail = False
            
            if len(sys.argv) > index:
                try:
                    numCores = int(sys.argv[index+1])
                except:
                    coreFail = True
                    pass
            else:
                coreFail = True
                
            if numCores <= 0:
                coreFail = True
            
            if coreFail:
                numCores = NUMCORES
                print("Cannot parse cores number. Defauling to " + str(NUMCORES) + " cores.\nPlease pass the number of cores after the '-c' flag as a positive integer.")
                
        runTests(plugins, numCores, exp)
    elif '-h' in sys.argv:
        showHelp()
    elif '-s' in sys.argv:
        index = sys.argv.index('-s')
        searchFiles(sys.argv[index + 1])
    else:
        showHelp()
