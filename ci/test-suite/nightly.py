#!/usr/bin/env python3
from tools import *

updateResults = {}
buildResults = {}

updateResults["FEBio"] = update("FEBio", FEBIODIR, "develop")
buildResults["FEBio"] = build("FEBio", FEBIODIR, FEBIOBUILDDIR, FEBIOPATH)
    
# updateResults["FEBioChem"] = update("FEBioChem", CHEMDIR, "master")
# buildResults["FEBioChem"] = build("FEBioChem", CHEMDIR, CHEMBUILDDIR, CHEMPATH)
# 
# updateResults["FEBioHeat"] = update("FEBioHeat", HEATDIR, "master")
# buildResults["FEBioHeat"] = build("FEBioHeat", HEATDIR, HEATBUILDDIR, HEATPATH)

# updateResults["FEBioStudio"] = update("FEBioStudio", FBSDIR, "develop")
# buildResults["FEBioStudio"] = build("FEBioStudio", FBSDIR, FBSBUILDDIR, FBSPATH)

# if OSNAME == "MACOS":
#     macPostBuild()

updateResults["TestSuite"] = update("TestSuite", TESTDIR)

# Only run the test suite if FEBio updated and built properly
runToday = updateResults["FEBio"][0] and buildResults["FEBio"][0]

if runToday:
    # Only run the plugin problems if they updated and built properly
    # runPlugins = updateResults["FEBioChem"][0] and buildResults["FEBioChem"][0] and updateResults["FEBioHeat"][0] and buildResults["FEBioHeat"][0]

    runPlugins = False
    print(runPlugins)
    testResults = runTests(runPlugins)
    exit()
else:
    testResults = (False,)

newFEBio = updateResults["FEBio"][0] and buildResults["FEBio"][0] and buildResults["FEBio"][1]
# newFBS = updateResults["FEBioStudio"][0] and buildResults["FEBioStudio"][0] and buildResults["FEBioStudio"][1]
# newChem = updateResults["FEBioChem"][0] and buildResults["FEBioChem"][0] and buildResults["FEBioChem"][1]
# newHeat = updateResults["FEBioHeat"][0] and buildResults["FEBioHeat"][0] and buildResults["FEBioHeat"][1]

mkDevRelease = False

uploadResults = {}

# if newFBS:
#     uploadResults["FEBioStudio"] = upload("FEBioStudio", FBSPATH, REMOTE_DEV_DIR + FBSREMOTEDIR)
#     
#     if uploadResults["FEBioStudio"]:
#         mkDevRelease = True
# 
# # Only upload FEBio and plugins if the tests passed. 
# if testResults[0]:
#     if newFEBio:
#         uploadResults["FEBio"] = upload("FEBio", FEBIOUPLOADPATH, REMOTE_DEV_DIR + FEBIOREMOTEDIR)
#         
#         if uploadResults["FEBio"]:
#             mkDevRelease = True
#         
#     if newChem:
#         uploadResults["FEBioChem"] = upload("FEBioChem", CHEMPATH, REMOTE_DEV_DIR + CHEMREMOTEDIR)
#         
#         if uploadResults["FEBioChem"]:
#             mkDevRelease = True
#         
#     if newHeat:
#         uploadResults["FEBioHeat"] = upload("FEBioHeat", HEATPATH, REMOTE_DEV_DIR + HEATREMOTEDIR)
#         
#         if uploadResults["FEBioHeat"]:
#             mkDevRelease = True
#     
# # Make the dev release
# if mkDevRelease:
#     devSucess = makeRelease()
# 
# # Construct Email
# success = True
# subject = ""
# message = ""
# 
# # Used to make sure that messages in the subject are concatenated cleanly
# def appendSubject(sub):
#     global subject
#     if not subject:
#         subject = sub
#     else:
#         subject += ", " + sub
# 
# # See if anything failed to update
# updatesFailed = False
# for result in updateResults.values():
#     if not result[0]:
#         success = False
#         updatesFailed = True
#         break
#         
# if updatesFailed:
#     appendSubject("Updates Failed")
#     message += "The following repositories failed to update:\n\n"
#     
#     for name, result in updateResults.items():
#         if not result[0]:
#             message += name + ":\n"
#             message += "\n\nSTDOUT:\n\n" + str(result[1].stdout, "Utf-8") + "\n\nSTDERR:\n\n" + str(result[1].stderr, "Utf-8")
#             message += "\n"
#             
#     message += "\n"
# 
# # See if anything failed to build
# buildsFailed = False
# for result in buildResults.values():
#     if not result[0]:
#         success = False
#         buildsFailed = True
#         break
#         
# if buildsFailed:
#     appendSubject("Builds Failed")
#     message += "The following failed to build:\n\n"
#     
#     for name, result in buildResults.items():
#         if not result[0]:
#             message += name + ":\n"
#             message += "\n\nSTDOUT:\n\n" + str(result[1].stdout, "Utf-8") + "\n\nSTDERR:\n\n" + str(result[1].stderr, "Utf-8")
#             message += "\n"
#             
#     message += "\n"
# 
# # Add test resutls
# if runToday:
#     if not testResults[0]:
#         success = False
#         
#     appendSubject(testResults[1])
#     message += testResults[2]
# else:
#     appendSubject("Tests not run")
#     message += "The test suite was not run.\n\n"
# 
# 
# # See if anything failed to upload
# uploadsFailed = False
# for result in uploadResults.values():
#     if not result:
#         success = False
#         uploadsFailed = True
#         break
#         
# if uploadsFailed:
#     appendSubject("Uploads Failed")
#     message += "The following failed to Upload:\n\n"
#     
#     for name, result in buildResults.items():
#         if not result:
#             message += name + ":\n"
#             
#     message += "\n"
#     
# if mkDevRelease:
#     if devSucess:
#         appendSubject("Dev Release Made")
#         message += "A development release was made with the following files:\n\t"
#         
#         if newFEBio and testResults[0] and uploadResults["FEBio"]:
#             message += "FEBio "
#         
#         if newFBS and uploadResults["FEBioStudio"]:
#             message += "FEBioStudio "
#         
#         if newChem and testResults[0] and uploadResults["FEBioChem"]:
#             message += "FEBioChem "
# 
#         if newHeat and testResults[0] and uploadResults["FEBioHeat"]:
#             message += "FEBioHeat"
#             
#     else:
#         success = False
#         appendSubject("Dev Release Failed")
#         message += "Failed to make dev release. makeDevRelease script failed to run on server."
#         
# sendEmail(success, subject, message, True)
