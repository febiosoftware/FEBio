import os, shutil, json

os.mkdir("plugins")

for dirName in os.listdir("pluginRepos"):
    dirPath = os.path.join("pluginRepos", dirName)
    
    if os.path.isdir(dirPath):
        newDir = os.path.join("plugins", dirName)
        os.mkdir(newDir)

        for root, dirs, files in os.walk(dirPath, followlinks=True):
            for name in files:
                if name.endswith(".dll") or name.endswith(".dylib") or name.endswith(".so"):
                    filename = name.split("/")[-1].split("\\")[-1]

                    shutil.copy2(os.path.join(root,name), os.path.join(newDir, filename))

def getVersion(path):
    if not path:
        return None

    major = minor = patch = None
    with open(path, "r") as f:
        for line in f:
            if "#define VERSION" in line:
                major = line.split()[-1].strip()
            elif "#define SUBVERSION" in line:
                minor = line.split()[-1].strip()
            elif "#define SUBSUBVERSION" in line:
                patch = line.split()[-1].strip()

    if major and minor and patch:
        return f"{major}.{minor}.{patch}"
    
    return None

versionInfo = {}

# FEBio version
versionInfo["febio"] = getVersion("febio4-sdk/include/FEBioLib/version.h")

# Find version for each plugin
for dirName in os.listdir("pluginRepos"):
    dirPath = os.path.join("pluginRepos", dirName)
    
    if os.path.isdir(dirPath):
        versionHeader = None
        for root, dirs, files in os.walk(dirPath):
            for name in files:
                if name.endswith("version.h"):
                    versionHeader = os.path.join(root, name)
                    break
                
            if versionHeader:
                break
        
        if versionHeader:
            versionInfo[dirName] = getVersion(versionHeader)
        else:
            print(f"Unable to find version header for {dirPath}")

with open("plugins/versions.json", "w") as f:
    json.dump(versionInfo, f)