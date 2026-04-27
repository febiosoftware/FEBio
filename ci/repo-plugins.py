import os, json

HOMEPATH = "/root"

try:
    with open("plugins/versions.json", "r") as file:
        versionInfo = json.load(file)

    febioVersion = versionInfo["febio"]
    del versionInfo["febio"]

    platform = os.environ["OS"]
    if platform == "Windows":
        osFlag = "-w"
    elif platform == "macOS":
        osFlag = "-m"
    else:
        osFlag = "-l"

    for name in versionInfo:
        os.system(f"scp plugins/{name}/* repo:{HOMEPATH}/pluginRepo/files/{name}/develop/stage/")
        os.system(f'ssh repo "python3 {HOMEPATH}/modelServer/plugins.py -d {name} {versionInfo[name]} {febioVersion} {osFlag}"')

except FileNotFoundError:
    print("Error: 'plugins/versions.json not found. ")
except json.JSONDecodeError:
    print("Error: Invalid JSON format in 'plugins/versions.json'.")