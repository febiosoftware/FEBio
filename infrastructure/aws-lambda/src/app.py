import sys
import subprocess

def handler(event, context):
    cmd = ["febio4", "-i", "/FEBio/TestSuite/Verify3/bi04.feb", "-o", "/FEBio/TestSuite/Verify3/bi04.log", "-p", "/FEBio/TestSuite/Verify3/bi04.xplt"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
