import sys

for fname in sys.argv[1:]:
    print(fname)
    lines = []
    with open(fname, "r") as f:
        lines = f.readlines()
        
    for line in range(len(lines)):
        if lines[line].startswith("Copyright (c) 20"):
            lines[line] = "Copyright (c) 2021 University of Utah, The Trustees of Columbia University in\n"
            break
            
    with open(fname, "w") as f:
        f.writelines(lines)
