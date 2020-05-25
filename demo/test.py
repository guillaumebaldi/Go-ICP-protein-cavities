import os
import glob

files = glob.glob("output/*.txt")
#files.sort(key=os.path.getmtime)

for filename in files:
    file = open(filename, "r")
    for i, line in enumerate(file):
        if i == 10:
            pos = line.find(' ')
            pos2 = line.find('\n')
            line = line[pos+1:pos2]
            print(filename, line)
        
    


