import os

similar = []
dissimilar = []

readmeFile = open("cavities_similar_BO1_clean.tsv", "r")
readmeFile2 = open("cavities_dissimilar_BO1_clean.tsv", "r")

while True:
    line = readmeFile.readline()
    if not line or line.isspace():
        break    
    data = line.split()
    source = data[2] + "_cavity6.mol2"
    target = data[3] + "_cavity6.mol2"
    similar.append(source)
    similar.append(target)

while True:
    line = readmeFile2.readline()
    if not line or line.isspace():
        break    
    data = line.split()
    source = data[2] + "_cavity6.mol2"
    target = data[3] + "_cavity6.mol2"
    dissimilar.append(source)
    dissimilar.append(target)

readmeFile.close()

path = 'output'
if not os.path.exists(path):
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

for i in range(0, len(similar)-1, 2):
    posSource = similar[i].find('.')
    posTarget = similar[i+1].find('.')
    similarSource = "cavities/" + similar[i][0:posSource] + ".mol2"
    similarTarget = "cavities/" + similar[i+1][0:posTarget] + ".mol2"
    sourceSimilarFile = open(similarSource, "r")
    for l in range(6):
        line = sourceSimilarFile.readline()
    pointsNumber = line.split(" ")[2]
    outputFile = "output/similar" + str(int(i/2+1)) + ".txt"
    sourceSimilarFile.close()
    os.system("./GoICP " + similarTarget + " " + similarSource + " " + pointsNumber + " config.txt " + outputFile + " " + str(int(i/2+1)))
    print()
    print("--------------------------------------")
    print()

for i in range(0, 0, 2):
#for i in range(0, len(dissimilar)-1, 2):
    posSource = dissimilar[i].find('.')
    posTarget = dissimilar[i+1].find('.')
    dissimilarSource = "cavities/" + dissimilar[i][0:posSource] + ".mol2"
    dissimilarTarget = "cavities/" + dissimilar[i+1][0:posTarget] + ".mol2"
    sourceDissimilarFile = open(dissimilarSource, "r")
    for l in range(6):
        line = sourceDissimilarFile.readline()
    pointsNumber = line.split(" ")[2]
    outputFile = "output/dissimilar" + str(int(i/2+1)) + ".txt"
    sourceDissimilarFile.close()
    os.system("./GoICP " + dissimilarTarget + " " + dissimilarSource + " " + pointsNumber + " config.txt " + outputFile + " " + str(int(i/2+1)))
    print()
    print("--------------------------------------")
    print()