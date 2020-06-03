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

#for i in range(0, 2, 2):
for i in range(0, len(similar)-1, 2):
    #print(similar[i] + "-" + similar[i+1]) 
    posSource = similar[i].find('.')
    posTarget = similar[i+1].find('.')
    similarSource = "cavitiesN/" + similar[i][0:posSource] + "_sim" + str(int(i/2+1)) + "N.xyz"
    similarTarget = "cavitiesN/" + similar[i+1][0:posTarget] + "_sim" + str(int(i/2+1)) + "N.xyz"
    sourceSimilarFile = open(similarSource, "r")
    pointsNumber = sourceSimilarFile.readline().replace("\n", "")
    outputFile = "output/similar" + str(int(i/2+1)) + ".txt"
    sourceSimilarFile.close()
    os.system(" ../GoICP " + similarTarget + " " + similarSource + " " + pointsNumber + " config.txt " + outputFile)

for i in range(0, 0, 2):
#for i in range(0, len(dissimilar)-1, 2):
    #print(similar[i] + "-" + similar[i+1]) 
    posSource = dissimilar[i].find('.')
    posTarget = dissimilar[i+1].find('.')
    dissimilarSource = "cavitiesN/" + dissimilar[i][0:posSource] + "_dis" + str(int(i/2+1)) + "N.xyz"
    dissimilarTarget = "cavitiesN/" + dissimilar[i+1][0:posTarget] + "_dis" + str(int(i/2+1)) + "N.xyz"
    sourceDissimilarFile = open(dissimilarSource, "r")
    pointsNumber = sourceDissimilarFile.readline().replace("\n", "")
    outputFile = "output/dissimilar" + str(int(i/2+1)) + ".txt"
    sourceDissimilarFile.close()
    os.system(" ../GoICP " + dissimilarTarget + " " + dissimilarSource + " " + pointsNumber + " config.txt " + outputFile)