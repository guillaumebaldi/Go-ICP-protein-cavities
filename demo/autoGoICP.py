import os

similar = []
dissimilar = []

readmeFile = open("README2.txt", "r")
count = 0
while True:
    line = readmeFile.readline()
    if (count > 10):
        if not line or line.isspace():
            break    
        index = line.find('\t')
        source = line[0:index].replace(" ","").replace("\t", "").replace("\n", "")
        target = line[index+1:].replace(" ", "").replace("\t", "").replace("\n", "")
        similar.append(source)
        similar.append(target)
    count += 1
    
readmeFile.readline()

while True:
    line = readmeFile.readline()
    if not line or line.isspace():
        break    
    index = line.find('\t')
    source = line[0:index].replace(" ","").replace("\t", "").replace("\n", "")
    target = line[index+1:].replace(" ", "").replace("\t", "").replace("\n", "")
    dissimilar.append(source)
    dissimilar.append(target)

readmeFile.close()

for i in range(0, len(similar)-1, 2):
#for i in range(10, len(similar)-1, 2):
    #if(i == 2 or i == 4) : continue
    print(similar[i] + "-" + similar[i+1]) 
    posSource = similar[i].find('.')
    posTarget = similar[i+1].find('.')
    similarSource = similar[i][0:posSource] + "_sim" + str(int(i/2+1)) + "N" + similar[i][posSource:]
    similarTarget = similar[i+1][0:posTarget] + "_sim" + str(int(i/2+1)) + "N" + similar[i+1][posTarget:]
    sourceSimilarFile = open(similarSource, "r")
    pointsNumber = sourceSimilarFile.readline().replace("\n", "")
    outputFile = "Esimilar" + str(int(i/2+1)) + ".txt"
    sourceSimilarFile.close()
    os.system(" ../GoICP " + similarTarget + " " + similarSource + " " + pointsNumber + " config.txt " + outputFile)


#for i in range(0, len(dissimilar)-1, 2):
for i in range(0, 0, 2):
    print(dissimilar[i] + "-" + dissimilar[i+1])
    posSource = dissimilar[i].find('.')
    posTarget = dissimilar[i+1].find('.')
    dissimilarSource = dissimilar[i][0:posSource] + "_dis" + str(int(i/2+1)) + "N" + dissimilar[i][posSource:]
    dissimilarTarget = dissimilar[i+1][0:posTarget] + "_dis" + str(int(i/2+1)) + "N" + dissimilar[i+1][posTarget:]
    sourceDissimilarFile = open(dissimilarSource, "r")
    pointsNumber = sourceDissimilarFile.readline().replace("\n", "")
    outputFile = "dissimilar" + str(int(i/2+1)) + ".txt"
    sourceDissimilarFile.close()
    os.system(" ../../build-Go-ICP-master-Desktop-Default/GoICP " + dissimilarTarget + " " + dissimilarSource + " " + pointsNumber + " config.txt " + outputFile)


