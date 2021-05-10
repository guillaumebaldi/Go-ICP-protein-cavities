/********************************************************************
Transformation class
- Read xyz point clouds files, normalize/centralize them, apply transformation
- Read xyzc point clouds files, normalize/centralize them, apply transformation
- Apply transformation on protein and RMSD computation
*********************************************************************/


#include "transformation.hpp"

Transformation::Transformation() {
    
}

/**
 * Physico-chemical property conversion between color code and name
 */
properties Transformation::string_to_prop(string type) {
    if(type == "OG") {
        return OG;
    }
    else if(type == "N") {
        return N;
    }
    else if(type == "O") {
        return O;
    }
    else if(type == "NZ") {
        return NZ;
    }
    else if(type == "CZ") {
        return CZ;
    }
    else if(type == "CA") {
        return CA;
    }
    else if(type == "DU") {
        return DU;
    }
    else if(type == "OD1") {
        return OD1;
    }
    else if(type == "C") {
        return C;
    }
    return OG;
}

/**
 * Read xyz point cloud
 */
vector<point> Transformation::readCloudFile(ifstream &file) {
    vector<point> cloud {};
    string line;
    if(file.is_open()) {
        point p {};
        file >> p.x >> p.y >> p.z;
        cloud.push_back(p);
        while(getline(file, line)) {
            point p {};
            file >> p.x >> p.y >> p.z;
            cloud.push_back(p);
        }
        cloud.pop_back();
    }
    file.close();
    return cloud;
}

/**
 * Normalize xyz point cloud
 */
void Transformation::normalizeCloud(vector<point> &cloud) {
    unsigned long size = cloud.size();
    double xMean, yMean, zMean = 0;
    for(point &p : cloud) {
        xMean += p.x;
        yMean += p.y;
        zMean += p.z;
    }
    xMean /= size;
    yMean /= size;
    zMean /= size;
    double maxNorm = 0;
    for(point &p : cloud) {
        p.x -= xMean;
        p.y -= yMean;
        p.z -= zMean;
        double norm = sqrt(pow(p.x,2) + pow(p.y,2) + pow(p.z,2));
        if(norm > maxNorm) {
            maxNorm = norm;
        }
    }
    for(point &p : cloud) {
        p.x /= maxNorm;
        p.y /= maxNorm;
        p.z /= maxNorm;
    }
}

/**
 * Write normalized xyz point cloud
 */
void Transformation::writeNormalizedCloudFile(ofstream &file, vector<point> &cloud) {
    if(file.is_open()) {
        normalizeCloud(cloud);
        file << cloud.size() << endl;
        for(point p : cloud) {
            ostringstream newLine;
            newLine << p.x << " " << p.y << " " << p.z << endl;
            file << newLine.str();
        }
    }
    file.close();
}

/**
 * Read output from Go-ICP
 */
void Transformation::readOutput(ifstream &file, double rot[3][3], double tra[3][1], double &time, double &error) {
    if(file.is_open()) {
        string s;
        file >> s;
        file >> time;
        file >> s >> s;
        for (int i = 0; i < 3; i++) {
            file >> rot[i][0] >> rot[i][1] >> rot[i][2];
        }
        file >> s >> s;
        for (int i = 0; i < 3; i++) {
            file >> tra[i][0];
        }
        file >> s;
        file >> error;
        file >> s;
        string comp;
        file >> comp;
    }
}

/**
 * Apply transformation on xyz point cloud
 */
vector<point> Transformation::applyTransformation(vector<point> cloud, string nameOutputFile, int number) {
    double rot[3][3] {};
    double tra[3][1] {};
    ifstream output (nameOutputFile);
    vector<point> newCloud;

    for (int i = 0; i < number; i++) {
        point p2 {};
        p2.x = rot[0][0] * cloud[i].x + rot[0][1] * cloud[i].y + rot[0][2] * cloud[i].z;
        p2.y = rot[1][0] * cloud[i].x + rot[1][1] * cloud[i].y + rot[1][2] * cloud[i].z;
        p2.z = rot[2][0] * cloud[i].x + rot[2][1] * cloud[i].y + rot[2][2] * cloud[i].z;
        p2.x += tra[0][0];
        p2.y += tra[1][0];
        p2.z += tra[2][0];
        newCloud.push_back(p2);
    }

    return newCloud;
}

/**
 * Get applied xyz point cloud
 */
void Transformation::writeAppliedFile(ofstream &file, string nameOutputFile, int number, vector<point> cloud) {
    vector<point> newCloud = applyTransformation(cloud, nameOutputFile, number);
    if(file.is_open()) {
        for(point p : newCloud) {
            ostringstream newLine;
            file << p.x << " " << p.y << " " << p.z << endl;
            file << newLine.str();
        }
    }
    file.close();
}

void Transformation::writePoints(ofstream &file, int number, vector<point> cloud) {
    if(file.is_open()) {
        for (int i = 0; i < number; i++) {
            ostringstream newLine;
            file << cloud[i].x << " " << cloud[i].y << " " << cloud[i].z << endl;
            file << newLine.str();
        }
    }
    file.close();
}


/**
 * Read .pcd point cloud names from a readme file
 */
void Transformation::readConfigProteinFile(ifstream &file, vector<string> &similar, vector<string> &dissimilar) {
    string line;
    if(file.is_open()) {
        for (int i = 0;i < 11; i++) {
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        while(getline(file, line)) {
            if(line.empty()) {
                break;
            }
            size_t posSource = line.find("\t");
            string s1 = line.substr(0, posSource);
            string s2 = line.substr(posSource);
            s1.erase(std::remove(s1.begin(),s1.end(),'\t'),s1.end());
            s1.erase(std::remove(s1.begin(),s1.end(),' '),s1.end());
            s2.erase(std::remove(s2.begin(),s2.end(),'\t'),s2.end());
            s2.erase(std::remove(s2.begin(),s2.end(),' '),s2.end());
            similar.push_back(s1);
            similar.push_back(s2);
        }
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        while(getline(file, line)) {
            if(line.empty()) {
                break;
            }
            size_t posSource = line.find("\t");
            string s1 = line.substr(0, posSource);
            string s2 = line.substr(posSource);
            s1.erase(std::remove(s1.begin(),s1.end(),'\t'),s1.end());
            s1.erase(std::remove(s1.begin(),s1.end(),' '),s1.end());
            s2.erase(std::remove(s2.begin(),s2.end(),'\t'),s2.end());
            s2.erase(std::remove(s2.begin(),s2.end(),' '),s2.end());
            dissimilar.push_back(s1);
            dissimilar.push_back(s2);
        }
    }
    file.close();
}

/**
 * Read xyzc point cloud from a .pcd file
 */
vector<point4D> Transformation::readPCDfile(ifstream &file) {
    vector<point4D> cloud {};
    string line;
    if(file.is_open()) {
        for (int i = 0;i < 10; i++) {
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        point4D p {};
        file >> p.x >> p.y >> p.z >> p.c;
        cloud.push_back(p);
        while(getline(file, line)) {
            point4D p {};
            file >> p.x >> p.y >> p.z >> p.c;
            cloud.push_back(p);
        }
    }
    file.close();
    return cloud;
}

/**
 * Read .mol point cloud names from a readme file
 */
void Transformation::readConfigMolFile(ifstream &file, vector<string> &cavities) {
    string line;
    if(file.is_open()) {
        while(getline(file, line)) {
            if(line.empty()) {
                break;
            }
            string delimiter = "\t";
            string temp = line.substr(line.find(delimiter)+1);
            temp = temp.substr(temp.find(delimiter)+1);
            string s1 = temp.substr(0, temp.find(delimiter)) + "_cavity6.mol2";
            temp = temp.substr(temp.find(delimiter)+1);
            string s2 = temp.substr(0, temp.find(delimiter)) + "_cavity6.mol2";;
            cavities.push_back(s1);
            cavities.push_back(s2);
        }
    }
    file.close();
}

/**
 * Read xyzc point cloud from a .mol file
 */
vector<point4D> Transformation::readMolFile(ifstream &file) {
    vector<point4D> cloud {};
    string line;
    bool header = true;
    string delimiter = "\t";
    if(file.is_open()) {
        while(getline(file, line)) {
	    if (line.find("@<TRIPOS>ATOM") != std::string::npos) {
    		header = false;
	    }
            if(header) {
                continue;
            }
            string a = "";
            string b = "";
            point4D p {};
            file >> a >> b >> p.x >> p.y >> p.z;
            p.c = string_to_prop(b);
            cloud.push_back(p);
        }
        cloud.pop_back();
    }

    return cloud;
}

/**
 * Normalize xyzc point cloud
 */
double Transformation::normalizeMolCloud(vector<point4D> &cloud, double &xMean, double &yMean, double &zMean) {
    unsigned long size = cloud.size();
    xMean = 0;
    yMean = 0;
    zMean = 0;
    for(point4D &p : cloud) {
        xMean += p.x;
        yMean += p.y;
        zMean += p.z;
    }
    xMean /= size;
    yMean /= size;
    zMean /= size;
    double maxNorm = 0;
    for(point4D &p : cloud) {
        p.x -= xMean;
        p.y -= yMean;
        p.z -= zMean;
        double norm = sqrt(pow(p.x,2) + pow(p.y,2) + pow(p.z,2));
        if(norm > maxNorm) {
            maxNorm = norm;
        }
    }
    return maxNorm;
}

/**
 * Write normalized xyzc point cloud
 */
void Transformation::writeNormalizedMolCloudFile(ofstream &file, vector<point4D> &cloud) {
    if(file.is_open()) {
        file << cloud.size() << endl;
        for(point4D p : cloud) {
            ostringstream newLine;
            newLine << p.x << " " << p.y << " " << p.z << " " << p.c << endl;
            file << newLine.str();
        }
    }
    file.close();
}

/**
 * Scale xyzc point cloud
 */
void Transformation::scaleCloud(vector<point4D> &cloud, double scale) {
    for(point4D &p : cloud) {
        p.x /= scale;
        p.y /= scale;
        p.z /= scale;
    }
}

/**
 * Apply transformation on xyzc point cloud
 */
void Transformation::applyMolTransformation(vector<point4D> &cloud, string nameOutputFile, int number, double rot[3][3], double tra[3][1], double &xTrans, double &yTrans, double &zTrans, double &time, double &error) {
    ifstream output (nameOutputFile);
    readOutput(output, rot, tra, time, error);
    xTrans = tra[0][0];
    yTrans = tra[1][0];
    zTrans = tra[2][0];
    for (int i = 0; i < number; i++) {
        double x = cloud[i].x;
        double y = cloud[i].y;
        double z = cloud[i].z;
        cloud[i].x = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z;
        cloud[i].y = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z;
        cloud[i].z = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z;
        cloud[i].x += tra[0][0];
        cloud[i].y += tra[1][0];
        cloud[i].z += tra[2][0];
    }
}

/**
 * Write applied xyzc point cloud
 */
void Transformation::writeAppliedMolFile(ofstream &file, string nameOutputFile, int number, vector<point4D> &cloud, double rot[3][3], double tra[3][1], double &xTrans, double &yTrans, double &zTrans, double &time, double &error) {
    applyMolTransformation(cloud, nameOutputFile, number, rot, tra, xTrans, yTrans, zTrans, time, error);
    if(file.is_open()) {
        for(point4D p : cloud) {
            ostringstream newLine;
            file << p.x << " " << p.y << " " << p.z << " " << p.c << endl;
            file << newLine.str();
        }
    }
    file.close();
}

/**
 * Write rescaled xyzc point cloud
 */
void Transformation::rescaleCloud(ofstream &file, vector<point4D> &cloud, double scale, double xMeanTemp, double yMeanTemp, double zMeanTemp, double time, double error, double xMean, double yMean, double zMean, double rot[3][3], double xTrans, double yTrans, double zTrans) {
    if(file.is_open()) {
        file << "Time: " << time << endl;
        file << "Rotation Matrix:" << endl << "   " << rot[0][0] << "   " << rot[0][1] << "   " << rot[0][2] << endl;
        file << "   " << rot[1][0] << "   " << rot[1][1] << "   " << rot[1][2] << endl;
        file << "   " << rot[2][0] << "   " << rot[2][1] << "   " << rot[2][2] << endl;
        double x, y, z = 0;
        x = -(rot[0][0] * xMean + rot[0][1] * yMean + rot[0][2] * zMean) + (scale * xTrans) + xMeanTemp;
        y = -(rot[1][0] * xMean + rot[1][1] * yMean + rot[1][2] * zMean) + (scale * yTrans) + yMeanTemp;
        z = -(rot[2][0] * xMean + rot[2][1] * yMean + rot[2][2] * zMean) + (scale * zTrans) + zMeanTemp;
        file << "Translation Vector:" << endl << "   " << x << endl << "   " << y << endl << "   " << z << endl;
        file << "Error: " << error << endl;
    }
    file.close();
}


/**
 * Get points with atoms concerned by RMSD computation
 */
vector<point4D> Transformation::getAtomBlock(ifstream &file) {
    vector<point4D> cloud {};
    string line;
    bool header = true;
    string delimiter = "\t";
    if(file.is_open()) {
        while(getline(file, line)) {
            if(line == "@<TRIPOS>ATOM") {
                    header = false;
            }
            if(header) {
                continue;
            }
            string a = "";
            string b = "";
            point4D p {};
            file >> a >> b >> p.x >> p.y >> p.z;
            p.c = string_to_prop(b);
            if(p.c == C || p.c == CA || p.c == N || p.c == O) {
                cloud.push_back(p);
            } 
        }
    }

    return cloud;
}

/**
 * RMSD computation
 */
float Transformation::computeRMSD(ifstream &alignedFile, ifstream &rotFile) {
    float rmsd = 0;
    vector<point4D> alignedPoints = getAtomBlock(alignedFile);
    vector<point4D> transformedPoints = getAtomBlock(rotFile);

    for(int i = 0; i < alignedPoints.size(); i++) {
        rmsd += pow((alignedPoints[i].x - transformedPoints[i].x), 2) + pow((alignedPoints[i].y - transformedPoints[i].y), 2) + pow((alignedPoints[i].z - transformedPoints[i].z), 2);
    }
    rmsd = sqrt(rmsd / alignedPoints.size());

    return rmsd;
}

/**
 * Apply transformation on protein file
 */
void Transformation::applyTransformationProtein(ofstream &output, string protein, int pair) {
    ifstream outputFile ("cavitiesR/similar" + to_string(pair) + ".txt");
    double time, error = 0;
    double rot[3][3] = {};
    double tra[3][1] = {};
    readOutput(outputFile, rot, tra, time, error);

    ifstream file(protein);
    vector<point4D> points = readMolFile(file);
    vector<point4D> appliedPoints;
    for(point4D p : points) {
        point4D np;
        double x = p.x;
        double y = p.y;
        double z = p.z;
        np.x = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z;
        np.y = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z;
        np.z = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z;
        np.x += tra[0][0];
        np.y += tra[1][0];
        np.z += tra[2][0];
        appliedPoints.push_back(np);
    }
    string line;
    int count = 0;
    bool header = true;
    bool end = false;
    ifstream proteinFile(protein);
    if(proteinFile.is_open()) {
        while(getline(proteinFile, line)) {
            if(line == "@<TRIPOS>ATOM") {
                output << line << endl;
                header = false;
            }
            if(header) {
                output << line << endl;
                continue;
            }
            if(end) {
                output << line << endl;
            }
            else {
                string s1 = "";
                string s2 = "";
                string s3 = "";
                string s4 = "";
                string s5 = "";
                string s6 = "";
                string s7 = "";
                string s8 = "";
                string s9 = "";
                proteinFile >> s1;
                if(s1 == "@<TRIPOS>BOND") {
                    output << s1 << endl;
                    end = true;
                }
                if(!end) {
                    proteinFile >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9;
                    s3 = to_string(appliedPoints[count].x);
                    s4 = to_string(appliedPoints[count].y);
                    s5 = to_string(appliedPoints[count].z);
                    output << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\t" << s5 << "\t" << s6 << "\t" << s7 << "\t" << s8 << "\t" << s9 << endl;
                    count++;
                } 
            }
        }
    }
    
    proteinFile.close();
    output.close();
}