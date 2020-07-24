/********************************************************************
Header File for Transformation class
*********************************************************************/

#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP

#include <iostream>
#include "math.h"
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <random>
#include <map>
#include <iterator>
#include "matrix.h"

using namespace std;

struct point {
    double x, y, z;
};

struct point4D {
    double x, y, z;
    int c;

    bool operator<(const point4D& a) const
    {
        return y < a.y;
    }
};

enum properties {OG = 8204959, N = 30894, O = 15219528, NZ = 15231913, CZ = 4646984, CA = 16741671, DU = 7566712, OD1 = 0, C = 1 };

class Transformation {
public:
    Transformation();
    properties string_to_prop(string type);

    //Fonctions on xyz point clouds files
    vector<point> readCloudFile(ifstream &file);
    void normalizeCloud(vector<point> &cloud);
    void writeNormalizedCloudFile(ofstream &file, vector<point> &cloud);
    void readOutput(ifstream &file, double rot[3][3], double tra[3][1], double &time, double &error);
    vector<point> applyTransformation(vector<point> cloud, string nameOutputFile, int number);
    void writeAppliedFile(ofstream &file, string nameOutputFile, int number, vector<point> cloud);
    void writePoints(ofstream &file, int number, vector<point> cloud);
    
    //Fonctions on xyzc point clouds files (cavities)
    void readConfigProteinFile(ifstream &file, vector<string> &similar, vector<string> &dissimilar);
    vector<point4D> readPCDfile(ifstream &file);
    void readConfigMolFile(ifstream &file, vector<string> &cavities);
    vector<point4D> readMolFile(ifstream &file);
    double normalizeMolCloud(vector<point4D> &cloud, double &xMean, double &yMean, double &zMean);
    void writeNormalizedMolCloudFile(ofstream &file, vector<point4D> &cloud);
    void scaleCloud(vector<point4D> &cloud, double scale);
    void applyMolTransformation(vector<point4D> &cloud, string nameOutputFile, int number, double rot[3][3], double tra[3][1], double &xTrans, double &yTrans, double &zTrans, double &time, double &error);
    void writeAppliedMolFile(ofstream &file, string nameOutputFile, int number, vector<point4D> &cloud, double rot[3][3], double tra[3][1], double &xTrans, double &yTrans, double &zTrans, double &time, double &error);
    void rescaleCloud(ofstream &file, vector<point4D> &cloud, double scale, double xMeanTemp, double yMeanTemp, double zMeanTemp, double time, double error, double xMean, double yMean, double zMean, double rot[3][3], double xTrans, double yTrans, double zTrans);
    
    //Fonctions of transformation application on protein and RMSD computation
    vector<point4D> getAtomBlock(ifstream &file);
    float computeRMSD(ifstream &alignedFile, ifstream &rotFile);
    void applyTransformationProtein(ofstream &output, string protein, int pair);
};

#endif