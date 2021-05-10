/********************************************************************
Main Function for point cloud registration with Go-ICP Algorithm
Last modified: Feb 13, 2014

"Go-ICP: Solving 3D Registration Efficiently and Globally Optimally"
Jiaolong Yang, Hongdong Li, Yunde Jia
International Conference on Computer Vision (ICCV), 2013

Copyright (C) 2013 Jiaolong Yang (BIT and ANU)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

/********************************************************************
Guillaume modifications:
- Get the cavities point clouds and normalize/centralize them 
- New config parameters
- Read c-FPFH files
- RMSD computation
*********************************************************************/

#include <time.h>
#include <iostream>
#include <fstream>
using namespace std;

#include "jly_goicp.h"
#include "ConfigMap.hpp"

#include "transformation.hpp"

#define DEFAULT_OUTPUT_FNAME "output.txt"
#define DEFAULT_CONFIG_FNAME "config.txt"
#define DEFAULT_MODEL_FNAME "model.txt"
#define DEFAULT_DATA_FNAME "data.txt"

bool useFPFH = false;

void parseInput(int argc, char **argv, string & modelFName, string & dataFName, int & NdDownsampled, string & configFName, string & outputFName, int & pair);
void readConfig(string FName, GoICP & goicp);
int loadPointCloud(string FName, int & N, POINT3D **  p);

int main(int argc, char** argv)
{
	int Nm, Nd, NdDownsampled;
	int pair = 1;
	clock_t  clockBegin, clockEnd;
	string modelFName, dataFName, configFName, outputFname;
	POINT3D * pModel, * pData;
	GoICP goicp;

	parseInput(argc, argv, modelFName, dataFName, NdDownsampled, configFName, outputFname, pair);
	readConfig(configFName, goicp);

	//////////////////////////////////////////////////////
	string proteinData = dataFName.substr(dataFName.find("/")+1, dataFName.find(".")-dataFName.find("/")-1);
	string proteinModel = modelFName.substr(modelFName.find("/")+1, modelFName.find(".")-modelFName.find("/")-1);
	
	Transformation t;
	//Source and target point clouds
    ifstream fileSource (dataFName);

    vector<point4D> cloudSource = t.readMolFile(fileSource);
    ifstream fileTarget (modelFName);
    vector<point4D> cloudTarget = t.readMolFile(fileTarget);
	//Normalization and centralization
    double xMean, yMean, zMean = 0;
    double xMeanTemp, yMeanTemp, zMeanTemp = 0;
    double xTrans, yTrans, zTrans = 0;
    double timeGO, error = 0;
    double rot[3][3] = {};
    double scaleSource = t.normalizeMolCloud(cloudSource, xMean, yMean, zMean);
    double scaleTarget = t.normalizeMolCloud(cloudTarget, xMeanTemp, yMeanTemp, zMeanTemp);
    double scale = scaleSource >= scaleTarget ? scaleSource : scaleTarget;
    t.scaleCloud(cloudSource, scale);
    t.scaleCloud(cloudTarget, scale);
    string sourceFileName = "cavitiesN/" + proteinData + "_sim" + to_string(pair) + "N.xyz";
    string targetFileName = "cavitiesN/" + proteinModel + "_sim" + to_string(pair) + "N.xyz";
    ofstream newFileSource(sourceFileName);
    ofstream newFileTarget(targetFileName);
    t.writeNormalizedMolCloudFile(newFileSource, cloudSource);
    t.writeNormalizedMolCloudFile(newFileTarget, cloudTarget);

	//////////////////////////////////////////////////////

	// Load model and data point clouds
	loadPointCloud(targetFileName, Nm, &pModel);
	loadPointCloud(sourceFileName, Nd, &pData);
	
	goicp.pModel = pModel;
	goicp.Nm = Nm;
	goicp.pData = pData;
	goicp.Nd = Nd;

	// Build Distance Transform
	cout << "Building Distance Transform..." << flush;
	clockBegin = clock();
	goicp.BuildDT();
	clockEnd = clock();
	cout << (double)(clockEnd - clockBegin)/CLOCKS_PER_SEC << "s (CPU)" << endl;

	// Run GO-ICP
	if(NdDownsampled > 0)
	{
		goicp.Nd = NdDownsampled; // Only use first NdDownsampled data points (assumes data points are randomly ordered)
	}
	cout << "Model ID: " << modelFName << " (" << goicp.Nm << "), Data ID: " << dataFName << " (" << goicp.Nd << ")" << endl;
	cout << "Registering..." << endl;
	clockBegin = clock();
	goicp.Register();
	clockEnd = clock();
	double time = (double)(clockEnd - clockBegin)/CLOCKS_PER_SEC;
	cout << "Optimal Rotation Matrix:" << endl;
	cout << goicp.optR << endl;
	cout << "Optimal Translation Vector:" << endl;
	cout << goicp.optT << endl;
	cout << "Finished in " << time << endl;
	cout << endl;

	ofstream ofile;
	ofile.open(outputFname.c_str(), ofstream::out);
    ofile << "Time: " << time << endl;
    ofile << "Rotation Matrix: " << endl << goicp.optR << endl;
    ofile << "Translation Vector: " << endl << goicp.optT << endl;
    ofile << "Error: " << goicp.optError << endl;
	//////////////////////////////////////////////////////
	//Add compatibilities to output file
	ofile << "Compatibilities: " << goicp.Nd-goicp.optComp << endl;
	//////////////////////////////////////////////////////
	ofile.close();

	//////////////////////////////////////////////////////
	//Application of transformation on cavity and rescaling
	//ofstream appliedFile("cavitiesA/" + proteinData + "_sim" + to_string(pair) + "A.mol2");
    	//t.writeAppliedMolFile(appliedFile, "output/similar" + to_string(pair) + ".txt", cloudSource.size(), cloudSource, rot, tra, xTrans, yTrans, zTrans, time, error);
    	xTrans = goicp.optT.val[0][0];
	yTrans = goicp.optT.val[1][0];
	zTrans = goicp.optT.val[2][0];
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			rot[i][j] = goicp.optR.val[i][j];
		}
	}
	ofstream rescaledFile(outputFname.substr(0, outputFname.find(".")) + "_rescaled.txt");
    	t.rescaleCloud(rescaledFile, cloudSource, scale, xMeanTemp, yMeanTemp, zMeanTemp, time, goicp.optError, xMean, yMean, zMean, rot, xTrans, yTrans, zTrans);

	//Application of transformation on protein and RMSD computation
	/*
	string protein = proteinData.substr(0, 6) + "_protein.mol2";
	ifstream proteinFile("chains/" + protein);
	ofstream transformedProteinFile("rot/rot_" + protein);
	t.applyTransformationProtein(transformedProteinFile, "chains/" + protein, pair);
	ifstream alignedFile("ref_proteins/" + proteinData.substr(0, 6) + "." + proteinModel.substr(0, 6) + "/aligned_" + protein);
	ifstream rotFile("rot/rot_" + protein);
	float rmsd = t.computeRMSD(alignedFile, rotFile);
	ofstream RMSDFile("resultsRMSD.txt", std::ios::app);
	if(RMSDFile.is_open()) {
		RMSDFile << pair << "\t" << proteinData.substr(0, 6) << "\t" << proteinModel.substr(0, 6) << "\t" << to_string(rmsd) << endl;
	}
	cout << "---> RMSD: " << rmsd << endl;
	*/
	//////////////////////////////////////////////////////
	
	delete(pModel);
	delete(pData);

	return 0;
}

void parseInput(int argc, char **argv, string & modelFName, string & dataFName, int & NdDownsampled, string & configFName, string & outputFName, int & pair)
{
	// Set default values
	modelFName = DEFAULT_MODEL_FNAME;
	dataFName = DEFAULT_DATA_FNAME;
	configFName = DEFAULT_CONFIG_FNAME;
	outputFName = DEFAULT_OUTPUT_FNAME;
	NdDownsampled = 0; // No downsampling

	//cout << endl;
	//cout << "USAGE:" << "./GOICP <MODEL FILENAME> <DATA FILENAME> <NUM DOWNSAMPLED DATA POINTS> <CONFIG FILENAME> <OUTPUT FILENAME>" << endl;
	//cout << endl;
	//////////////////////////////////////////////////////
	//New program argument with the pair number
	if(argc > 6) 
	{
		pair = atoi(argv[6]);
	}
	//////////////////////////////////////////////////////
 	if(argc > 5)
	{
		outputFName = argv[5];
	}
	if(argc > 4)
	{
		configFName = argv[4];
	}
	if(argc > 3)
	{
		NdDownsampled = atoi(argv[3]);
	}
	if(argc > 2)
	{
		dataFName = argv[2];
	}
	if(argc > 1)
	{
		modelFName = argv[1];
	}

	cout << "INPUT:" << endl;
	cout << "(modelFName)->(" << modelFName << ")" << endl;
	cout << "(dataFName)->(" << dataFName << ")" << endl;
	cout << "(NdDownsampled)->(" << NdDownsampled << ")" << endl;
	cout << "(configFName)->(" << configFName << ")" << endl;
	cout << "(outputFName)->(" << outputFName << ")" << endl;
	cout << "(pair)->(" << pair << ")" << endl;
	cout << endl;
}

void readConfig(string FName, GoICP & goicp)
{
	// Open and parse the associated config file
	ConfigMap config(FName.c_str());

	goicp.MSEThresh = config.getF("MSEThresh");
	goicp.initNodeRot.a = config.getF("rotMinX");
	goicp.initNodeRot.b = config.getF("rotMinY");
	goicp.initNodeRot.c = config.getF("rotMinZ");
	goicp.initNodeRot.w = config.getF("rotWidth");
	goicp.initNodeTrans.x = config.getF("transMinX");
	goicp.initNodeTrans.y = config.getF("transMinY");
	goicp.initNodeTrans.z = config.getF("transMinZ");
	goicp.initNodeTrans.w = config.getF("transWidth");
	goicp.trimFraction = config.getF("trimFraction");

	//////////////////////////////////////////////////////
	//New parameters
	goicp.regularization = config.getF("regularization");
	goicp.regularizationNeighbors = config.getF("regularizationNeighbors");
	goicp.regularizationFPFH = config.getF("regularizationFPFH");
	goicp.cfpfh = config.getI("cfpfh");
	if (goicp.cfpfh != 0) useFPFH = true;
	goicp.norm = config.getI("norm");
	goicp.ponderation = config.getI("ponderation");
	//////////////////////////////////////////////////////

	// If < 0.1% trimming specified, do no trimming
	if(goicp.trimFraction < 0.001)
	{
		goicp.doTrim = false;
	}
	goicp.dt.SIZE = config.getI("distTransSize");
	goicp.dt.expandFactor = config.getF("distTransExpandFactor");

	cout << "CONFIG:" << endl;
	config.print();
	//cout << "(doTrim)->(" << goicp.doTrim << ")" << endl;
	cout << endl;
}

int loadPointCloud(string FName, int & N, POINT3D ** p)
{
	int i;
	//////////////////////////////////////////////////////
	//c-FPFH files
	ifstream ifile, cfpfhfile;

	string fileName = "cfpfh/" + FName.substr(FName.find("/")+1, FName.find_last_of("_") - FName.find("/") - 1) + ".cfpfh";
	cfpfhfile.open(fileName.c_str(), ifstream::in);
	//////////////////////////////////////////////////////
	ifile.open(FName.c_str(), ifstream::in);

	if(!ifile.is_open())
	{
		cout << "Unable to open point file '" << FName << "'" << endl;
		exit(-1);
	}
	ifile >> N; // First line has number of points to follow
	*p = (POINT3D *)malloc(sizeof(POINT3D) * N);

	//////////////////////////////////////////////////////
	//read c-FPFH files
	if(!cfpfhfile.is_open())
	{
		cout << "Unable to open fpfh file '" << fileName << "'" << endl;
		exit(-1);
	}
	for(i = 0; i < N; i++)
	{
        ifile >> (*p)[i].x >> (*p)[i].y >> (*p)[i].z >> (*p)[i].c;
		
		for(int j = 0; j < 41; j++) {
			float bin = 0;
			cfpfhfile >> bin;
			(*p)[i].cfpfh.push_back(bin);
		}
	}
	cfpfhfile.close();
	//////////////////////////////////////////////////////
	ifile.close();

	return 0;
}