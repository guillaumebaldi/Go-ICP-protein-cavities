/********************************************************************
Header File for Go-ICP Class
Last modified: Apr 21, 2014

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
- Add of property, neighbors number, density and c-FPFH to POINT3D structure
- Structures TRANSCOMPATIBILITIES and TRANSFPFH to keep in memory already tested translation for incompatibilities number and FPFH difference
- Fonctions and variables on parameters
*********************************************************************/

#ifndef JLY_GOICP_H
#define JLY_GOICP_H

#include <queue>
using namespace std;

#include "jly_icp3d.hpp"
#include "jly_3ddt.h"
#include <map>
#include "transformation.hpp"
//#include <eigen3/Eigen/Eigenvalues>

#define PI 3.1415926536
#define SQRT3 1.732050808

typedef struct _POINT3D
{
    float x, y, z;
	//////////////////////////////////////////////////////
    int c;
	int neighbors;
	float density;
	vector<float> cfpfh;
	//////////////////////////////////////////////////////
}POINT3D;

typedef struct _ROTNODE
{
	float a, b, c, w;
	float ub, lb;
	int l;
	friend bool operator < (const struct _ROTNODE & n1, const struct _ROTNODE & n2)
	{
		if(n1.lb != n2.lb)
			return n1.lb > n2.lb;
		else
			return n1.w < n2.w;
			//return n1.ub > n2.ub;
	}
	
}ROTNODE;

typedef struct _TRANSNODE
{
	float x, y, z, w;
	float ub, lb;
	friend bool operator < (const struct _TRANSNODE & n1, const struct _TRANSNODE & n2)
	{
		if(n1.lb != n2.lb)
			return n1.lb > n2.lb;
		else
			return n1.w < n2.w;
			//return n1.ub > n2.ub;
	}
}TRANSNODE;

/********************************************************/



/********************************************************/

#define MAXROTLEVEL 20
#define ROUND(x) (int((x)+0.5))

//////////////////////////////////////////////////////
struct TRANSCOMPATIBILITIES
{
    float x, y ,z;
	int c;
};

struct TRANSFPFH
{
    float x, y ,z;
	float fpfh;
};
//////////////////////////////////////////////////////

class GoICP
{
public:
	int Nm, Nd;
	POINT3D * pModel, * pData;

	ROTNODE initNodeRot;
	TRANSNODE initNodeTrans;

	DT3D dt;

	ROTNODE optNodeRot;
	TRANSNODE optNodeTrans;

	GoICP();
	float Register();
	void BuildDT();

	float MSEThresh;
	float SSEThresh;
	float icpThresh;

	float optError;
	Matrix optR;
	Matrix optT;

	clock_t clockBegin;
    ICP3D<float> icp3d;
	float trimFraction;
	int inlierNum;
	bool doTrim;

	//////////////////////////////////////////////////////
	float regularization;
	float regularizationNeighbors;
	float regularizationFPFH;
	int optComp;
	int norm;
	int ponderation;
	int cfpfh;
	//int * compBNB;
	//////////////////////////////////////////////////////

private:
	//temp variables
	float * normData;
	float * minDis;
	float** maxRotDis;
	float * maxRotDisL;
	POINT3D * pDataTemp;
	POINT3D * pDataTempICP;

	float * M_icp;
	float * D_icp;

	float ICP(Matrix& R_icp, Matrix& t_icp);
	float InnerBnB(float* maxRotDisL, TRANSNODE* nodeTransOut);
	float OuterBnB();
	void Initialize();
	void Clear();

	//////////////////////////////////////////////////////
	//Fonctions on physico-chemical properties incompatibilities
    int countCompatibilities(bool init);
	int checkCompatibilities(float x, float y, float z);
    void updateCompatibilities();
    map<properties, vector<properties>> compatibilities;
    void assignCellColor();
    bool checkCompatibility(int i, float x, float y, float z, bool update);
    void closestEmptyCell(int x, int y, int z);
	void sortDistances(); 
	bool checkProperty(properties source, properties target, int cx, int cy, int cz, double x, double y, double z);

	//Fonctions on direct neighbors
	bool isNeighbor(float radius, POINT3D p, POINT3D p2);
	double calculateMean(vector<POINT3D> points, char d);
	void centralizePoints(vector<POINT3D> &points, double &xMean, double &yMean, double &zMean);
	void calculateCovarianceMatrix(vector<POINT3D> points, double matrix[9]);
	void calculateEigen(double matrix[9], vector<double> &eigenvalues);
	double computePlanarity(double lambda1, double lambda2, double lambda3);
	double computeScattering(double lambda1, double lambda3);
	int nearestNeighbor(float x, float y, float z, int cx, int cy, int cz);
	void assignNeighbors();
	int compareNeighbors(bool icp, float x, float y, float z);
	int compareNeighborsV2(bool icp, float x, float y, float z);
	int compareNeighborsV3(bool icp, float x, float y, float z);

	//Fonctions on neighbors ponderation and density
	void assignNeighborsPrio();
	void neighborsWeights();
	void neighborsDensity();
	float computeDensityDifference(bool icp, int point, float x, float y, float z);
	float sumDensities(float x, float y, float z);

	//Fonctions on c-FPFH
	float computeFPFHDifference(bool icp, int point, float x, float y, float z);
	float sumFPFH(float x, float y, float z);

	float * weights;
	//false : BNB is the last operation which decreased the error; true : ICP is the last operation which decreased the error
	bool lastICP = false;
	int maxComp;
	//Error decomposition
	float geomError, incompError, FPFHError = 0;
	float geomErrorTemp, incompErrorTemp, FPFHErrorTemp = 0;
	//////////////////////////////////////////////////////
};


/********************************************************/

#endif
