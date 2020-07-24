/********************************************************************
Implementation of Go-ICP Algorithm
Last modified: Jun 18, 2014

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
- Modifications to ICP and BnB to add parameters
- Fonctions on parameters
*********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string>


//using namespace std;

#include "jly_goicp.h"
#include "jly_sorting.hpp"

GoICP::GoICP()
{
	initNodeRot.a = -PI;
	initNodeRot.b = -PI;
	initNodeRot.c = -PI;
	initNodeRot.w = 2*PI;
	initNodeRot.l = 0;

	initNodeRot.lb = 0;
	initNodeTrans.lb = 0;

	doTrim = true;

	//////////////////////////////////////////////////////
	//Compatibilities between physico-chemical properties
    /*compatibilities.insert(make_pair(CA, vector<properties>{CA, CZ}));
    compatibilities.insert(make_pair(CZ, vector<properties>{CZ, CA}));
    compatibilities.insert(make_pair(N, vector<properties>{N, NZ, OG}));
    compatibilities.insert(make_pair(NZ, vector<properties>{NZ, N, OG}));
    compatibilities.insert(make_pair(O, vector<properties>{O, OD1, OG}));
    compatibilities.insert(make_pair(OD1, vector<properties>{OD1, O, OG}));
    compatibilities.insert(make_pair(OG, vector<properties>{OG, N, O, OD1, NZ}));
    compatibilities.insert(make_pair(DU, vector<properties>{DU}));*/
	compatibilities.insert(make_pair(CA, vector<properties>{CA}));
    compatibilities.insert(make_pair(CZ, vector<properties>{CZ}));
    compatibilities.insert(make_pair(N, vector<properties>{N}));
    compatibilities.insert(make_pair(NZ, vector<properties>{NZ}));
    compatibilities.insert(make_pair(O, vector<properties>{O}));
    compatibilities.insert(make_pair(OD1, vector<properties>{OD1}));
    compatibilities.insert(make_pair(OG, vector<properties>{OG}));
    compatibilities.insert(make_pair(DU, vector<properties>{DU}));
	//compBNB = (int*)malloc(sizeof(int)*Nd);
	//////////////////////////////////////////////////////
}

// Build Distance Transform
void GoICP::BuildDT()
{
	double* x = (double*)malloc(sizeof(double)*Nm);
	double* y = (double*)malloc(sizeof(double)*Nm);
	double* z = (double*)malloc(sizeof(double)*Nm);
	for(int i = 0; i < Nm; i++)
	{
		x[i] = pModel[i].x;
		y[i] = pModel[i].y;
		z[i] = pModel[i].z;
	}
	dt.Build(x, y, z, Nm);
	//////////////////////////////////////////////////////
	//Property of DT cells and direct neighbors of points
    assignCellColor();
	if(regularizationNeighbors > 0) assignNeighbors();
	//////////////////////////////////////////////////////
	delete(x);
	delete(y);
	delete(z);
}

// Run ICP and calculate sum squared L2 error
float GoICP::ICP(Matrix& R_icp, Matrix& t_icp)
{
	int i;
	float error, dis;

	icp3d.Run(D_icp, Nd, R_icp, t_icp); // data cloud, # data points, rotation matrix, translation matrix

	// Transform point cloud and use DT to determine the L2 error
	error = 0;

	//////////////////////////////////////////////////////
	//float dens = 0;
	float fpfh = 0;
	//////////////////////////////////////////////////////

	for(i = 0; i < Nd; i++)
	{
		POINT3D& p = pData[i];
		pDataTempICP[i].x = R_icp.val[0][0]*p.x+R_icp.val[0][1]*p.y+R_icp.val[0][2]*p.z + t_icp.val[0][0];
		pDataTempICP[i].y = R_icp.val[1][0]*p.x+R_icp.val[1][1]*p.y+R_icp.val[1][2]*p.z + t_icp.val[1][0];
		pDataTempICP[i].z = R_icp.val[2][0]*p.x+R_icp.val[2][1]*p.y+R_icp.val[2][2]*p.z + t_icp.val[2][0];
        int x, y, z = 0;
		if(!doTrim)
		{
			//////////////////////////////////////////////////////
			//Add weights and norm consideration
            dis = weights[i] * dt.Distance(pDataTempICP[i].x, pDataTempICP[i].y, pDataTempICP[i].z, x, y, z);
			if(norm == 2) error += dis * dis;
            if(norm == 1) error += dis;
			//////////////////////////////////////////////////////
		}
		else
		{
            minDis[i] = dt.Distance(pDataTempICP[i].x, pDataTempICP[i].y, pDataTempICP[i].z, x, y, z);
		}

		//////////////////////////////////////////////////////
		//c-FPFH and density
		//dens += computeDensityDifference(true, i, 0, 0 ,0);
		if(cfpfh != 0) fpfh += computeFPFHDifference(true, i, 0, 0 ,0);
		//error = error * computeDensityDifference(true, i, 0, 0 ,0);
		//////////////////////////////////////////////////////
	}
	//////////////////////////////////////////////////////
	//c-FPFH computation and add of parameters to error function according to their regularization
	fpfh = fpfh / Nd;
	geomErrorTemp = error;
	if(regularizationNeighbors > 0) {
		int neighbors = compareNeighbors(true, 0, 0, 0);
		cout << "ICP : " << neighbors << endl;
		error += regularizationNeighbors * (neighbors*neighbors);
	}
	if(regularization > 0) {
		int incomp = countCompatibilities(true);
    	error += regularization * (incomp*incomp);
		incompErrorTemp = abs(error - geomErrorTemp);
		//error += regularization * (dens * dens);
	}
	if(regularizationFPFH > 0) {
		error += regularizationFPFH * (fpfh * fpfh);
		FPFHErrorTemp = abs(error - incompErrorTemp - geomErrorTemp);
	}
	//////////////////////////////////////////////////////

	if(doTrim)
	{
		//qsort(minDis, Nd, sizeof(float), cmp);
		//myqsort(minDis, Nd, inlierNum);
		intro_select(minDis,0,Nd-1,inlierNum-1);
		for(i = 0; i < inlierNum; i++)
		{
			error += minDis[i]*minDis[i];
		}
	}
	cout << endl;
	return error;
}

void GoICP::Initialize()
{
	int i, j;
	float sigma, maxAngle;

	// Precompute the rotation uncertainty distance (maxRotDis) for each point in the data and each level of rotation subcube

	// Calculate L2 norm of each point in data cloud to origin
	normData = (float*)malloc(sizeof(float)*Nd);
	for(i = 0; i < Nd; i++)
	{
		normData[i] = sqrt(pData[i].x*pData[i].x + pData[i].y*pData[i].y + pData[i].z*pData[i].z);
	}

	maxRotDis = new float*[MAXROTLEVEL];
	for(i = 0; i < MAXROTLEVEL; i++)
	{
		maxRotDis[i] = (float*)malloc(sizeof(float*)*Nd);

		sigma = initNodeRot.w/pow(2.0,i)/2.0; // Half-side length of each level of rotation subcube
		maxAngle = SQRT3*sigma;

		if(maxAngle > PI)
			maxAngle = PI;
		for(j = 0; j < Nd; j++)
			maxRotDis[i][j] = 2*sin(maxAngle/2)*normData[j];
	}

	// Temporary Variables
	minDis = (float*)malloc(sizeof(float)*Nd);
	pDataTemp = (POINT3D *)malloc(sizeof(POINT3D)*Nd);
	pDataTempICP = (POINT3D *)malloc(sizeof(POINT3D)*Nd);

	// ICP Initialisation
	// Copy model and data point clouds to variables for ICP
	M_icp = (float*)calloc(3*Nm,sizeof(float));
	D_icp = (float*)calloc(3*Nd,sizeof(float));
	for(i = 0, j = 0; i < Nm; i++)
	{
		M_icp[j++] = pModel[i].x;
		M_icp[j++] = pModel[i].y;
		M_icp[j++] = pModel[i].z;
	}
	for(i = 0, j = 0; i < Nd; i++)
	{
		D_icp[j++] = pData[i].x;
		D_icp[j++] = pData[i].y;
		D_icp[j++] = pData[i].z;
	}

	// Build ICP kdtree with model dataset
	icp3d.Build(M_icp,Nm);
	icp3d.err_diff_def = MSEThresh/10000;
	icp3d.trim_fraction = trimFraction;
	icp3d.do_trim = doTrim;

	// Initialise so-far-best rotation and translation nodes
	optNodeRot = initNodeRot;
	optNodeTrans = initNodeTrans;
	// Initialise so-far-best rotation and translation matrices
	optR = Matrix::eye(3);
	optT = Matrix::ones(3,1)*0;

	// For untrimmed ICP, use all points, otherwise only use inlierNum points
	if(doTrim)
	{
		// Calculate number of inlier points
		inlierNum = (int)(Nd * (1 - trimFraction));
	}
	else
	{
		inlierNum = Nd;
	}

	//////////////////////////////////////////////////////
	//Assignation of neighbors to discard points with the least numbers, and weights ponderation
	//assignNeighborsPrio();
	weights = (float*)malloc(sizeof(float)*Nd);
	for (size_t i = 0; i < Nd; i++)
	{
		weights[i] = 1;
	}
	if(ponderation == 1) neighborsWeights();
	//neighborsDensity();
	//////////////////////////////////////////////////////

	SSEThresh = MSEThresh * inlierNum;
}

void GoICP::Clear()
{
	delete(pDataTemp);
	delete(pDataTempICP);
	delete(normData);
	delete(minDis);
	for(int i = 0; i < MAXROTLEVEL; i++)
	{
		delete(maxRotDis[i]);
	}
	delete(maxRotDis);
	delete(M_icp);
	delete(D_icp);
	//delete(compBNB);
}

// Inner Branch-and-Bound, iterating over the translation space
float GoICP::InnerBnB(float* maxRotDisL, TRANSNODE* nodeTransOut)
{
    int i, j, c;
	float transX, transY, transZ;
	float lb, ub, optErrorT;
	float dis, maxTransDis;
	TRANSNODE nodeTrans, nodeTransParent;
	priority_queue<TRANSNODE> queueTrans;

	// Set optimal translation error to overall so-far optimal error
	// Investigating translation nodes that are sub-optimal overall is redundant
	optErrorT = optError;

	// Push top-level translation node into the priority queue
	queueTrans.push(initNodeTrans);

	//////////////////////////////////////////////////////
	//Keep in memory the incompabilities and c-FPFH difference for the already tested translations to save time
	vector<TRANSCOMPATIBILITIES> storedCompatibilities;
	vector<TRANSFPFH> storedFPFH;
	//////////////////////////////////////////////////////
	
	while(1)
	{
		if(queueTrans.empty()) {
			break;
		}

		nodeTransParent = queueTrans.top();
		queueTrans.pop();

        if(optErrorT-nodeTransParent.lb < SSEThresh)
		{
			break;
		}

		nodeTrans.w = nodeTransParent.w/2;
		maxTransDis = SQRT3/2.0*nodeTrans.w;

		for(j = 0; j < 8; j++)
		{
			nodeTrans.x = nodeTransParent.x + (j&1)*nodeTrans.w ;
			nodeTrans.y = nodeTransParent.y + (j>>1&1)*nodeTrans.w ;
			nodeTrans.z = nodeTransParent.z + (j>>2&1)*nodeTrans.w ;

			transX = nodeTrans.x + nodeTrans.w/2;
			transY = nodeTrans.y + nodeTrans.w/2;
			transZ = nodeTrans.z + nodeTrans.w/2;

            int comp = 0;
            int notComp = 0;

			/*int countScattering = 0;
			int countPlanarity = 0;
			int countArete = 0;*/
			
			// For each data point, calculate the distance to it's closest point in the model cloud
			for(i = 0; i < Nd; i++)
			{
				//////////////////////////////////////////////////////
				//METHOD1: points always compatible and incompatible according to 8 corners
                /*int itmp, ctmp = 0;
                for(c = 0; c < 8; c++) {
                    if(checkCompatibility(i, pDataTemp[i].x + corners[c].x, pDataTemp[i].y + corners[c].y, pDataTemp[i].z + corners[c].z)) {
                        ctmp++;
                    }
                    else {
                        itmp++;
                    }
                }
                if(ctmp == 8) {
                    comp++;
                }
                if(itmp == 8) {
                    notComp++;
                }*/
				//////////////////////////////////////////////////////
				
				// Find distance between transformed point and closest point in model set ||R_r0 * x + t0 - y||
				// pDataTemp is the data points rotated by R0
                //////////////////////////////////////////////////////
				//Add weights, density
				int x, y, z = 0;
                minDis[i] = weights[i] * dt.Distance(pDataTemp[i].x + transX, pDataTemp[i].y + transY, pDataTemp[i].z + transZ, x, y, z);
				//minDis[i] = minDis[i] * computeDensityDifference(false, i, pDataTemp[i].x + transX, pDataTemp[i].y + transY, pDataTemp[i].z + transZ);
				//////////////////////////////////////////////////////

				// Subtract the rotation uncertainty radius if calculating the rotation lower bound
				// maxRotDisL == NULL when calculating the rotation upper bound
				if(maxRotDisL)
					minDis[i] -= maxRotDisL[i];

				if(minDis[i] < 0)
				{
					minDis[i] = 0;
				}
			}

			if(doTrim)
			{
				// Sort by distance
				//qsort(minDis, Nd, sizeof(float), cmp);
				//myqsort(minDis, Nd, inlierNum);
				intro_select(minDis,0,Nd-1,inlierNum-1);
			}

			// For each data point, find the incremental upper and lower bounds
			ub = 0;
			for(i = 0; i < inlierNum; i++)
            {
				//////////////////////////////////////////////////////
				//Norm distinction
				if(norm == 2) ub += minDis[i] * minDis[i];
            	if(norm == 1) ub += minDis[i];
				//////////////////////////////////////////////////////
			}

			lb = 0;
			for(i = 0; i < inlierNum; i++)
			{
				// Subtract the translation uncertainty radius
				dis = minDis[i] - maxTransDis;
				if(dis > 0) {
					//////////////////////////////////////////////////////
					//Norm distinction
					if(norm == 2) lb += dis * dis;
            		if(norm == 1) lb += dis;
					//////////////////////////////////////////////////////
				}
			}

			//////////////////////////////////////////////////////
			//METHOD2: number of incompatibilities for translation by center of cube to upper and lower bounds
			/*int incomp = countCompatibilities(false);
            ub += regularization * ((Nd-comp)*(Nd-comp));
            ub += regularization * ((incomp)*(incomp));
            lb += regularization * ((notComp)*(notComp));
            lb += regularization * ((incomp)*(incomp));

			int neighbors = compareNeighbors(false, transX, transY, transZ);
			ub += regularizationNeighbors * (neighbors*neighbors);
            lb += regularizationNeighbors * (neighbors*neighbors);*/

			//METHOD3: highest number of incompatibilities to upper bound and lowest to lower bound according to 8 corners
			//Same with direct neighbors, density and c-FPFH
			int minNeighbors, maxNeighbors, minIncomp, maxIncomp = 0;
			//float minDens, maxDens = 0;
			float minFPFH, maxFPFH = 0;
			//Check for the 8 corners of translation cube. Start with the first corner to have base values of max and min.
			//At each corner, we check if the translation has already be made, if it's the case we take the stored results, else we compute and we store
			if(regularization > 0 || regularizationNeighbors > 0 || ( regularizationFPFH > 0 && cfpfh != 0 ) ) {	
				float xIndex = nodeTrans.x;
				float yIndex = nodeTrans.y;
				float zIndex = nodeTrans.z;
				if(regularization > 0) {
					int incompFirst = 0;
					auto it = find_if(storedCompatibilities.begin(), storedCompatibilities.end(), [&xIndex, &yIndex, &zIndex](const TRANSCOMPATIBILITIES& obj) { if (obj.x == xIndex && obj.y == yIndex && obj.z == zIndex) return true; return false; });
    				if(it != storedCompatibilities.end()) {
 						incompFirst = it->c;
					}
					else {
						TRANSCOMPATIBILITIES t;
						t.x = xIndex;
						t.y = yIndex;
						t.z = zIndex;
						t.c = checkCompatibilities(xIndex, yIndex, zIndex);
    				    storedCompatibilities.push_back(t);
						incompFirst = t.c;
					}
					maxIncomp = incompFirst;
					minIncomp = incompFirst;

					/*float densFirst = sumDensities(xIndex, yIndex, zIndex);
					maxDens = densFirst;
					minDens = densFirst;*/			
				}
				if(regularizationNeighbors > 0) {
					int neighborsFirst = compareNeighbors(false, xIndex, yIndex, zIndex);
					maxNeighbors = neighborsFirst;
					minNeighbors = neighborsFirst;
				}
				if(regularizationFPFH > 0) {
					int FPFHFirst = 0;
					auto it = find_if(storedFPFH.begin(), storedFPFH.end(), [&xIndex, &yIndex, &zIndex](const TRANSFPFH& obj) { if (obj.x == xIndex && obj.y == yIndex && obj.z == zIndex) return true; return false; });
    				if(it != storedFPFH.end()) {
 						FPFHFirst = it->fpfh;
					}
					else {
						TRANSFPFH t;
						t.x = xIndex;
						t.y = yIndex;
						t.z = zIndex;
						t.fpfh = sumFPFH(xIndex, yIndex, zIndex);
    				    storedFPFH.push_back(t);
						FPFHFirst = t.fpfh;
					}
					maxFPFH = FPFHFirst;
					minFPFH = FPFHFirst;
				}
				for(c = 1; c < 8; c++) {
					xIndex = nodeTrans.x + (c&1)*nodeTrans.w;
					yIndex = nodeTrans.y + (c>>1&1)*nodeTrans.w;
					zIndex = nodeTrans.z + (c>>2&1)*nodeTrans.w;
					if(regularizationNeighbors > 0) {
						int tempNeighbors = compareNeighbors(false, xIndex, yIndex, zIndex);
						if(tempNeighbors > maxNeighbors) maxNeighbors = tempNeighbors;
						if(tempNeighbors < minNeighbors) minNeighbors = tempNeighbors;
					}
					if(regularizationFPFH > 0) {
						int tempFPFH = 0;
						auto it2 = find_if(storedFPFH.begin(), storedFPFH.end(), [&xIndex, &yIndex, &zIndex](const TRANSFPFH& obj) { if (obj.x == xIndex && obj.y == yIndex && obj.z == zIndex) return true; return false; });
    					if(it2 != storedFPFH.end()) {
							tempFPFH = it2->fpfh;
						}
						else {
							TRANSFPFH t;
							t.x = xIndex;
							t.y = yIndex;
							t.z = zIndex;
							t.fpfh = sumFPFH(xIndex, yIndex, zIndex);
    					    storedFPFH.push_back(t);
							tempFPFH = t.fpfh;
    					}
						if(tempFPFH > maxFPFH) maxFPFH = tempFPFH;
						if(tempFPFH < minFPFH) minFPFH = tempFPFH;
					}
					if(regularization > 0) {
						int tempComp = 0;
						auto it2 = find_if(storedCompatibilities.begin(), storedCompatibilities.end(), [&xIndex, &yIndex, &zIndex](const TRANSCOMPATIBILITIES& obj) { if (obj.x == xIndex && obj.y == yIndex && obj.z == zIndex) return true; return false; });
    					if(it2 != storedCompatibilities.end()) {
							tempComp = it2->c;
						}
						else {
							TRANSCOMPATIBILITIES t;
							t.x = xIndex;
							t.y = yIndex;
							t.z = zIndex;
							t.c = checkCompatibilities(xIndex, yIndex, zIndex);
    					    storedCompatibilities.push_back(t);
							tempComp = t.c;
    					}
						if(tempComp > maxIncomp) maxIncomp = tempComp;
						if(tempComp < minIncomp) minIncomp = tempComp;
						
						/*float tempDens = sumDensities(xIndex, yIndex, zIndex);
						if(tempDens > maxDens) maxDens = tempDens;
						if(tempDens < minDens) minDens = tempDens;*/			
					}	
				}
				//Add the values with regularization to the bounds
				if(regularization > 0) {
					ub += regularization * (maxIncomp*maxIncomp);
					lb += regularization * (minIncomp*minIncomp);
					/*ub += regularization * (maxDens*maxDens);
					lb += regularization * (minDens*minDens);*/
				}
				if(regularizationNeighbors > 0) {
					ub += regularizationNeighbors * (maxNeighbors*maxNeighbors);
					lb += regularizationNeighbors * (minNeighbors*minNeighbors);
				}
				if(regularizationFPFH > 0) {
					ub += regularizationFPFH * (maxFPFH*maxFPFH);
					lb += regularizationFPFH * (minFPFH*minFPFH);
				}
			}	
			//////////////////////////////////////////////////////		

			// If upper bound is better than best, update optErrorT and optTransOut (optimal translation node)
			if(ub < optErrorT)
			{          
				//////////////////////////////////////////////////////
				//Decomposition of error
				incompErrorTemp = regularization * (maxIncomp*maxIncomp);
				FPFHErrorTemp = regularizationFPFH * (maxFPFH*maxFPFH);
				geomErrorTemp = abs(ub - incompErrorTemp - FPFHErrorTemp);
				//////////////////////////////////////////////////////
				
				optErrorT = ub;
				if(nodeTransOut)
					*nodeTransOut = nodeTrans;
			}
			// Remove subcube from queue if lb is bigger than optErrorT
			if(lb >= optErrorT)
			{
				//discard
				continue;
			}
			nodeTrans.ub = ub;
			nodeTrans.lb = lb;
			queueTrans.push(nodeTrans);
		}
	}
	return optErrorT;
}

// Outer Branch-and-Bound, iterating over the rotation space
float GoICP::OuterBnB()
{
	int i, j;
	ROTNODE nodeRot, nodeRotParent;
	TRANSNODE nodeTrans;
	float v1, v2, v3, t, ct, ct2,st, st2;
	float tmp121, tmp122, tmp131, tmp132, tmp231, tmp232;
	float R11, R12, R13, R21, R22, R23, R31, R32, R33;
	float lb, ub, error, dis;
    clock_t clockBeginICP, clockBegin;
	priority_queue<ROTNODE> queueRot;

    clockBegin = clock();

	// Calculate Initial Error
	optError = 0;

	//////////////////////////////////////////////////////
	//Add weights
    int x, y, z = 0;
	for(i = 0; i < Nd; i++)
	{
        minDis[i] = weights[i] * dt.Distance(pData[i].x, pData[i].y, pData[i].z, x, y, z);
	}
	//////////////////////////////////////////////////////
	
	if(doTrim)
	{
		// Sort by distance
		//qsort(minDis, Nd, sizeof(float), cmp);
		//myqsort(minDis, Nd, inlierNum);
		intro_select(minDis,0,Nd-1,inlierNum-1);
	}

	//////////////////////////////////////////////////////
	//Distinction according to norm, add to initial error the different parameters with their regularization
	for(i = 0; i < inlierNum; i++)
	{
		if(norm == 2) optError += minDis[i] * minDis[i];
        if(norm == 1) optError += minDis[i];
    }
	if(regularization > 0) optError += regularization * (Nd*Nd);
	if(regularizationFPFH > 0) optError += regularizationFPFH * (100*8 * 100*8);
	if(regularizationNeighbors > 0) optError += regularizationNeighbors * (Nd*6*Nd*6);
	//////////////////////////////////////////////////////
	cout << "Error*: " << optError << " (Init)" << endl;

	Matrix R_icp = optR;
	Matrix t_icp = optT;

	// Run ICP from initial state
	clockBeginICP = clock();
	error = ICP(R_icp, t_icp);
    icp3d.optPoints = icp3d.points;
	if(error < optError)
	{
		optError = error;
		//////////////////////////////////////////////////////
		//ICP is the last operation to reduce error
		lastICP = true;
		geomError = geomErrorTemp;
		incompError = incompErrorTemp;
		FPFHError = FPFHErrorTemp;
		//////////////////////////////////////////////////////
		optR = R_icp;
		optT = t_icp;
		//////////////////////////////////////////////////////
		//Count incompatibilities of initial ICP
		optComp = countCompatibilities(false);
		//maxComp = Nd-optComp;
        cout << "Begin ICP Compatibilities: " << Nd-optComp << endl;
		//////////////////////////////////////////////////////
		cout << "Error*: " << error << " (ICP " << (double)(clock()-clockBeginICP)/CLOCKS_PER_SEC << "s)" << endl;
		cout << "Geometrical error: " << geomError << " / Incompatibilities error: " << incompError << " / c-FPFH error: " << FPFHError << endl;
		cout << "ICP-ONLY Rotation Matrix:" << endl;
		cout << R_icp << endl;
		cout << "ICP-ONLY Translation Vector:" << endl;
		cout << t_icp << endl;
		cout << endl;
	}

	// Push top-level rotation node into priority queue
	queueRot.push(initNodeRot);

	// Keep exploring rotation space until convergence is achieved
	long long count = 0;
	while(1)
	{
		if(queueRot.empty())
		{
		  	cout << "Rotation Queue Empty" << endl;
		  	cout << "End Compatibilities: " << Nd-optComp << endl;
		  	cout << "Error*: " << optError << ", LB: " << lb << endl;
			cout << "Geometrical error: " << geomError << " / Incompatibilities error: " << incompError << " / c-FPFH error: " << FPFHError << endl;
		  break;
		}

		// Access rotation cube with lowest lower bound...
		nodeRotParent = queueRot.top();
		// ...and remove it from the queue
		queueRot.pop();

        // Exit if the optError is less than or equal to the lower bound plus a small epsilon
        if((optError-nodeRotParent.lb) <= SSEThresh)
		{
            cout << "Threshold reached" << endl;
            cout << "End Compatibilities: " << Nd-optComp << endl;
			cout << "Error*: " << optError << ", LB: " << nodeRotParent.lb << ", epsilon: " << SSEThresh << endl;
			cout << "Geometrical error: " << geomError << " / Incompatibilities error: " << incompError << " / c-FPFH error: " << FPFHError << endl;
			break;
		}

        if(count>0 && count%300 == 0) {
            cout << "LB=" << nodeRotParent.lb << "  L=" << nodeRotParent.l << "  Time:" << (double)(clock()-clockBegin)/CLOCKS_PER_SEC << "s" << endl;
            /*cout << "Rotation Matrix:" << endl;
            cout << optR << endl;
            cout << "Translation Vector:" << endl;
            cout << optT << endl;*/
        }
		count ++;
		
		// Subdivide rotation cube into octant subcubes and calculate upper and lower bounds for each
		nodeRot.w = nodeRotParent.w/2;
		nodeRot.l = nodeRotParent.l+1;
		// For each subcube,
		for(j = 0; j < 8; j++)
		{
		  // Calculate the smallest rotation across each dimension
			nodeRot.a = nodeRotParent.a + (j&1)*nodeRot.w ;
			nodeRot.b = nodeRotParent.b + (j>>1&1)*nodeRot.w ;
			nodeRot.c = nodeRotParent.c + (j>>2&1)*nodeRot.w ;

            //cout << "ROTATION " << j << " - " << count << endl;
            //cout << nodeRot.a << " " << nodeRot.b << " " << nodeRot.c << " " << nodeRot.l << " " << nodeRot.w << " " << nodeRot.lb << " " << nodeRot.ub << endl;

			// Find the subcube centre
			v1 = nodeRot.a + nodeRot.w/2;
			v2 = nodeRot.b + nodeRot.w/2;
			v3 = nodeRot.c + nodeRot.w/2;

			// Skip subcube if it is completely outside the rotation PI-ball
			if(sqrt(v1*v1+v2*v2+v3*v3)-SQRT3*nodeRot.w/2 > PI)
			{
				continue;
			}

			// Convert angle-axis rotation into a rotation matrix
			t = sqrt(v1*v1 + v2*v2 + v3*v3);
			if(t > 0)
			{
				v1 /= t;
				v2 /= t;
				v3 /= t;

				ct = cos(t);
				ct2 = 1 - ct;
				st = sin(t);
				st2 = 1 - st;

				tmp121 = v1*v2*ct2; tmp122 = v3*st;
				tmp131 = v1*v3*ct2; tmp132 = v2*st;
				tmp231 = v2*v3*ct2; tmp232 = v1*st;

				R11 = ct + v1*v1*ct2;		R12 = tmp121 - tmp122;		R13 = tmp131 + tmp132;
				R21 = tmp121 + tmp122;		R22 = ct + v2*v2*ct2;		R23 = tmp231 - tmp232;
				R31 = tmp131 - tmp132;		R32 = tmp231 + tmp232;		R33 = ct + v3*v3*ct2;

				// Rotate data points by subcube rotation matrix
				for(i = 0; i < Nd; i++)
				{
					POINT3D& p = pData[i];
					pDataTemp[i].x = R11*p.x + R12*p.y + R13*p.z;
					pDataTemp[i].y = R21*p.x + R22*p.y + R23*p.z;
					pDataTemp[i].z = R31*p.x + R32*p.y + R33*p.z;
				}
			}
			// If t == 0, the rotation angle is 0 and no rotation is required
			else
			{
				memcpy(pDataTemp, pData, sizeof(POINT3D)*Nd);
			}

			// Upper Bound
			// Run Inner Branch-and-Bound to find rotation upper bound
			// Calculates the rotation upper bound by finding the translation upper bound for a given rotation,
			// assuming that the rotation is known (zero rotation uncertainty radius)
			ub = InnerBnB(NULL /*Rotation Uncertainty Radius*/, &nodeTrans);

            // If the upper bound is the best so far, run ICP
			if(ub < optError)
			{
				//////////////////////////////////////////////////////
				//BnB is the last operation to reduce error
				lastICP = false;
				//////////////////////////////////////////////////////

				// Update optimal error and rotation/translation nodes
				optError = ub;
				optNodeRot = nodeRot;
				optNodeTrans = nodeTrans;

				optR.val[0][0] = R11; optR.val[0][1] = R12; optR.val[0][2] = R13;
				optR.val[1][0] = R21; optR.val[1][1] = R22; optR.val[1][2] = R23;
				optR.val[2][0] = R31; optR.val[2][1] = R32; optR.val[2][2] = R33;
				optT.val[0][0] = optNodeTrans.x+optNodeTrans.w/2;
				optT.val[1][0] = optNodeTrans.y+optNodeTrans.w/2;
				optT.val[2][0] = optNodeTrans.z+optNodeTrans.w/2;

				//////////////////////////////////////////////////////
				updateCompatibilities();

				geomError = geomErrorTemp;
				incompError = incompErrorTemp;
				FPFHError = FPFHErrorTemp;
				//////////////////////////////////////////////////////

                cout << "BNB Compatibilities: " << Nd-optComp << endl;
                cout << "Error*: " << optError << " (BNB " << (double)(clock()-clockBegin)/CLOCKS_PER_SEC << "s)" << endl;
				cout << "Geometrical error: " << geomError << " / Incompatibilities error: " << incompError << " / c-FPFH error: " << FPFHError << endl;
                cout << "BNB Rotation Matrix:" << endl;
                cout << optR << endl;
                cout << "BNB Translation Vector:" << endl;
                cout << optT << endl;

				// Run ICP
				clockBeginICP = clock();
				R_icp = optR;
				t_icp = optT;
				error = ICP(R_icp, t_icp);
				//Our ICP implementation uses kdtree for closest distance computation which is slightly different from DT approximation, 
				//thus it's possible that ICP failed to decrease the DT error. This is no big deal as the difference should be very small.
				if(error < optError)
				{
					//////////////////////////////////////////////////////
					lastICP = true;
					//////////////////////////////////////////////////////

					optError = error;
					optR = R_icp;
					optT = t_icp;

					//////////////////////////////////////////////////////
                    icp3d.optPoints = icp3d.points;
					optComp = countCompatibilities(false);

					geomError = geomErrorTemp;
					incompError = incompErrorTemp;
					FPFHError = FPFHErrorTemp;
					//////////////////////////////////////////////////////

                    cout << "ICP Compatibilities: " << Nd-optComp << endl;
                    cout << "Error*: " << error << " (ICP " << (double)(clock() - clockBeginICP)/CLOCKS_PER_SEC << "s)" << endl;
					cout << "Geometrical error: " << geomError << " / Incompatibilities error: " << incompError << " / c-FPFH error: " << FPFHError << endl;
                    cout << "ICP Rotation Matrix:" << endl;
                    cout << R_icp << endl;
                    cout << "ICP Translation Vector:" << endl;
                    cout << t_icp << endl;
					cout << endl;
				}

				// Discard all rotation nodes with high lower bounds in the queue
				priority_queue<ROTNODE> queueRotNew;
				while(!queueRot.empty())
				{
					ROTNODE node = queueRot.top();
					queueRot.pop();
					if(node.lb < optError)
						queueRotNew.push(node);
					else
						break;
				}
				queueRot = queueRotNew;
			}

            // Lower Bound
			// Run Inner Branch-and-Bound to find rotation lower bound
			// Calculates the rotation lower bound by finding the translation upper bound for a given rotation,
			// assuming that the rotation is uncertain (a positive rotation uncertainty radius)
			// Pass an array of rotation uncertainties for every point in data cloud at this level
			lb = InnerBnB(maxRotDis[nodeRot.l], NULL /*Translation Node*/);
            // If the best error so far is less than the lower bound, remove the rotation subcube from the queue
			if(lb >= optError)
			{
				continue;
			}

			// Update node and put it in queue
			nodeRot.ub = ub;
			nodeRot.lb = lb;
			queueRot.push(nodeRot);
		}
    }

	return optError;
}

float GoICP::Register()
{
	Initialize();
	OuterBnB();
	Clear();

	return optError;
}

/**
 * Count incompatibilties from ICP
 */
int GoICP::countCompatibilities(bool init) {
    int countComp = 0;
    int countNotComp = 0;
    for(int i = 0; i < Nd; i++) {
        properties propertySource;
        properties propertyTarget;
		//On points transformed during ICP
        if(init) {
            propertySource = (properties)pData[icp3d.points[i].id_data].c;
            propertyTarget = (properties)pModel[icp3d.points[i].id_model].c;
        }
		//On points transformed by an ICP which reduced the error
        else {
            propertySource = (properties)pData[icp3d.optPoints[i].id_data].c;
            propertyTarget = (properties)pModel[icp3d.optPoints[i].id_model].c;
        }
        if(find(compatibilities[propertySource].begin(), compatibilities[propertySource].end(), propertyTarget) != compatibilities[propertySource].end()) {
            countComp++;
        }
        else {
            countNotComp++;
        }
    }
    return countNotComp;
}

/**
 * Count incompatibilties from BNB
 */
int GoICP::checkCompatibilities(float x, float y, float z) {
	int countNotComp = 0;
	for (int i = 0; i < Nd; i++)
	{
		if (!checkCompatibility(i, pDataTemp[i].x + x, pDataTemp[i].y + y, pDataTemp[i].z + z, false)) {
			countNotComp++;
		}
	}
	return countNotComp;
}

/**
 * Update the number of incompatibilities from BNB which reduces the error
 */
void GoICP::updateCompatibilities() {
    int countNotComp = 0;
    for (int i = 0; i < Nd; i++) {
        float x, y, z;
        x = optR.val[0][0]*pData[i].x + optR.val[0][1]*pData[i].y + optR.val[0][2]*pData[i].z + optT.val[0][0];
        y = optR.val[1][0]*pData[i].x + optR.val[1][1]*pData[i].y + optR.val[1][2]*pData[i].z + optT.val[1][0];
        z = optR.val[2][0]*pData[i].x + optR.val[2][1]*pData[i].y + optR.val[2][2]*pData[i].z + optT.val[2][0];
        if (!checkCompatibility(i, x, y, z, true)) {
            countNotComp++;
        }
    }
	//cout << endl;
	optComp = countNotComp;
}

/**
 * Assign the property color of each DT cell 
 */
void GoICP::assignCellColor() {
	for(int i = 0; i < dt.SIZE; i++) {
		for(int j = 0; j < dt.SIZE; j++) {
			for(int k = 0; k < dt.SIZE; k++) {
				//If DT cell contains at least one point, check property of the points in the cell
				if(dt.cellPoints[k][j][i].c == -1) {
					int prop = pModel[dt.cellPoints[k][j][i].points[0]].c;
					dt.cellPoints[k][j][i].c = prop;
           	 		for(int p = 1; p < dt.cellPoints[k][j][i].points.size(); p++) {
                		if(pModel[dt.cellPoints[k][j][i].points[p]].c != prop) {
                			dt.cellPoints[k][j][i].c = -1;
							break;
                		}
            		}
				}
			}
		}
	}
}

/**
 * Check if the transformed point is compatible with the property of the closest occupied cell
 */
bool GoICP::checkCompatibility(int i, float x, float y, float z, bool update) {
    bool comp = false;
    int rx = ROUND((x - dt.xMin) * dt.scale);
    int ry = ROUND((y - dt.yMin) * dt.scale);
    int rz = ROUND((z - dt.zMin) * dt.scale);
    if(rx < 0) rx = 0;
    if(rx>=dt.SIZE) rx = dt.SIZE - 1;
    if(ry < 0) ry = 0;
    if(ry>=dt.SIZE) ry = dt.SIZE - 1;
    if(rz < 0) rz = 0;
    if(rz>=dt.SIZE) rz = dt.SIZE - 1;
    int cx = dt.emptyCells[rz][ry][rx].cx;
    int cy = dt.emptyCells[rz][ry][rx].cy;
    int cz = dt.emptyCells[rz][ry][rx].cz;
    properties propertySource = (properties)pData[i].c;
    properties propertyTarget = (properties)dt.cellPoints[cz][cy][cx].c;
	if(update) {
		/*string s = to_string(i) + "-" + to_string(dt.cellPoints[cz][cy][cx].points[0]);
		if(dt.cellPoints[cz][cy][cx].points.size() > 1) {
			for(int j = 1; j < dt.cellPoints[cz][cy][cx].points.size(); j++) {
				s = s + "/" + to_string(dt.cellPoints[cz][cy][cx].points[j]);
			}
		}*/
		//compBNB[i] = nearestNeighbor(x, y, z, cx, cy, cz);
		//cout << s << " ---> Nearest : " << compBNB[i] << endl;
	}
	comp = checkProperty(propertySource, propertyTarget, cx, cy, cz, x, y, z);
	
    /*if(propertyTarget == -1) {
        for(int p : dt.cellPoints[cz][cy][cx].points) {
            if(pModel[p].c == propertySource) {
                comp = true;
                break;
            }
        }
    }
    else if(find(compatibilities[propertySource].begin(), compatibilities[propertySource].end(), propertyTarget) != compatibilities[propertySource].end()) {
        comp = true;
    }*/

	/*properties propertyTargetLeft, propertyTargetRight, propertyTargetUp, propertyTargetDown, propertyTargetFront, propertyTargetBack;
	if(comp == false && cx > 0 && dt.cellPoints[cz][cy][cx-1].x > -1) {
		propertyTargetLeft = (properties)dt.cellPoints[cz][cy][cx-1].c;
		comp = checkProperty(propertySource, propertyTargetLeft, cx-1, cy, cz);	
	}
	if(comp == false && cx < dt.SIZE && dt.cellPoints[cz][cy][cx+1].x > -1) {
		propertyTargetRight = (properties)dt.cellPoints[cz][cy][cx+1].c;
		comp = checkProperty(propertySource, propertyTargetRight, cx+1, cy, cz);
	}
	if(comp == false && cy > 0 && dt.cellPoints[cz][cy-1][cx].x > -1) {
		propertyTargetDown = (properties)dt.cellPoints[cz][cy-1][cx].c;
		comp = checkProperty(propertySource, propertyTargetDown, cx, cy-1, cz);
	}
	if(comp == false && cy < dt.SIZE && dt.cellPoints[cz][cy+1][cx].x > -1) {
		propertyTargetUp = (properties)dt.cellPoints[cz][cy+1][cx].c;
		comp = checkProperty(propertySource, propertyTargetUp, cx, cy+1, cz);
	}
	if(comp == false && cz > 0 && dt.cellPoints[cz-1][cy][cx].x > -1) {
		propertyTargetFront = (properties)dt.cellPoints[cz-1][cy][cx].c;
		comp = checkProperty(propertySource, propertyTargetFront, cx, cy, cz-1);
	}
	if(comp == false && cz < dt.SIZE && dt.cellPoints[cz+1][cy][cx].x > -1) {
		propertyTargetBack = (properties)dt.cellPoints[cz+1][cy][cx].c;
		comp = checkProperty(propertySource, propertyTargetBack, cx, cy, cz+1);
	}*/
    
    return comp;
}

void GoICP::closestEmptyCell(int x, int y, int z) {
    cout << dt.emptyCells[z][y][x].cz << "-" << dt.emptyCells[z][y][x].cy << "-" << dt.emptyCells[z][y][x].cx;
}

/**
 * Sort distances of the matching points after a transformation
 */
void GoICP::sortDistances() {
	sort(minDis, minDis + Nd);
	//inlierNum = 0.75 * Nd;

	float threshold = minDis[(int)(inlierNum*0.9)];

	for(int i = inlierNum-1; i >= 0; i--) {
		if(minDis[i] == threshold) {
			inlierNum = i+1;
			break;
		}
	}

}

/**
 * Check if the two properties are compatible
 */
bool GoICP::checkProperty(properties source, properties target, int cx, int cy, int cz, double x, double y, double z) {
	bool comp = false;
	if(target == -1) {
        for(int p : dt.cellPoints[cz][cy][cx].points) {
            if(pModel[p].c == source) {
                comp = true;
                break;
            }
        }
		/*double minD = 100;
		int ind = 0;
		for(int p : dt.cellPoints[cz][cy][cx].points) {
			double distance = sqrt(pow((x - pModel[p].x),2) + pow((y - pModel[p].y),2) + pow((z - pModel[p].z),2));
			if (distance < minD) { 
				ind = p;
				minD = distance;
			}
		}
		target = (properties)pModel[ind].c;*/
    }
    else if(find(compatibilities[source].begin(), compatibilities[source].end(), target) != compatibilities[source].end()) {
        comp = true;
    }
	return comp;
}

/**
 * Check if 2 points are neighbors
 */
bool GoICP::isNeighbor(float radius, POINT3D p, POINT3D p2) {
    double distance = sqrt(pow((p2.x - p.x),2) + pow((p2.y - p.y),2) + pow((p2.z - p.z),2));
    if( distance < sqrt(radius) ) {
        return true;
    }
    return false;
}

////////// Planarity - Scattering method //////////
double GoICP::calculateMean(vector<POINT3D> points, char d){
    double mean = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
        if(d == 'x') {
            mean += points[i].x;
        }
        else if(d == 'y') {
            mean += points[i].y;
        }
        else if(d == 'z') {
            mean += points[i].z;
        }
    }
    mean /= points.size();
    return mean;
}

void GoICP::centralizePoints(vector<POINT3D> &points, double &xMean, double &yMean, double &zMean) {
    xMean = calculateMean(points, 'x');
    yMean = calculateMean(points, 'y');
    zMean = calculateMean(points, 'z');
    for (size_t i = 0; i < points.size(); i++)
    {
        points[i].x -= xMean;
        points[i].y -= yMean;
        points[i].z -= zMean;
    }
}

void GoICP::calculateCovarianceMatrix(vector<POINT3D> points, double matrix[9]) {
    double xMean = calculateMean(points, 'x');
    double yMean = calculateMean(points, 'y');
    double zMean = calculateMean(points, 'z');
    double totalXX = 0;
    double totalYY = 0;
    double totalZZ = 0;
    double totalXY = 0;
    double totalXZ = 0;
    double totalYZ = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
        totalXX += (points[i].x - xMean) * (points[i].x - xMean);
        totalYY += (points[i].y - yMean) * (points[i].y - yMean);
        totalZZ += (points[i].z - zMean) * (points[i].z - zMean);
        totalXY += (points[i].x - xMean) * (points[i].y - yMean);
        totalXZ += (points[i].x - xMean) * (points[i].z - zMean);
        totalYZ += (points[i].y - yMean) * (points[i].z - zMean);
    }
    totalXX /= points.size()-1;
    totalYY /= points.size()-1;
    totalZZ /= points.size()-1;
    totalXY /= points.size()-1;
    totalXZ /= points.size()-1;
    totalYZ /= points.size()-1;
    matrix[0] = totalXX;
    matrix[1] = totalXY;
    matrix[2] = totalXZ;
    matrix[3] = totalXY;
    matrix[4] = totalYY;
    matrix[5] = totalYZ;
    matrix[6] = totalXZ;
    matrix[7] = totalYZ;
    matrix[8] = totalZZ;
}

void GoICP::calculateEigen(double matrix[9], vector<double> &eigenvalues) {
    /*Eigen::Matrix<double, 3, 3> A;
    A << matrix[0], matrix[1], matrix[2], matrix[3], matrix[4], matrix[5], matrix[6], matrix[7], matrix[8];
    Eigen::EigenSolver<Eigen::Matrix<double, 3, 3> > s(A);
    eigenvalues.push_back(real(s.eigenvalues()(0)));
    eigenvalues.push_back(real(s.eigenvalues()(1)));
    eigenvalues.push_back(real(s.eigenvalues()(2)));
    for(int i = eigenvalues.size()-1; i >= 0 ; i--) {
        for (size_t j = 1; j <= i; j++)
        {
            if(eigenvalues[j-1] < eigenvalues[j]) {
                double tmp = eigenvalues[j-1];
                eigenvalues[j-1] = eigenvalues[j];
                eigenvalues[j] = tmp;
            }
        } 
    }*/
}

double GoICP::computePlanarity(double lambda1, double lambda2, double lambda3) {
    return (lambda2-lambda3)/lambda1;
}

double GoICP::computeScattering(double lambda1, double lambda3) {
    return lambda3/lambda1;
}

////////// Direct neighbors method //////////
int GoICP::nearestNeighbor(float x, float y, float z, int cx, int cy, int cz) {
	double minD = 100;
	int ind = 0;
	for(int p : dt.cellPoints[cz][cy][cx].points) {
		double distance = sqrt(pow((x - pModel[p].x),2) + pow((y - pModel[p].y),2) + pow((z - pModel[p].z),2));
		if (distance < minD) { 
			ind = p;
			minD = distance;
		}
	}
	return ind;
}

void GoICP::assignNeighbors() {
	int maxN = 0;
	int ind = 0;
	for (size_t i = 0; i < Nd; i++)
	{
		int count = 0;
		for (size_t j = 0; j < Nd; j++)
        { 
            if(j == i) {
                continue;
            }
            if(isNeighbor(0.050, pData[i], pData[j])) {
                count++;
            }
        }
		if (count>maxN) {
			maxN = count;
			ind = i;
		}
		pData[i].neighbors = count;
	}
	for (size_t i = 0; i < Nm; i++)
	{
		int count = 0;
		for (size_t j = 0; j < Nm; j++)
        { 
            if(j == i) {
                continue;
            }
            if(isNeighbor(0.050, pModel[i], pModel[j])) {
                count++;
            }
        }
		pModel[i].neighbors = count;
	}
}

int GoICP::compareNeighbors(bool icp, float x, float y, float z) {
	int sum = 0;

	if(icp) {
		for(int i = 0; i < Nd; i++) {
			int neighborsSource = pData[icp3d.points[i].id_data].neighbors;
			int neighborsTarget = pModel[icp3d.points[i].id_model].neighbors;
			int diff = abs(neighborsSource - neighborsTarget);
			sum += diff;
		}
	}
	else {
		for(int i = 0; i < Nd; i++) {
			double ax = pDataTemp[i].x + x;
			double ay = pDataTemp[i].y + y;
			double az = pDataTemp[i].z + z;
			int rx = ROUND((ax - dt.xMin) * dt.scale);
    		int ry = ROUND((ay - dt.yMin) * dt.scale);
    		int rz = ROUND((az - dt.zMin) * dt.scale);
    		if(rx < 0) rx = 0;
    		if(rx>=dt.SIZE) rx = dt.SIZE - 1;
    		if(ry < 0) ry = 0;
    		if(ry>=dt.SIZE) ry = dt.SIZE - 1;
    		if(rz < 0) rz = 0;
    		if(rz>=dt.SIZE) rz = dt.SIZE - 1;
    		int cx = dt.emptyCells[rz][ry][rx].cx;
    		int cy = dt.emptyCells[rz][ry][rx].cy;
    		int cz = dt.emptyCells[rz][ry][rx].cz;
			//compBNB[i] = nearestNeighbor(ax, ay, az, cx, cy, cz);

			int neighborsSource = pData[i].neighbors;
			int neighborsTarget = pModel[nearestNeighbor(ax, ay, az, cx, cy, cz)].neighbors;
			int diff = abs(neighborsSource - neighborsTarget);
			sum += diff;
		}
	}

	return sum;
}

int GoICP::compareNeighborsV2(bool icp, float x, float y, float z) {
	int sum = 0;

	if(icp) {
		for(int i = 0; i < Nd; i++) {
			int neighborsSource = pData[icp3d.points[i].id_data].neighbors;
			int neighborsTarget = pModel[icp3d.points[i].id_model].neighbors;
			int diff = abs(neighborsSource - neighborsTarget);
			if(diff > 3) sum += diff;
		}
	}
	else {
		for(int i = 0; i < Nd; i++) {
			double ax = pDataTemp[i].x + x;
			double ay = pDataTemp[i].y + y;
			double az = pDataTemp[i].z + z;
			int rx = ROUND((ax - dt.xMin) * dt.scale);
    		int ry = ROUND((ay - dt.yMin) * dt.scale);
    		int rz = ROUND((az - dt.zMin) * dt.scale);
    		if(rx < 0) rx = 0;
    		if(rx>=dt.SIZE) rx = dt.SIZE - 1;
    		if(ry < 0) ry = 0;
    		if(ry>=dt.SIZE) ry = dt.SIZE - 1;
    		if(rz < 0) rz = 0;
    		if(rz>=dt.SIZE) rz = dt.SIZE - 1;
    		int cx = dt.emptyCells[rz][ry][rx].cx;
    		int cy = dt.emptyCells[rz][ry][rx].cy;
    		int cz = dt.emptyCells[rz][ry][rx].cz;
			//compBNB[i] = nearestNeighbor(ax, ay, az, cx, cy, cz);

			int neighborsSource = pData[i].neighbors;
			int neighborsTarget = pModel[nearestNeighbor(ax, ay, az, cx, cy, cz)].neighbors;
			int diff = abs(neighborsSource - neighborsTarget);
			if(diff > 3) sum += diff;
		}
	}

	return sum;
}

int GoICP::compareNeighborsV3(bool icp, float x, float y, float z) {
	int sum = 0;

	if(icp) {
		for(int i = 0; i < Nd; i++) {
			int neighborsSource = pData[icp3d.points[i].id_data].neighbors;
			int neighborsTarget = pModel[icp3d.points[i].id_model].neighbors;
			if(neighborsSource == 0 || neighborsSource == 1 || neighborsSource == 2) {
				if(neighborsTarget == 3 || neighborsTarget == 4) {
					sum += 1;
				}
				else if(neighborsTarget == 5 || neighborsTarget == 6) {
					sum += 2;
				}
			}
			else if(neighborsSource == 3 || neighborsSource == 4) {
				if(neighborsTarget != 3 && neighborsTarget != 4) {
					sum += 1;
				}
			}
			else if(neighborsSource == 5 || neighborsSource == 6) {
				if(neighborsTarget == 0 || neighborsTarget == 1 || neighborsTarget == 2) {
					sum += 2;
				}
				else if(neighborsTarget == 3 || neighborsTarget == 4) {
					sum += 1;
				}
			}
		}
	}
	else {
		for(int i = 0; i < Nd; i++) {
			double ax = pDataTemp[i].x + x;
			double ay = pDataTemp[i].y + y;
			double az = pDataTemp[i].z + z;
			int rx = ROUND((ax - dt.xMin) * dt.scale);
    		int ry = ROUND((ay - dt.yMin) * dt.scale);
    		int rz = ROUND((az - dt.zMin) * dt.scale);
    		if(rx < 0) rx = 0;
    		if(rx>=dt.SIZE) rx = dt.SIZE - 1;
    		if(ry < 0) ry = 0;
    		if(ry>=dt.SIZE) ry = dt.SIZE - 1;
    		if(rz < 0) rz = 0;
    		if(rz>=dt.SIZE) rz = dt.SIZE - 1;
    		int cx = dt.emptyCells[rz][ry][rx].cx;
    		int cy = dt.emptyCells[rz][ry][rx].cy;
    		int cz = dt.emptyCells[rz][ry][rx].cz;
			//compBNB[i] = nearestNeighbor(ax, ay, az, cx, cy, cz);

			int neighborsSource = pData[i].neighbors;
			int neighborsTarget = pModel[nearestNeighbor(ax, ay, az, cx, cy, cz)].neighbors;
			if(neighborsSource == 0 || neighborsSource == 1 || neighborsSource == 2) {
				if(neighborsTarget == 3 || neighborsTarget == 4) {
					sum += 1;
				}
				else if(neighborsTarget == 5 || neighborsTarget == 6) {
					sum += 2;
				}
			}
			else if(neighborsSource == 3 || neighborsSource == 4) {
				if(neighborsTarget != 3 && neighborsTarget != 4) {
					sum += 1;
				}
			}
			else if(neighborsSource == 5 || neighborsSource == 6) {
				if(neighborsTarget == 0 || neighborsTarget == 1 || neighborsTarget == 2) {
					sum += 2;
				}
				else if(neighborsTarget == 3 || neighborsTarget == 4) {
					sum += 1;
				}
			}
		}
	}

	return sum;
}


/**
 * Check all source points neighbors and then discard the points with the least neighbors
 */
void GoICP::assignNeighborsPrio() {
	int maxN = 0;
	//int ind = 0;
	float distance = 0.035;
	while(maxN < 19) {
		for (size_t i = 0; i < Nd; i++)
		{
			int count = 0;
			for (size_t j = 0; j < Nd; j++)
    	    { 
    	        if(j == i) {
    	            continue;
    	        }
    	        if(isNeighbor(distance, pData[i], pData[j])) {
    	            count++;
    	        }
    	    }
			if (count>maxN) {
				maxN = count;
				//ind = i;
			}
			pData[i].neighbors = count;
		}
		distance = distance + 0.001;
	}

	sort(pData, pData + Nd, [](const POINT3D &a, const POINT3D &b){return a.neighbors > b.neighbors;});

	float threshold = pData[(int)(Nd*0.75)].neighbors;
	for(int i = Nd-1; i >= 0; i--) {
		if(pData[i].neighbors >= threshold) {
			Nd = i+1;
			break;
		}
	}
	inlierNum = Nd;
}

/**
 * Computes source points weights according to their number of neighbours
 */
void GoICP::neighborsWeights() {
	int maxN = 0;
	int minN = 100;
	//Arbitrary distance determined for minimum 3 nearest neighbors radius
	float distance = 0.035;
	//While points have least of 19 points, that means we are at less of 3 nearest neighbors radius
	while(maxN < 19) {
		for (size_t i = 0; i < Nd; i++)
		{
			int count = 0;
			//int countProperty = 0;
			for (size_t j = 0; j < Nd; j++)
    	    { 
    	        if(j == i) {
    	            continue;
    	        }
    	        if(isNeighbor(distance, pData[i], pData[j])) {
    	            count++;
					/*if(pData[i].c == pData[j].c) {
						countProperty += 2;
					}*/
    	        }
    	    }
			//pData[i].neighbors = count + countProperty;
			pData[i].neighbors = count;
			/*if ( (count + countProperty) >maxN) {
				maxN = count + countProperty;
			}*/
			if(count > maxN) {
				maxN = count;
			}
			if(count < minN) {
				minN = count;
			}
		}
		distance = distance + 0.001;
	}
	if(minN == 0) minN = 1;
	for (size_t i = 0; i < Nd; i++)
	{
		if(pData[i].neighbors == 0) pData[i].neighbors = 1;
		//float f = (1 - ( (float)pData[i].neighbors / (float)maxN ) ) * 2;
		float f = ((float)minN / (float)pData[i].neighbors) * 2;
		weights[i] += f; 
	}
}

/**
 * Computes source and target points weights and property density according to their number of neighbours
 */
void GoICP::neighborsDensity() {
	int maxN = 0;
	int minN = 100;
	float distance = 0.035;
	while(maxN < 19) {
		for (size_t i = 0; i < Nd; i++)
		{
			int count = 0;
			int countProperty = 0;
			for (size_t j = 0; j < Nd; j++)
    	    { 
    	        if(j == i) {
    	            continue;
    	        }
    	        if(isNeighbor(distance, pData[i], pData[j])) {
    	            count++;
					if(pData[i].c == pData[j].c) {
						countProperty++;
					}
    	        }
    	    }
			pData[i].neighbors = count;
			pData[i].density = (float)countProperty / (float)count;
			if(count > maxN) {
				maxN = count;
			}
			if(count < minN) {
				minN = count;
			}
		}
		distance = distance + 0.001;
	}
	if(minN == 0) minN = 1;
	for (size_t i = 0; i < Nd; i++)
	{
		if(pData[i].neighbors == 0) pData[i].neighbors = 1;
		//float f = (1 - ( (float)pData[i].neighbors / (float)maxN ) ) * 2;
		float f = ((float)minN / (float)pData[i].neighbors) * 2;
		weights[i] += f; 
	}
	maxN = 0;
	distance = 0.035;
	while(maxN < 19) {
		for (size_t i = 0; i < Nm; i++)
		{
			int count = 0;
			int countProperty = 0;
			for (size_t j = 0; j < Nm; j++)
    	    { 
    	        if(j == i) {
    	            continue;
    	        }
    	        if(isNeighbor(distance, pModel[i], pModel[j])) {
    	            count++;
					if(pModel[i].c == pModel[j].c) {
						countProperty++;
					}
    	        }
    	    }
			pModel[i].neighbors = count;
			pModel[i].density = (float)countProperty / (float)count;
			if(count > maxN) {
				maxN = count;
			}
		}
		distance = distance + 0.001;
	}
}

/**
 * Computes density difference
 */
float GoICP::computeDensityDifference(bool icp, int point, float x, float y, float z) {
	float d = 0;
	//If icp: compute difference between 2 corresponding points
	if(icp) {
		d =  ((abs(pData[icp3d.points[point].id_data].density - pModel[icp3d.points[point].id_model].density)) );
	}
	//If bnb: apply transformation on point to check closest occupied DT cell and compute smallest difference with points in the cell 
	else {
		int rx = ROUND((x - dt.xMin) * dt.scale);
    	int ry = ROUND((y - dt.yMin) * dt.scale);
    	int rz = ROUND((z - dt.zMin) * dt.scale);
    	if(rx < 0) rx = 0;
    	if(rx>=dt.SIZE) rx = dt.SIZE - 1;
    	if(ry < 0) ry = 0;
    	if(ry>=dt.SIZE) ry = dt.SIZE - 1;
    	if(rz < 0) rz = 0;
    	if(rz>=dt.SIZE) rz = dt.SIZE - 1;
    	int cx = dt.emptyCells[rz][ry][rx].cx;
    	int cy = dt.emptyCells[rz][ry][rx].cy;
    	int cz = dt.emptyCells[rz][ry][rx].cz;
		float minD = 100;
		for(int p : dt.cellPoints[cz][cy][cx].points) {
            if(  ((abs(pData[point].density - pModel[p].density) ) ) < minD) {
				minD =   ((abs(pData[point].density - pModel[p].density) ) );
			}
        }
		d = minD;
	}

	return d;
}

/**
 * Computes sum of density differences for BnB
 */
float GoICP::sumDensities(float x, float y, float z) {
	float sum = 0;
	for (int i = 0; i < Nd; i++)
	{
		sum += computeDensityDifference(false, i, pDataTemp[i].x + x, pDataTemp[i].y + y, pDataTemp[i].z + z);
	}
	return sum;
}

/**
 * Computes c-FPFH difference
 */
float GoICP::computeFPFHDifference(bool icp, int point, float x, float y, float z) {
	float d = 0;
	//if icp: compute difference between 2 corresponding points
	if(icp) {
		if(cfpfh == 1) {
			for(int i = 0; i < 41; i++) {
				d += abs(pData[icp3d.points[point].id_data].cfpfh[i] - pModel[icp3d.points[point].id_model].cfpfh[i]);
			}	
		}
		else if(cfpfh == 2) {
			for(int i = 0; i < 33; i++) {
				d += abs(pData[icp3d.points[point].id_data].cfpfh[i] - pModel[icp3d.points[point].id_model].cfpfh[i]);
			}	
		}
		else if(cfpfh == 3) {
			for(int i = 33; i < 41; i++) {
				d += abs(pData[icp3d.points[point].id_data].cfpfh[i] - pModel[icp3d.points[point].id_model].cfpfh[i]);
			}	
		}
	}
	//if BnB: apply transformation on point to check closest occupied DT cell and compute smallest difference with points in the cell 
	else {
		int rx = ROUND((x - dt.xMin) * dt.scale);
    	int ry = ROUND((y - dt.yMin) * dt.scale);
    	int rz = ROUND((z - dt.zMin) * dt.scale);
    	if(rx < 0) rx = 0;
    	if(rx>=dt.SIZE) rx = dt.SIZE - 1;
    	if(ry < 0) ry = 0;
    	if(ry>=dt.SIZE) ry = dt.SIZE - 1;
    	if(rz < 0) rz = 0;
    	if(rz>=dt.SIZE) rz = dt.SIZE - 1;
    	int cx = dt.emptyCells[rz][ry][rx].cx;
    	int cy = dt.emptyCells[rz][ry][rx].cy;
    	int cz = dt.emptyCells[rz][ry][rx].cz;
		float minD = 1000000000;
		for(int p : dt.cellPoints[cz][cy][cx].points) {
			float diff = 0;

			if(cfpfh == 1) {
				for(int i = 0; i < 41; i++) {
					diff += (abs(pData[point].cfpfh[i] - pModel[p].cfpfh[i]));
				}	
			}
			else if(cfpfh == 2) {
				for(int i = 0; i < 33; i++) {
					diff += (abs(pData[point].cfpfh[i] - pModel[p].cfpfh[i]));
				}	
			}
			else if(cfpfh == 3) {
				for(int i = 33; i < 41; i++) {
					diff += (abs(pData[point].cfpfh[i] - pModel[p].cfpfh[i]));
				}	
			}
			
            if(diff < minD) {
				minD = diff;
			}
        }
		d = minD;
	}

	return d;
}

/**
 * Computes sum of c-FPFH differences for BnB
 */
float GoICP::sumFPFH(float x, float y, float z) {
	float sum = 0;
	for (int i = 0; i < Nd; i++)
	{
		sum += computeFPFHDifference(false, i, pDataTemp[i].x + x, pDataTemp[i].y + y, pDataTemp[i].z + z);
	}
	sum = sum / Nd;
	return sum;
}
