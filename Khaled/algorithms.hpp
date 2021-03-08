#ifndef _ALGORITHMS_H_
#define _ALGORITHMS_H_

#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <iterator>
#include <ctime>
#include <cmath>
#include <map>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <stack>
#include <queue>
#include <string>
#include <sstream>
#include <random>
#include <limits>
#include <cassert>
#include <algorithm>
#include <parallel/algorithm>
#include <parallel/numeric>
//#include <parallel/random_shuffle>

#include "../CSR.h"

#include "../utility.h"

#include "../sample/commonutility.hpp"
#include "../sample/Coordinate.hpp"

#include "../sample/MortonCode.hpp"
#include "../sample/BarnesHut.hpp"

using namespace std;

#define VALUETYPE double
#define INDEXTYPE int
#define CACHEBLOCK 4
#define MAXMIN 3.0
#define MAX_SECTORS 9
#define t 0.99
#define PI 3.14159265358979323846
static int PROGRESS = 0;

#pragma omp declare reduction(plus:Coordinate<VALUETYPE>:omp_out += omp_in) initializer(omp_priv = Coordinate<VALUETYPE>(0.0, 0.0))

class vertex{
	public:
		INDEXTYPE v;
		VALUETYPE left;
		VALUETYPE right;
		vertex(INDEXTYPE a, VALUETYPE l, VALUETYPE r){
			v = a;
			left = l;
			right = r;
		}
		~vertex(){}
};

class algorithms{
	public:
		CSR<INDEXTYPE, VALUETYPE> graph;
		Coordinate<VALUETYPE>  *nCoordinates, *prevCoordinates;
		VALUETYPE *sectors;
		VALUETYPE K = 1.0, C = 1.0, Shift=1.0, init;
		VALUETYPE minX, minY, maxX, maxY, W, threshold, learningrate, maxV;
		string filename;
		string outputdir;
		string labelfile;
		string initfile;
		vector<string> labels;
		INDEXTYPE rootnode = -1;
		INDEXTYPE *childcount;
		INDEXTYPE *visitednodes;
	public:
	algorithms(CSR<INDEXTYPE, VALUETYPE> &A_csr, string input, string outputd, string lfile, int init, double weight, double lr, string ifile){
		graph.make_empty();
		graph = A_csr;
		outputdir = outputd;
		labelfile = lfile;
		initfile = ifile;
		//cout << initfile << endl;
		W = weight;
		filename = input;
		learningrate = lr;
		nCoordinates = static_cast<Coordinate<VALUETYPE> *> (::operator new (sizeof(Coordinate<VALUETYPE>[A_csr.rows])));
		prevCoordinates = static_cast<Coordinate<VALUETYPE> *> (::operator new (sizeof(Coordinate<VALUETYPE>[A_csr.rows])));
		this->init = init;
		childcount = static_cast<INDEXTYPE *> (::operator new (sizeof(INDEXTYPE[A_csr.rows])));
		visitednodes = static_cast<INDEXTYPE *> (::operator new (sizeof(INDEXTYPE[A_csr.rows])));
		for(INDEXTYPE i = 0; i < graph.rows; i++){
			visitednodes[i] = 0;
		}
	}
	~algorithms(){
		//delete nCoordinates;
		//delete prevCoordinates;
	}
	void randInit(){
		#pragma omp parallel for schedule(static)
		for(INDEXTYPE i = 0; i < graph.rows; i++){
			VALUETYPE x, y;
			//x = -1.0 + 2.0 * rand()/(RAND_MAX+1.0);
			//y = -1.0 + 2.0 * rand()/(RAND_MAX+1.0);
			do{
				x = -1.0 + 2.0 * rand()/(RAND_MAX+1.0);
				y = -1.0 + 2.0 * rand()/(RAND_MAX+1.0);
			}while(x * x + y * y > 1.0);
			x = x * MAXMIN;
			y = y * MAXMIN;
			
			nCoordinates[i] = Coordinate <VALUETYPE>(x, y, i);
		}
	}
	void countNumOfChildren(INDEXTYPE root){
		childcount[root] = 1;
		for(INDEXTYPE i = graph.rowptr[root]; i < graph.rowptr[root+1]; i++){
			if(!visitednodes[graph.colids[i]]){
				visitednodes[graph.colids[i]] = 1;
				countNumOfChildren(graph.colids[i]);
				childcount[root] += childcount[graph.colids[i]];
			}
		}
	}
	bool comp(Coordinate<VALUETYPE> A, Coordinate<VALUETYPE> B){
		return A.y > B.y;
	}
	void readlabels(){
		cout << "Reading Label File:" << labelfile << endl;
		ifstream infile(labelfile);
                if(infile.is_open()){
                        int index = 0;
                        string line;
                        while(getline(infile, line)){
				short int l = (short int)line.length();
                                nCoordinates[index].z = l;
				//printf("%d: length = %d\n", index, nCoordinates[index].z);
                                index++;
                        }
                }
                infile.close();
	}
	INDEXTYPE findRoot(){
		VALUETYPE degcen[graph.rows] = {0.0};
#if NV > 1500		
		 #pragma omp parallel for schedule(static)
#endif
		for(INDEXTYPE i = 0; i < graph.rows; i++){
			INDEXTYPE dist = 0;
			stack<INDEXTYPE> STACK;
			INDEXTYPE visited[graph.rows] = {0};
			STACK.push(i);
			STACK.push(0);
			visited[i] = 1;
			while(!STACK.empty()){
				INDEXTYPE d = STACK.top();
				STACK.pop();
				INDEXTYPE q = STACK.top();
				STACK.pop();
				dist += d;
				//printf("Node=%d, target = %d, dist = %d\n", i, q, dist);
				for(INDEXTYPE j = graph.rowptr[q]; j < graph.rowptr[q+1]; j++){
					if(!visited[graph.colids[j]]){
						STACK.push(graph.colids[j]);
						STACK.push(d+1);
						visited[graph.colids[j]] = 1;
					}
				}
			}
			degcen[i] = (1.0 * graph.rows) / (1.0 * dist);
			//printf("FRoot:%d = %lf, dist = %d\n", i, degcen[i], dist);
		}
		INDEXTYPE ret = 0;
		VALUETYPE val = degcen[0];
		for(INDEXTYPE i = 1; i < graph.rows; i++){
			if(val < degcen[i]){
				val = degcen[i];
				ret = i;
			}
		}
		return ret;
	}

	void initDFS(INDEXTYPE root){
                int visited[graph.rows] = {0};
                stack <vertex> STACK;
		double scalefactor = 1.0;
		minX = minY = 99.0;//numeric_limits<double>::max();
		maxX = maxY = 0.0;//numeric_limits<double>::min();
		maxV = maxX;
		vertex ROOT(root, 0.0, 360.0);
                STACK.push(ROOT);
		double radi = 50;
                visited[root] = 1;
                nCoordinates[root] = Coordinate <VALUETYPE>(0.0, 0.0);
		while(!STACK.empty()){
                        auto node = STACK.top();
                        STACK.pop();
			//printf("Heree..id(%d) = %lf, %lf\n", node.v, node.left, node.right);
			INDEXTYPE NN = graph.rowptr[node.v+1] - graph.rowptr[node.v];
                        //printf("Heree..deg(%d) = %d\n", node.v, NN);
			if(NN > 0){
				radi = NN;
                                VALUETYPE deg = (node.right - node.left) / childcount[node.v];
                                VALUETYPE degree = node.left + deg / 2.0;
                                for(INDEXTYPE n = graph.rowptr[node.v]; n < graph.rowptr[node.v+1]; n++){
                                        if(visited[graph.colids[n]] == 0){
						auto ldeg = degree;
						auto rdeg = degree + deg * childcount[graph.colids[n]];
						//printf("Node ID:%d, degree = %lf, ldeg = %lf, rdeg = %lf\n", graph.colids[n], degree, ldeg, rdeg);
                                                nCoordinates[graph.colids[n]] = Coordinate <VALUETYPE>(nCoordinates[node.v].getX() + radi*cos(PI*(degree)/180.0), nCoordinates[node.v].getY() + radi*sin(PI*(degree)/180.0));
						visited[graph.colids[n]] = 1;
						vertex temp = vertex(graph.colids[n], ldeg, rdeg);
                                                STACK.push(temp);
                                                degree += deg * childcount[graph.colids[n]];
                                        }
                                }
                        }
                }
		scalefactor = 1.0;//2.0 * MAXMIN / max(maxX - minX, maxY - minY);
		//#pragma omp parallel for schedule(static)
		for(int i = 0; i < graph.rows; i++){
			nCoordinates[i] = nCoordinates[i] * scalefactor;
			maxX = max(maxX, nCoordinates[i].x);
			minX = min(minX, nCoordinates[i].x);
			maxY = max(maxY, nCoordinates[i].y);
			minY = min(minY, nCoordinates[i].y);
			maxV = max(maxV, fabs(maxX));
			maxV = max(maxV, fabs(maxY));
			maxV = max(maxV, fabs(minX));
			maxV = max(maxV, fabs(minY));
		}
        }
	INDEXTYPE randIndex(INDEXTYPE max_num, INDEXTYPE min_num){
        	const INDEXTYPE newIndex = (rand() % (max_num - min_num)) + min_num;
        	return newIndex;
	}
	VALUETYPE scaleX(VALUETYPE x){
		if(x > maxX) return maxX;
		else if(x < minX) return minX;
		else return x;
	}
	VALUETYPE scaleY(VALUETYPE y){
                if(y > maxY) return maxY;
                else if(y < minY) return minY;
                else return y;
        }
	VALUETYPE getAngle(Coordinate<VALUETYPE> v){
		VALUETYPE arad = atan2(v.y, v.x);
		return arad;
	}

	VALUETYPE scaleXcircular(VALUETYPE x, VALUETYPE a){
		VALUETYPE tempX = maxV * cos(a);
                if(x > tempX) return tempX;
                else if(x < -tempX) return -tempX;
                else return x;
        }
	VALUETYPE scaleYcircular(VALUETYPE y, VALUETYPE a){
		VALUETYPE tempY = maxV * sin(a);
                if(y > tempY) return tempY;
                else if(y < tempY) return -tempY;
                else return y;
        }
	void fileInitialization()
	{
		VALUETYPE x, y;
                INDEXTYPE i;
		FILE *infile;
		infile = fopen(initfile.c_str(), "r");
		if(infile == NULL){
			cout << "ERROR in input coordinates file!\n" << endl;
			exit(1);
		}else{
			int index = 0;
			char line[256];
			VALUETYPE x, y;
			INDEXTYPE i;
			while(fgets(line, 256, infile)){
				sscanf(line, "%lf %lf %d", &x, &y, &i);
				nCoordinates[index] = Coordinate <VALUETYPE>(x, y, i); 
				index++;
			}
		}
		fclose(infile);
	}
	Coordinate<VALUETYPE> calcAttraction(INDEXTYPE i, INDEXTYPE j){
                Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
		for(INDEXTYPE n = graph.rowptr[i]; n < graph.rowptr[i+1]; n++){
			f = f + (this->nCoordinates[graph.colids[n]] - this->nCoordinates[i]) * (W * (this->nCoordinates[graph.colids[n]] - this->nCoordinates[i]).getMagnitude()) - calcRepulsion(i, graph.colids[n]);
		}
                return f;
        }

	Coordinate<VALUETYPE> calcRepulsion(INDEXTYPE i, INDEXTYPE n){
		Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
                if((this->nCoordinates[n] - this->nCoordinates[i]).getMagnitude2() > 0)
                	f = f - (this->nCoordinates[n] - this->nCoordinates[i]) * (1.0 /(this->nCoordinates[n] - this->nCoordinates[i]).getMagnitude2());
		return f;
	}
	
	VALUETYPE updateStepLength(VALUETYPE STEP, VALUETYPE ENERGY, VALUETYPE ENERGY0){
		if(ENERGY < ENERGY0){
                        PROGRESS = PROGRESS + 1;
                        if(PROGRESS >= 5){
                                PROGRESS = 0;
                                STEP = STEP / t;
                        }
                }else{
                        PROGRESS = 0;
                        STEP = t * STEP;
                }
		return STEP;
	}
	
	vector<VALUETYPE> cacheBlockingminiBatchForceDirectedAlgorithmSD(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE, int flag = 0){
		INDEXTYPE LOOP = 0;
                INDEXTYPE blocky = 512, blockx = 2;
                VALUETYPE start, end, ENERGY, ENERGY0;
                VALUETYPE STEP = 1.0;
                vector<VALUETYPE> result;
                vector<INDEXTYPE> indices;
                vector<int> kindex(graph.rows, 0);
                for(INDEXTYPE i = 0; i < graph.rows; i++) indices.push_back(i);
                ENERGY0 = ENERGY = numeric_limits<VALUETYPE>::max();
                omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
                if(flag == 0){
                if(init == 0){
                        initDFS(0);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }}else{
                        STEP = pow(0.999, 4 * ITERATIONS);
                }
		while(LOOP < ITERATIONS){
			ENERGY0 = ENERGY;
                        ENERGY = 0;
                        #pragma omp parallel for simd
                        for(INDEXTYPE k = 0; k < graph.rows; k++){
                                prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
                        }
                        for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static)
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i += 1){
                                        if(i >= graph.rows)continue;
					//#pragma omp simd
					for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1]; j += 1){
						int v = graph.colids[j];
						if(((this->nCoordinates[v] - this->nCoordinates[i]).getMagnitude2()) > 0)
						prevCoordinates[i] += (this->nCoordinates[v] - this->nCoordinates[i]) * (W * (this->nCoordinates[v] - this->nCoordinates[i]).getMagnitude()) + (this->nCoordinates[v] - this->nCoordinates[i]) * (1.0 / ((this->nCoordinates[v] - this->nCoordinates[i]).getMagnitude2()));
						else
						prevCoordinates[i] += (this->nCoordinates[v] - this->nCoordinates[i]) * (W * (this->nCoordinates[v] - this->nCoordinates[i]).getMagnitude());
					}
					Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
					//#pragma omp simd
                                        for(INDEXTYPE j = 0; j < i; j += 1){
						if(((this->nCoordinates[j] - this->nCoordinates[i]).getMagnitude2()) > 0)
						f += (this->nCoordinates[j] - this->nCoordinates[i]) * (1.0 / ((this->nCoordinates[j] - this->nCoordinates[i]).getMagnitude2()));
                                        }
					//#pragma omp simd
					for(INDEXTYPE j = i+1; j < graph.rows; j += 1){
						if(((this->nCoordinates[j] - this->nCoordinates[i]).getMagnitude2()) > 0)
						f += (this->nCoordinates[j] - this->nCoordinates[i]) * (1.0 / ((this->nCoordinates[j] - this->nCoordinates[i]).getMagnitude2()));
					}
					prevCoordinates[i] = prevCoordinates[i] - f;
                                }
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        nCoordinates[i] = nCoordinates[i] + prevCoordinates[i].getUnitVector() * STEP;
                                        ENERGY += prevCoordinates[i].getMagnitude2();
                                }
                        }
			STEP = STEP * 0.999;
                        LOOP++;
                }
                end = omp_get_wtime();
                if(flag == 0){
                cout << "Cache BlockingSD Minibatch Size:" << BATCHSIZE  << endl;
                cout << "Cache BlockingSD Minbatch Energy:" << ENERGY << endl;
                cout << "Cache BlockingSD Minibatch Parallel Wall time required:" << end - start << endl;
                writeToFile("CACHESDMINB"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(LOOP));
                }
                result.push_back(ENERGY);
                result.push_back(end - start);
		return result;
	}
	bool checkTriples(Coordinate<VALUETYPE> A, Coordinate<VALUETYPE> B, Coordinate<VALUETYPE> C){
		return (C.y-A.y) * (B.x-A.x) >= (B.y-A.y) * (C.x-A.x);
	}
	bool checkIntersection(Coordinate<VALUETYPE> A, Coordinate<VALUETYPE> B, Coordinate<VALUETYPE> C, Coordinate<VALUETYPE> D){
		return checkTriples(A,C,D) != checkTriples(B,C,D) and checkTriples(A,B,C) != checkTriples(A,B,D);
	}

	bool checkOverlap(Coordinate<VALUETYPE> A, Coordinate<VALUETYPE> B){
		VALUETYPE offsetA = A.z * 0.45;//A.z * 1.15;
		VALUETYPE offsetB = B.z * 0.45;///2.0;//B.z * 1.15;
		VALUETYPE h = 1.15;
    		if((A.x - offsetA < B.x + offsetB && B.x - offsetB < A.x + offsetA) && (A.y + h > B.y - h && A.y - h < B.y + h))
			return true;
		else
			return false;
	}

	bool checkCrossing(INDEXTYPE LOOP){
		bool flag = false;
		#pragma omp parallel
		{
			INDEXTYPE tid = omp_get_thread_num();
			INDEXTYPE nthreads = omp_get_num_threads();
			INDEXTYPE perthreadwork = graph.rows / nthreads;
			INDEXTYPE starti = tid * perthreadwork;
 			INDEXTYPE endi = (tid + 1) * perthreadwork;
			if(tid == nthreads - 1) endi = max(endi, graph.rows);
			for(INDEXTYPE i = starti; i < endi && flag == false; i++){
                		for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1] && flag == false; j++){
                        		for(INDEXTYPE k = 0; k < graph.rows && flag == false; k++){
                                		for(INDEXTYPE l = graph.rowptr[k]; l < graph.rowptr[k+1]; l++){
                                        		if(i == k || graph.colids[j] == graph.colids[l]) continue;
                                                	if(i == graph.colids[l] || k == graph.colids[j]) continue;
                                                	if(checkIntersection(nCoordinates[i], nCoordinates[graph.colids[j]], nCoordinates[k], nCoordinates[graph.colids[l]])){
                                                        	flag = true;//return true;//exit(0);
                                                	}
                                        	}
                                	}
                        	}
                	}
		}
		return flag;
	}
	bool checkCrossingSeq(INDEXTYPE LOOP){
                bool flag = false;
                for(INDEXTYPE i = 0; i < graph.rows && flag == false; i++){
                	for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1] && flag == false; j++){
                        	for(INDEXTYPE k = 0; k < graph.rows && flag == false; k++){
                                	for(INDEXTYPE l = graph.rowptr[k]; l < graph.rowptr[k+1]; l++){
                                        	if(i == k || graph.colids[j] == graph.colids[l]) continue;
                                                if(i == graph.colids[l] || k == graph.colids[j]) continue;
                                                if(checkIntersection(nCoordinates[i], nCoordinates[graph.colids[j]], nCoordinates[k], nCoordinates[graph.colids[l]])){                                                    
                                                	flag = true;//return true;//exit(0);
                                                }
                                        }
                                }
                        }
                }
                return flag;
        }
	bool hasEdgeCrossing(INDEXTYPE LOOP, INDEXTYPE n, Coordinate<VALUETYPE> N){
		bool flag = false;
		#pragma omp parallel
		{
			INDEXTYPE tid = omp_get_thread_num();
                        INDEXTYPE nthreads = omp_get_num_threads();
                        INDEXTYPE perthreadwork = graph.rows / nthreads;
                        INDEXTYPE starti = tid * perthreadwork;
                        INDEXTYPE endi = (tid + 1) * perthreadwork;
                        if(tid == nthreads - 1) endi = max(endi, graph.rows);
			for(int i = starti; i < endi && flag == false; i++){
				for(int j = graph.rowptr[i]; j < graph.rowptr[i+1]; j++){
					for(int f = graph.rowptr[n]; f < graph.rowptr[n+1]; f++){
						if(n == i || graph.colids[f] == graph.colids[j]) continue;
						if(graph.colids[f] == i || n == graph.colids[j]) continue;
						if(checkIntersection(N, nCoordinates[graph.colids[f]], nCoordinates[i], nCoordinates[graph.colids[j]])){
							//printf("Test:(%d %d), (%d, %d)\n", n, graph.colids[f], i, graph.colids[j]);
							flag = true;
						}
					}
				}
			}
		}
		return flag;
	}
	void rescaleLayout(INDEXTYPE maxloop = 25, INDEXTYPE BATCHSIZE=128, VALUETYPE STEP=0.6){
		bool flag = true;
		int loop = 0;
		VALUETYPE dedgelength = 40.0;
		printf("Running rescaling function...\n");
		while(flag){
			flag = false;
			for(INDEXTYPE k = 0; k < graph.rows; k++){
                                prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
                        }
 			for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static)
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i += 1){
                                        if(i >= graph.rows)continue;
                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
					for(INDEXTYPE j = 0; j < graph.rows; j++){
						if(i == j) continue;
						/*for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1]; j++){
                                                	INDEXTYPE colj = graph.colids[j];
                                                	auto forceDiff = nCoordinates[colj] - nCoordinates[i];
                                                	auto attrc = forceDiff.getMagnitude2();
                                                	if(sqrt(attrc) > dedgelength)
                                                        	f += forceDiff * sqrt(attrc);
						}*/
						if(checkOverlap(nCoordinates[i], nCoordinates[j])){
							VALUETYPE dist = (this->nCoordinates[j] - this->nCoordinates[i]).getMagnitude2();
							if(dist > 0.0)
							f = f - (this->nCoordinates[j] - this->nCoordinates[i]) * (1.0 / (dist));	
							flag = true;
						}
						
					}
					prevCoordinates[i] += f;
				}
				//#pragma omp parallel for
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        auto p = Coordinate <VALUETYPE>(0.0, 0.0);
                                        p.x = nCoordinates[i].x;
                                        p.y = nCoordinates[i].y;
					if(prevCoordinates[i].getMagnitude2() > 0.0)
					nCoordinates[i] = nCoordinates[i] + prevCoordinates[i].getUnitVector() * STEP;
                                        if(hasEdgeCrossing(loop, i, nCoordinates[i])){
                                                nCoordinates[i] = p;
                                        }
					maxV = max(maxV, fabs(nCoordinates[i].x));
                                        maxV = max(maxV, fabs(nCoordinates[i].y));
                                	prevCoordinates[i] = Coordinate <VALUETYPE>(0.0, 0.0);
				}
			}
			if(loop > maxloop) break;
			loop++;
		}
		printf("Final Loop: %d\n", loop);
	}
	vector<VALUETYPE> cacheBlockingminiBatchForceDirectedAlgorithm2(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE, int flag = 0){
                INDEXTYPE LOOP = 0;
                INDEXTYPE blocky = 512, blockx = 2;
		VALUETYPE dedgelength = 10.0;
                VALUETYPE start, end, ENERGY, ENERGY0;
                VALUETYPE STEP = learningrate;
                vector<VALUETYPE> result;
                vector<INDEXTYPE> indices;
                vector<int> kindex(graph.rows, 0);
                VALUETYPE delta = 10.0;
		if(labelfile.size() > 0) readlabels();
                printf("Calling findRoot...\n");
                INDEXTYPE root = findRoot();
                printf("Root = %d\n", root);
                printf("Calling ChildCount..\n");
                visitednodes[root] = 1;
                countNumOfChildren(root);
		if(delta < 1.0){
                        printf("Normalizing...\n");
                        for(INDEXTYPE k = 0; k < graph.rows; k++){
                                for(INDEXTYPE j = graph.rowptr[k]; j < graph.rowptr[k+1]; j++){
                                        delta += (nCoordinates[k] - nCoordinates[graph.colids[j]]).getMagnitude();
                                }
                        }
                        delta = delta / graph.rows;
                }
                ENERGY0 = numeric_limits<VALUETYPE>::max();
                ENERGY = 0;
                INDEXTYPE *twoends = static_cast<INDEXTYPE *> (::operator new (sizeof(INDEXTYPE[graph.rows*2])));
                VALUETYPE gamma = 4 * delta;
                omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
                printf("Number of threads = %d\n", omp_get_num_threads());
                if(flag == 0){
                if(init == 0){
                        initDFS(root);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }}else{
                        STEP = pow(0.999, 4 * ITERATIONS);
                }
                VALUETYPE angle = getAngle(Coordinate <VALUETYPE>(-10.0, -10.0));
                printf("P1(0, %lf), P2(10.0, 10.0) = Angle: %lf\n", maxY, angle * 180.0 / PI);
                printf("After initialization...\n");
                if(!checkCrossing(-1))printf("No Crossing After initialization...\n");
		while(LOOP < ITERATIONS){
                        ENERGY0 = ENERGY;
                        ENERGY = 0;
                        #pragma omp parallel for simd proc_bind(close)
                        for(INDEXTYPE k = 0; k < graph.rows; k++){
                                kindex[k] = graph.rowptr[k];
                                prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
                        }
                        for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static)
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i += 1){
                                        if(i >= graph.rows)continue;
                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
					for(INDEXTYPE j = 0; j < graph.rows; j += 1){
                                                VALUETYPE dist = (this->nCoordinates[j] - this->nCoordinates[i]).getMagnitude2();
                                                if(dist > 0){
                                                        f = f - (this->nCoordinates[j] - this->nCoordinates[i]) * ((delta * delta) / (dist));
                                                }
						if(checkOverlap(nCoordinates[i], nCoordinates[j]) && dist > 0){
							f = f - (this->nCoordinates[j] - this->nCoordinates[i]) * ((delta * delta) / (dist));
						}
                                        }
					
					for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1]; j += 1){
						VALUETYPE dist = (this->nCoordinates[j] - this->nCoordinates[i]).getMagnitude();
                                                if(dist < dedgelength && dist > 0.0)
							f = f - (nCoordinates[graph.colids[j]] - nCoordinates[i]) * ((1.0 / delta) * ( (delta * delta) / (dist)));
                                        	else if(dist > 0.0)
							f += (nCoordinates[graph.colids[j]] - nCoordinates[i]) * ((0.05 / delta) * dist);
					}
					
					prevCoordinates[i] += f;
				}
				for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
					auto p = Coordinate <VALUETYPE>(0.0, 0.0);
					p.x = nCoordinates[i].x;
					p.y = nCoordinates[i].y;
					nCoordinates[i].x = nCoordinates[i].x + STEP * prevCoordinates[i].x;
                                        nCoordinates[i].y = nCoordinates[i].y + STEP * prevCoordinates[i].y;
					//if(hasEdgeCrossing(LOOP, i, nCoordinates[i])){
					if(checkCrossing(LOOP)){
						//printf("i:%d = (%lf, %lf), j:(%lf, %lf)\n", i, p.x, p.y, nCoordinates[i].x, nCoordinates[i].y);
						nCoordinates[i] = p;
					}
					maxV = max(maxV, fabs(nCoordinates[i].x));
					maxV = max(maxV, fabs(nCoordinates[i].y));
				}
				if(checkCrossing(LOOP)){
					printf("Dead End!\n");
					exit(0);
				}
			}
			//STEP = STEP * 0.9;
			LOOP++;
		}
		end = omp_get_wtime();
		rescaleLayout(105, BATCHSIZE, STEP);
                if(flag == 0){
                cout << "Cache Blocking Minibatch Size:" << BATCHSIZE  << endl;
                cout << "Cache Blocking Minbatch Energy:" << ENERGY << endl;
                cout << "Cache Blocking Minibatch Parallel Wall time required:" << end - start << endl;
                writeToFile("Batch2PrEd"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(LOOP));
                }
                result.push_back(ENERGY);
                result.push_back(end - start);
                return result;
	}

        bool doOverlap(VALUETYPE l1_x, VALUETYPE l1_y, VALUETYPE r1_x,  VALUETYPE r1_y, VALUETYPE l2_x, VALUETYPE l2_y, VALUETYPE r2_x, VALUETYPE r2_y) 
        { 
          // If one rectangle is on left side of other 
          if (l1_x >= r2_x || l2_x >= r1_x) 
            return false; 
  
          // If one rectangle is above other 
          if (l1_y <= r2_y || l2_y <= r1_y) 
            return false; 
  
          return true; 
        }

        void postProcessing()
        {
          printf("inside postProcessing\n");
          VALUETYPE x_mx = nCoordinates[0].x, x_mn = nCoordinates[0].x, y_mx = nCoordinates[0].y, y_mn = nCoordinates[0].y;
          VALUETYPE total_area, scale, total_len = 0;
          for(INDEXTYPE i=0;i<graph.rows;i++)
          {
            if(x_mx<nCoordinates[i].x)x_mx = nCoordinates[i].x;
            if(x_mn>nCoordinates[i].x)x_mn = nCoordinates[i].x;
            if(y_mx<nCoordinates[i].y)y_mx = nCoordinates[i].y;
            if(y_mn>nCoordinates[i].y)y_mn = nCoordinates[i].y;
            total_len += nCoordinates[i].z;
          }
          VALUETYPE cntr_x = (x_mx-x_mn)/2;
          VALUETYPE cntr_y = (y_mx-y_mn)/2;
          for(INDEXTYPE i=0;i<graph.rows;i++)
          {
            nCoordinates[i].x = nCoordinates[i].x-cntr_x;
            nCoordinates[i].y = nCoordinates[i].y-cntr_y;
          }
          for(INDEXTYPE i=0;i<graph.rows;i++)
          {
            //nCoordinates[i].x = nCoordinates[i].x*2500/(x_mx-x_mn);
            //nCoordinates[i].y = nCoordinates[i].y*2500/(y_mx-y_mn);
            nCoordinates[i].x = nCoordinates[i].x*500/(x_mx-x_mn);
            nCoordinates[i].y = nCoordinates[i].y*500/(y_mx-y_mn);
          }
          x_mx = nCoordinates[0].x, x_mn = nCoordinates[0].x, y_mx = nCoordinates[0].y, y_mn = nCoordinates[0].y;
          for(INDEXTYPE i=0;i<graph.rows;i++)
          {
            if(x_mx<nCoordinates[i].x)x_mx = nCoordinates[i].x;
            if(x_mn>nCoordinates[i].x)x_mn = nCoordinates[i].x;
            if(y_mx<nCoordinates[i].y)y_mx = nCoordinates[i].y;
            if(y_mn>nCoordinates[i].y)y_mn = nCoordinates[i].y;
          }
          total_area = (x_mx-x_mn)*(y_mx-y_mn);
          //scale = 0.01*total_area/(0.6*total_len);
          //scale = 0.1*total_area/(0.6*total_len);
          scale = sqrt(0.1*total_area/(0.6*total_len));
          //cout << "width " << scale*.6 << " height " << scale << endl;
          INDEXTYPE count = 0, count_solved = 0;
          for(INDEXTYPE i=0;i<graph.rows;i++)
          {
            for(INDEXTYPE j=i+1;j<graph.rows;j++)
            {
              // check overlap
              VALUETYPE l1_x = nCoordinates[i].x, l1_y = nCoordinates[i].y+scale, r1_x = nCoordinates[i].x+nCoordinates[i].z*scale*.6, r1_y = nCoordinates[i].y, l2_x = nCoordinates[j].x, l2_y = nCoordinates[j].y+scale, r2_x = nCoordinates[j].x+nCoordinates[j].z*scale*.6, r2_y = nCoordinates[j].y;
              if(doOverlap(l1_x, l1_y, r1_x, r1_y, l2_x, l2_y, r2_x, r2_y))
              {
                INDEXTYPE overlapped_labels [] = {i, j};
                bool overlap_crossing_free = false;
                for(INDEXTYPE k=0;k<2;k++)
                {
                  INDEXTYPE ind = overlapped_labels[k];
                  // continue this for 100 times
                  //for(INDEXTYPE l=0;l<1200;l++)
                  //for(INDEXTYPE l=0;l<600;l++)
                  for(INDEXTYPE l=0;l<300;l++)
                  //The following one is the default
                  //for(INDEXTYPE l=0;l<100;l++)
                  //for(INDEXTYPE l=0;l<1;l++)
                  {
                    // consider a small square of 300x300
                    // random sample a number from -150 to 150, once for x coord, once for y coord
                    //VALUETYPE box_size = 1201;
                    //VALUETYPE box_size = 301;
                    // The following one is the default one
                    VALUETYPE box_size = 601;
                    //VALUETYPE box_size = 100;
                    //VALUETYPE box_size = 50;
                    //VALUETYPE box_size = 10;
                    VALUETYPE shift_x = (((double) rand() / RAND_MAX) * box_size) - (box_size/2);
                    VALUETYPE shift_y = (((double) rand() / RAND_MAX) * box_size) - (box_size/2);
                    VALUETYPE prev_x = nCoordinates[ind].x;
                    VALUETYPE prev_y = nCoordinates[ind].y;
                    nCoordinates[ind].x += shift_x;
                    nCoordinates[ind].y += shift_y;
                    // check whether it removes overlap without new crossing
                    bool overlap_free = true;
                    for(INDEXTYPE m=0;m<graph.rows;m++)
                    {
                      if(m==ind)continue;
                      l1_x = nCoordinates[ind].x, l1_y = nCoordinates[ind].y+scale, r1_x = nCoordinates[ind].x+nCoordinates[ind].z*scale*.6, r1_y = nCoordinates[ind].y, l2_x = nCoordinates[m].x, l2_y = nCoordinates[m].y+scale, r2_x = nCoordinates[m].x+nCoordinates[m].z*scale*.6, r2_y = nCoordinates[m].y;
                      if(doOverlap(l1_x, l1_y, r1_x, r1_y, l2_x, l2_y, r2_x, r2_y))
                      {
                        overlap_free = false;
                        break;
                      }
                    }
                    if(overlap_free)
                    {
                      bool introducesCrossing = hasEdgeCrossing(0, ind, nCoordinates[ind]);
                      if(overlap_free && (!introducesCrossing))
                      //if(overlap_free)
                      {
                        overlap_crossing_free = true;
                        cout << "Found an overlap and crossing free coordinate for [" << overlapped_labels[0] << ", " << overlapped_labels[1] << ']' << endl;
                        count_solved += 1;
                      }
                    }
                    if(overlap_crossing_free)break;
                    else
                    {
                      nCoordinates[ind].x = prev_x;
                      nCoordinates[ind].y = prev_y;
                    }
                  }
                  if(overlap_crossing_free)break;
                  else
                  {
                    cout << "Could not found an overlap and crossing free coordinate for [" << overlapped_labels[0] << ", " << overlapped_labels[1] << ']' << endl;
                  }
                }
                count += 1;
              }
            }
          }
          cout << "Post-processing- overlaps: " << count << " removed: " << count_solved << endl;
        }

	vector<VALUETYPE> cacheBlockingminiBatchForceDirectedAlgorithm(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE, INDEXTYPE ns, VALUETYPE lr, VALUETYPE lrforlo, INDEXTYPE ITER){
        	INDEXTYPE LOOP = 0, DIM = 2;
        	Coordinate<VALUETYPE>  *samples;
		samples = static_cast<Coordinate<VALUETYPE> *> (::operator new (sizeof(Coordinate<VALUETYPE>[ns])));
		VALUETYPE STEP = learningrate;
                vector<VALUETYPE> result;
		VALUETYPE start, end, ENERGY, ENERGY0;
                vector<INDEXTYPE> indices;
		for(INDEXTYPE i = 0; i < graph.rows; i++) indices.push_back(i);
		if(labelfile.size() > 0) readlabels();
                printf("From Option 2: Calling findRoot...\n");
                INDEXTYPE root = findRoot();
                printf("Root = %d\n", root);
                printf("Calling ChildCount..\n");
                visitednodes[root] = 1;
                countNumOfChildren(root);
                ENERGY0 = numeric_limits<VALUETYPE>::max();
                ENERGY = 0;
                omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
                printf("Number of threads = %d\n", omp_get_num_threads());
		initDFS(root);
                VALUETYPE angle = getAngle(Coordinate <VALUETYPE>(-10.0, -10.0));
                //printf("P1(0, %lf), P2(10.0, 10.0) = Angle: %lf\n", maxY, angle * 180.0 / PI);
                printf("After initialization...\n");
                if(!checkCrossing(-1))printf("No Crossing After initialization...\n");

		#pragma omp parallel for simd proc_bind(close)
		for(INDEXTYPE k = 0; k < graph.rows; k++){
                	prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
                }
		VALUETYPE dedgelength = 20.0;
	        while(LOOP < ITERATIONS){
        	        ENERGY0 = ENERGY;
                	ENERGY = 0;
                	for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                        	INDEXTYPE baseindex = b * BATCHSIZE;
                        	for(INDEXTYPE s = 0; s < ns; s++){
                                	INDEXTYPE rindex = randIndex(graph.rows-1, 0);
                                	samples[s] = nCoordinates[rindex];
                        	}
                        	#pragma omp parallel for schedule(static)
				for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
					Coordinate<VALUETYPE> forceDiff = Coordinate <VALUETYPE>(0.0, 0.0);
                                        INDEXTYPE k = graph.rowptr[indices[i]];
					
					for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1]; j++){
						INDEXTYPE colj = graph.colids[j];
						forceDiff = nCoordinates[colj] - nCoordinates[i];
						auto attrc = forceDiff.getMagnitude2();
						if(sqrt(attrc) > dedgelength)
							f += forceDiff * sqrt(attrc);
					}
                                        for(INDEXTYPE j = 0; j < ns; j++){
                                        	forceDiff = samples[j] - nCoordinates[i];
						auto repuls = forceDiff.getMagnitude2();
						if(repuls > 0.0) f -= forceDiff * (1.0 / repuls); 
					}
                                        prevCoordinates[indices[i]] += f;
                                }
				//#pragma omp parallel for
                		for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
					auto p = Coordinate <VALUETYPE>(0.0, 0.0);
                                        p.x = nCoordinates[i].x;
                                        p.y = nCoordinates[i].y;
                                        nCoordinates[indices[i]] += prevCoordinates[indices[i]].getUnitVector() * STEP;
					if(hasEdgeCrossing(LOOP, i, nCoordinates[indices[i]])){
						nCoordinates[indices[i]] = p;
					}else{
                                        	ENERGY += prevCoordinates[indices[i]].getMagnitude2();
                                	}
					maxV = max(maxV, fabs(nCoordinates[i].x));
                                        maxV = max(maxV, fabs(nCoordinates[i].y));
					prevCoordinates[indices[i]] = Coordinate <VALUETYPE>(0.0, 0.0);
				}
			}
                	STEP = STEP * 0.999;
                	LOOP++;
			//printf("Loop Count: %d\n", LOOP);
        	}
		if(checkCrossing(LOOP)){
                	printf("(after forceupdate) Dead End!\n");
                }else{
			printf("No edge-crossing after force-updates!\n");
		}
		rescaleLayout(ITER, BATCHSIZE, lrforlo);
		if(checkCrossing(LOOP)){
                        printf("(after label attachment) Dead End!\n");
                }
                postProcessing();
        	end = omp_get_wtime();
        	cout << "BatchPrEL Parallel Wall time required:" << end - start << endl;
        	result.push_back(end - start);
        	writeToFile("BatchPrEL"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(LOOP));
		return result;
	}
	vector<VALUETYPE> cacheBlockingminiBatchForceDirectedAlgorithmConverged(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE, int flag = 0){
                INDEXTYPE LOOP = 0;
                INDEXTYPE blocky = 512, blockx = 2;
                VALUETYPE start, end, ENERGY, ENERGY0;
                VALUETYPE STEP = 1.0;
                vector<VALUETYPE> result;
                vector<INDEXTYPE> indices;
                vector<int> kindex(graph.rows, 0);
		ENERGY0 = numeric_limits<VALUETYPE>::max();
                ENERGY = 0;
                omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
                if(flag == 0){
                if(init == 0){
                        initDFS(0);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }}else{
                        STEP = pow(0.999, 4 * ITERATIONS);
                }
                while((fabs(ENERGY0 - ENERGY)/ENERGY0 > threshold)){
                        ENERGY0 = ENERGY;
                        ENERGY = 0;
                        #pragma omp parallel for simd proc_bind(close)
                        for(INDEXTYPE k = 0; k < graph.rows; k++){
                                kindex[k] = graph.rowptr[k];
                                prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
                        }
                        for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static)
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i += blockx){
                                        if(i >= graph.rows)continue;
                                        for(INDEXTYPE j = 0; j < graph.rows; j += blocky){
                                                for(INDEXTYPE bi = 0; bi < blockx && i + bi < (b + 1) * BATCHSIZE; bi++){
                                                        if(i+bi >= graph.rows) break;
                                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
                                                        for(INDEXTYPE bj = 0; bj < blocky && j + bj < graph.rows; bj++){
                                                                if(j + bj == graph.colids[kindex[bi+i]]){
                                                                        f += (nCoordinates[j+bj] - nCoordinates[i+bi]) * (W * (nCoordinates[j+bj] - nCoordinates[i+bi]).getMagnitude());
                                                                        if(kindex[bi+i] < graph.rowptr[i+bi+1] - 1){
                                                                                kindex[bi+i]++;
                                                                        }
                                                                }else{
                                                                        VALUETYPE dist = (this->nCoordinates[j+bj] - this->nCoordinates[i+bi]).getMagnitude2();
                                                                        if(dist > 0)
                                                                        {
                                                                                f = f - (this->nCoordinates[j+bj] - this->nCoordinates[i+bi]) * (1.0 / (dist));
                                                                        }
                                                                }
                                                        }
                                                        prevCoordinates[i+bi] += f;
                                                }
                                        }
                                }
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        nCoordinates[i] = nCoordinates[i] + prevCoordinates[i].getUnitVector() * STEP;
                                        ENERGY += prevCoordinates[i].getMagnitude2();
                                }
                        }
                        STEP = STEP * 0.999;
                        LOOP++;
                }
		end = omp_get_wtime();
                if(flag == 0){
                cout << "Cache Blocking (converged) Minibatch Size:" << BATCHSIZE  << endl;
                cout << "Cache Blocking (converged)  Minbatch Energy:" << ENERGY << endl;
                cout << "Cache Blocking (converged) Minibatch Parallel Wall time required:" << end - start << endl;
                writeToFile("CACHEMINB"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(LOOP));
                }
                result.push_back(ENERGY);
                result.push_back(end - start);
                return result;
        }


	vector<VALUETYPE> LinLogcacheBlockingminiBatchForceDirectedAlgorithm(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE){
                INDEXTYPE LOOP = 0;
                INDEXTYPE blocky = 512, blockx = 2;
                VALUETYPE start, end, ENERGY, ENERGY0;
                VALUETYPE STEP = 1.0;
                vector<VALUETYPE> result;
                vector<INDEXTYPE> indices;
                vector<int> kindex(graph.rows, 0);
                for(INDEXTYPE i = 0; i < graph.rows; i++) indices.push_back(i);
                ENERGY0 = ENERGY = numeric_limits<VALUETYPE>::max();
                omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
                if(init == 0){
                        initDFS(0);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }
		while(LOOP < ITERATIONS){
                        ENERGY0 = ENERGY;
                        ENERGY = 0;
                        #pragma omp parallel for simd schedule(static)
                        for(INDEXTYPE k = 0; k < graph.rows; k++){
                                kindex[k] = graph.rowptr[k];
				prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
                        }
                        for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static)
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i += blockx){
                                        if(i >= graph.rows)continue;
                                        for(INDEXTYPE j = 0; j < graph.rows; j += blocky){
                                                for(INDEXTYPE bi = 0; bi < blockx && i + bi < (b + 1) * BATCHSIZE; bi++){
                                                        if(i+bi >= graph.rows) break;
                                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
                                                        for(INDEXTYPE bj = 0; bj < blocky && j + bj < graph.rows; bj++){
                                                                if(j + bj == graph.colids[kindex[bi+i]]){
                                                                        f += (nCoordinates[j+bj] - nCoordinates[i+bi]) * log2(1.0 + W * (nCoordinates[j+bj] - nCoordinates[i+bi]).getMagnitude());
                                                                        if(kindex[bi+i] < graph.rowptr[i+bi+1] - 1){
                                                                                kindex[bi+i]++;
                                                                        }
                                                                }else{
									VALUETYPE dist = (this->nCoordinates[j+bj] - this->nCoordinates[i+bi]).getMagnitude2();
                                                                        if(dist > 0)
                                                                        {
                                                                                f = f - (this->nCoordinates[j+bj] - this->nCoordinates[i+bi]) * (1.0 / dist);
                                                                        }
                                                                }
                                                        }
                                                        prevCoordinates[i+bi] += f;
                                                }
                                        }
                                }
				for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        nCoordinates[i] = nCoordinates[i] + prevCoordinates[i].getUnitVector() * STEP;
                                        ENERGY += prevCoordinates[i].getMagnitude2();
                                }
                        }
                        STEP = STEP * t;
                        LOOP++;
                }
                end = omp_get_wtime();
                cout << "LinLog Batch Size:" << BATCHSIZE << endl;
                cout << "LinLog Minbatch Energy:" << ENERGY << endl;
                result.push_back(ENERGY);
                cout << "LinLog Minibatch Parallel Wall time required:" << end - start << endl;
                result.push_back(end - start);
                writeToFile("LLCACHEMINB"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(ITERATIONS));
                return result;
        }

	vector<VALUETYPE> FAcacheBlockingminiBatchForceDirectedAlgorithm(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE){
                INDEXTYPE LOOP = 0;
                INDEXTYPE blocky = 512, blockx = 2;
                VALUETYPE start, end, ENERGY, ENERGY0;
                VALUETYPE STEP = 1.0;
                vector<VALUETYPE> result;
                vector<INDEXTYPE> indices;
                vector<int> kindex(graph.rows, 0);
                for(INDEXTYPE i = 0; i < graph.rows; i++) indices.push_back(i);
                ENERGY0 = ENERGY = numeric_limits<VALUETYPE>::max();
                omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
                if(init == 0){
                        initDFS(0);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }
                while(LOOP < ITERATIONS){
                        ENERGY0 = ENERGY;
                        ENERGY = 0;
                        #pragma omp parallel for simd schedule(static)
                        for(INDEXTYPE k = 0; k < graph.rows; k++){
                                kindex[k] = graph.rowptr[k];
                                prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
                        }
                        for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static)
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i += blockx){
                                        if(i >= graph.rows)continue;
                                        for(INDEXTYPE j = 0; j < graph.rows; j += blocky){
                                                for(INDEXTYPE bi = 0; bi < blockx && i + bi < (b + 1) * BATCHSIZE; bi++){
                                                        if(i+bi >= graph.rows) break;
                                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
                                                        for(INDEXTYPE bj = 0; bj < blocky && j + bj < graph.rows; bj++){
                                                                if(j + bj == graph.colids[kindex[bi+i]]){
                                                                        f += (nCoordinates[j+bj] - nCoordinates[i+bi]) * W;
                                                                        if(kindex[bi+i] < graph.rowptr[i+bi+1] - 1){
                                                                                kindex[bi+i]++;
                                                                        }
                                                                }else{
                                                                        VALUETYPE dist = (this->nCoordinates[j+bj] - this->nCoordinates[i+bi]).getMagnitude2();
									if(dist > 0)
                                                                        {
                                                                                f = f - (this->nCoordinates[j+bj] - this->nCoordinates[i+bi]) * (1.0 / dist);
                                                                        }
                                                                }
                                                        }
                                                        prevCoordinates[i+bi] += f;
                                                }
                                        }
                                }
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        nCoordinates[i] = nCoordinates[i] + prevCoordinates[i].getUnitVector() * STEP;
                                        ENERGY += prevCoordinates[i].getMagnitude2();
                                }
                        }
                        STEP = STEP * t;
                        LOOP++;
                }
                end = omp_get_wtime();
                cout << "FA Batch Size:" << BATCHSIZE << endl;
                cout << "FA Minbatch Energy:" << ENERGY << endl;
                result.push_back(ENERGY);
                cout << "FA Minibatch Parallel Wall time required:" << end - start << endl;
                result.push_back(end - start);
                writeToFile("FACACHEMINB"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(ITERATIONS));
                return result;
        }
		
	vector<VALUETYPE> BarnesHutApproximation(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE, VALUETYPE TH, int flag = 0){
                INDEXTYPE LOOP = 0;
		VALUETYPE start, end, ENERGY, ENERGY0;
		vector<VALUETYPE> result;
		VALUETYPE STEP = 1.0;
		ENERGY0 = ENERGY = numeric_limits<VALUETYPE>::max();
		Coordinate<VALUETYPE> *tempCoordinates = static_cast<Coordinate<VALUETYPE> *> (::operator new (sizeof(Coordinate<VALUETYPE>[graph.rows])));
                omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
		if(flag == 0){
		if(init == 0){
                        initDFS(0);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }}
		while(LOOP < ITERATIONS){
			ENERGY0 = ENERGY;
                        ENERGY = 0;
			#pragma omp parallel for simd schedule(static)
                        for(INDEXTYPE k = 0; k < graph.rows; k++){
				prevCoordinates[k] = Coordinate <VALUETYPE>(0.0, 0.0);
				tempCoordinates[k] = nCoordinates[k];
                        }
			//VALUETYPE s = omp_get_wtime();
			BarnesHut bh(tempCoordinates, graph.rows, TH);
			//VALUETYPE e = omp_get_wtime();
			//printf("BH Time: %lf\n", e - s);
			for(INDEXTYPE b = 0; b < (int)ceil( 1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static) 
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
					#pragma omp simd
                                        for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1]; j++){
						f = f + (nCoordinates[graph.colids[j]] - nCoordinates[i]) * (W * (nCoordinates[graph.colids[j]] - nCoordinates[i]).getMagnitude());
					}
					f = f - bh.calcRepForce(nCoordinates[i]);
					prevCoordinates[i] = f;
				}
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        nCoordinates[i] = nCoordinates[i] + prevCoordinates[i].getUnitVector() * STEP;
                                        ENERGY += prevCoordinates[i].getMagnitude2();
                                }
				
			}
			//STEP = updateStepLength(STEP, ENERGY, ENERGY0);
                        STEP = STEP * 0.999;
			LOOP++;
		}
		end = omp_get_wtime();
		if(flag == 0){
                cout << "Barnes-Hut Minibatch Size:" << BATCHSIZE << endl;
                cout << "Barnes-Hut Minbatch Energy:" << ENERGY << endl;
                result.push_back(ENERGY);
                cout << "Barnes-Hut Minibatch Parallel Wall time required:" << end - start << endl;
                result.push_back(end - start);
                writeToFile("BHMINB"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(ITERATIONS));
		}
                return result;
	}

	vector<VALUETYPE> approxForceDirectedAlgorithm(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE){
		INDEXTYPE LOOP = 0, approxITER = (int)(ITERATIONS * 0.8);
                VALUETYPE start, end, ENERGY, ENERGY0;
                VALUETYPE STEP = 1.0;
                vector<VALUETYPE> result;
		vector<INDEXTYPE> indices;
                for(INDEXTYPE i = 0; i < graph.rows; i++) indices.push_back(i);
                ENERGY0 = ENERGY = numeric_limits<VALUETYPE>::max();
		omp_set_num_threads(NUMOFTHREADS);
                start = omp_get_wtime();
                if(init == 0){
                        initDFS(0);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }
		while(LOOP < ITERATIONS){
                        ENERGY0 = ENERGY;
                        ENERGY = 0;
                        INDEXTYPE j;
                        INDEXTYPE k;
                        for(INDEXTYPE b = 0; b < (int)ceil(1.0 * graph.rows / BATCHSIZE); b += 1){
                                #pragma omp parallel for schedule(static)   
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        Coordinate<VALUETYPE> f = Coordinate <VALUETYPE>(0.0, 0.0);
					stack <int> STACKnode;
					if(LOOP > approxITER){
						INDEXTYPE k = graph.rowptr[indices[i]];
                                        	for(INDEXTYPE j = 0; j < graph.rows; j++){
                                                	if(j == graph.colids[k] && k < graph.nnz){
                                                        	f += (nCoordinates[j] - nCoordinates[indices[i]]) * (W * (nCoordinates[j] - nCoordinates[indices[i]]).getMagnitude());
                                                        	if(k < graph.rowptr[indices[i]+1]-1){
                                                                	k++;
                                                        	}
                                                	}else{
								VALUETYPE dist = (this->nCoordinates[j] - this->nCoordinates[indices[i]]).getMagnitude();
								if(dist > 0)
                                                        	{
                                                                	f = f - (this->nCoordinates[j] - this->nCoordinates[indices[i]]) * (1.0 / dist);
                                                        	}
                                                	}
                                        	}
					}
					else{
						unordered_map<int, int> neighbors;
                                        	neighbors.insert(pair<int, int>(indices[i], indices[i]));
						if(i < graph.rows - 1){
                        				for(INDEXTYPE j = graph.rowptr[i]; j < graph.rowptr[i+1]; j++){
								f += (nCoordinates[graph.colids[j]] - nCoordinates[indices[i]]) * (W * (nCoordinates[graph.colids[j]] - nCoordinates[indices[i]]).getMagnitude());
								STACKnode.push(graph.colids[j]);
								neighbors.insert(pair<int, int>(graph.colids[j], indices[i]));
                        				}
                				}else{
                        				for(INDEXTYPE j = graph.rowptr[i]; j < graph.nnz; j++){
								f += (nCoordinates[graph.colids[j]] - nCoordinates[indices[i]]) * (W * (nCoordinates[graph.colids[j]] - nCoordinates[indices[i]]).getMagnitude());
                                                        	STACKnode.push(graph.colids[j]);
                                                       		neighbors.insert(pair<int, int>(graph.colids[j], indices[i]));
                        				}
                				}
						int countNodes = 200;
						while(!STACKnode.empty()){
							int currentn = STACKnode.top();
							STACKnode.pop();
							if(currentn < graph.rows - 1){
                                                 	       for(INDEXTYPE n = graph.rowptr[currentn]; n < graph.rowptr[currentn+1]; n++){
                                                        	        if(neighbors.count(graph.colids[n]) < 1){
                                                                	        f += calcRepulsion(indices[i], graph.colids[n]);
                                                                        	STACKnode.push(graph.colids[n]);
                                                                        	countNodes--;
                                                                	}
                                                        	}
                                                	}else{
                                                        	for(INDEXTYPE n = graph.rowptr[currentn]; n < graph.nnz; n++){
                                                                	if(neighbors.count(graph.colids[n]) < 1){
                                                                        	f += calcRepulsion(indices[i], graph.colids[n]);
                                                                        	STACKnode.push(graph.colids[n]);
                                                                        	countNodes--;
                                                                	}
                                                        	}
                                                	}
							if(countNodes <= 0)break;
						}
                                	}
					prevCoordinates[indices[i]] = f;
				}
                                for(INDEXTYPE i = b * BATCHSIZE; i < (b + 1) * BATCHSIZE; i++){
                                        if(i >= graph.rows) continue;
                                        nCoordinates[indices[i]] = nCoordinates[indices[i]] + prevCoordinates[indices[i]].getUnitVector() * STEP;
                                	ENERGY += prevCoordinates[indices[i]].getMagnitude2();
				}
                        }
			STEP = STEP * t;
                        LOOP++;
                }
                end = omp_get_wtime();
		cout << "Greedy Approximation" << endl;
                cout << "Greedy Approximation Energy:" << ENERGY << endl;
                cout << "Greedy Approximation Wall time required:" << end - start << endl;
                result.push_back(ENERGY);
                result.push_back(end - start);
		writeToFile("GAPPROX"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(ITERATIONS));
                return result;
        }	
	
	vector<VALUETYPE> approxCacheBlockBH(INDEXTYPE ITERATIONS, INDEXTYPE NUMOFTHREADS, INDEXTYPE BATCHSIZE){
                INDEXTYPE approxITER = (int)(ITERATIONS * 0.8);
		VALUETYPE start, end;
		vector<VALUETYPE> result;
		start = omp_get_wtime();
		if(init == 0){
                        initDFS(0);
                }else if(init == 2){
                        fileInitialization();
                }
                else{
                        randInit();
                }
		BarnesHutApproximation(approxITER, NUMOFTHREADS, BATCHSIZE, 1.2, 1);
		result = cacheBlockingminiBatchForceDirectedAlgorithm(ITERATIONS-approxITER, NUMOFTHREADS, BATCHSIZE, 6, 0.4, 1, 400);	
		end = omp_get_wtime();
		result[1] = end - start;
		cout << "80%% - 20%% BH - CB" << endl;
		cout << "BH-CACHE Approximation Energy:" << result[0] << endl;
                cout << "BH-CACHE Approximation Wall time required:" << end - start << endl;
		writeToFile("BCAPPROX"+ to_string(BATCHSIZE)+"PARAOUT" + to_string(ITERATIONS));
		return result;
		
	}

	void print(){
		for(INDEXTYPE i = 0; i < graph.rows; i++){
                	cout << "Node:" << i << ", X:" << nCoordinates[i].getX() << ", Y:" << nCoordinates[i].getY()<< endl;
        	}
		cout << endl;
	}
	void writeRepulsiveForce(vector<Coordinate<VALUETYPE> > &repulse, string f){
		ofstream output;
		output.open(f);
		for(INDEXTYPE i = 0; i < graph.rows; i++){
			output << repulse[i].getMagnitude2() << "\t" << repulse[i].getX() << "\t" << repulse[i].getY() << endl;
		}
		output.close();
	}
	void writeToFileBH(Coordinate<VALUETYPE> *tCoordinates, string f){
		stringstream  data(filename);
                string lasttok;
                while(getline(data,lasttok,'/'));
                filename = outputdir + lasttok + f + ".txt";
                ofstream output;
                output.open(filename);
                cout << "Creating output file in following directory:" << filename << endl;
                for(INDEXTYPE i = 0; i < graph.rows; i++){
                        output << tCoordinates[i].getX() <<"\t"<< tCoordinates[i].getY() << "\t" << i+1 << endl;
                }
                output.close();
	}
	void writeToFile(string f){
		stringstream  data(filename);
    		string lasttok;
    		while(getline(data,lasttok,'/'));
		filename = outputdir + lasttok + f + ".txt";
		ofstream output;
		output.open(filename);
		cout << "Creating output file in following directory:" << filename << endl;
		for(INDEXTYPE i = 0; i < graph.rows; i++){
			output << nCoordinates[i].getX() <<"\t"<< nCoordinates[i].getY() << "\t" << i+1 << endl;
		}
		output.close();
	}
};
#endif
