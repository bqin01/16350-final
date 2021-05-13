/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include <utility>
#include <queue>
#include <algorithm>
#include <chrono>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ROBOT_IN	prhs[1]
#define	GOAL_IN     prhs[2]
#define PAYLOAD_IN  prhs[3]
#define FUEL_IN 	prhs[4]
#define C_IN			prhs[5]
#define D_IN			prhs[6]
#define P_IN			prhs[7]
#define F_IN			prhs[8]
#define FCUR_IN		prhs[9]
#define FMAX_IN		prhs[10]

/* Notable types */
#define ROBOT 0
#define TARGET 1
#define PAYLOAD 2
#define FUEL 3

/* Output Arguments */
#define	ACTION_OUT	plhs[0]

//access to the map is shifted to account for 0-based indexing in the map, whereas
//1-based indexing in matlab (so, robotpose and goalpose are 1-indexed)
#define GETMAPINDEX(X, Y, XSIZE, YSIZE) ((Y-1)*XSIZE + (X-1))

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define NUMOFDIRS 8

const double INF = 2000000000;

int temp = 0;
bool initEdges = false;
double epsilon = 4.0;
bool payloaded = false;

int ABS(int x)
{
	return (x<0) ? -x : x;
}

typedef struct notable
{
	int type;
	int x;
	int y;
	notable(int atype, int ax, int ay){
		type = atype;
		x = ax;
		y = ay;
	}
} notable_t;

typedef struct node
{
	int g;
	double prior;
	bool closed;
	int x;
	int y;
	int parent_x;
	int parent_y;
} node_t;

node_t* nodes[5002][5002];


int N_NOTABLES = 2;
notable_t* notables[12];
double notable_adjmat[24][24];
std::pair<double,int> true_adjmat[24][24];

int CURRENT_TARGET = 0;

int halfeuclidH (int x1, int x2, int y1, int y2)
{
	return 0.5*sqrt((x1-x2)*(x1-x2)*1.0 + (y1-y2)*(y1-y2)*1.0);
}

void compPrior(node_t* n, int goalX, int goalY)
{
	n->prior = n->g + epsilon * halfeuclidH(n->x, goalX, n->y, goalY);
}

node_t* nodePair(std::pair<int,int> p)
{
	return nodes[p.first][p.second];
}

void initNode(std::pair<int,int> coords)
{
	int i = coords.first;
	int j = coords.second;
	nodes[i][j] = new node_t();
 	nodes[i][j]->g = INF;
 	nodes[i][j]->prior = INF;
	nodes[i][j]->closed = false;
	nodes[i][j]->x = i;
	nodes[i][j]->y = j;
	nodes[i][j]->parent_x = -1;
	nodes[i][j]->parent_y = -1;
}
void deleteNodes(std::queue< std::pair<int,int> > *toDelete)
{
	while(!((*toDelete).empty()))
	{
		std::pair<int,int> coord = (*toDelete).front();
		delete nodes[coord.first][coord.second];
		(*toDelete).pop();
	}
}
struct CustomCompare {
public:
	bool operator() (const std::pair<int,int>& n1, const std::pair<int,int>& n2)
	{
		return nodePair(n1)->prior > nodePair(n2)->prior;
	}
};

struct dijkstraCompare {
public:
	bool operator() (const std::pair<int,double>& n1, const std::pair<int,double>& n2)
	{
		return n1.second > n2.second;
	}
};

std::pair<int,int> bestAction(int rX, int rY, int gX, int gY)
{
	if(nodes[gX][gY] == NULL) return std::make_pair(0,0);
	else if (nodes[gX][gY]->parent_x == rX && nodes[gX][gY]->parent_y == rY)
	{
		return std::make_pair(gX-rX,gY-rY);
	}else{
		return bestAction(rX,rY,nodes[gX][gY]->parent_x,nodes[gX][gY]->parent_y);
	}
}

std::pair<int,int> A_star_micro(double* map, int x_size, int y_size, int robotposeX, int robotposeY, int goalposeX, int goalposeY, const double* decreasingEpsilons, int nEpsilon, int max_ms)
{
	for(int i = 0; i < 5002; i++)
		for(int j = 0; j < 5002; j++)
			nodes[i][j] = NULL;
	int dX[NUMOFDIRS] = {1, -1, 0,  0, 1,  1, -1, -1};
	int dY[NUMOFDIRS] = {0,  0, 1, -1, 1, -1,  1, -1};
	std::priority_queue< std::pair<int,int>, std::vector< std::pair<int,int> >, CustomCompare> OPEN;
	std::priority_queue< std::pair<int,int>, std::vector< std::pair<int,int> >, CustomCompare> INCONS;
	std::queue< std::pair<int,int> > initializedNodes;

	OPEN.push(std::make_pair(robotposeX, robotposeY));
	initNode(std::make_pair(robotposeX, robotposeY));
	initializedNodes.push(std::make_pair(robotposeX, robotposeY));
	nodes[robotposeX][robotposeY]->g = 0;
	compPrior(nodes[robotposeX][robotposeY],goalposeX,goalposeY);
	int currentIter = 0;
	int timeDifference = 0;

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	while(currentIter <= nEpsilon && timeDifference <= max_ms){
		epsilon = decreasingEpsilons[currentIter];
		while(!INCONS.empty())
		{
			OPEN.push(INCONS.top());
			INCONS.pop();
		}
		double threshold = INF;
		while(!OPEN.empty() && threshold > nodePair(OPEN.top())->prior)
		{
			std::pair<int,int> curS = OPEN.top();
			OPEN.pop();
			nodePair(curS)->closed = true;
			int current_g = nodePair(curS)->g;
			for(int i = 0; i < NUMOFDIRS; i++)
			{
				int newX = curS.first + dX[i];
				int newY = curS.second + dY[i];
				if(newX >= 1 && newX <= x_size && newY >= 1 && newY <= y_size && (int)map[GETMAPINDEX(newX,newY,x_size,y_size)] >= 0){
					if(nodes[newX][newY] == NULL){
						initNode(std::make_pair(newX,newY));
						initializedNodes.push(std::make_pair(newX,newY));
					}
					node_t* nextNode = nodes[newX][newY];
					if(nextNode->g > 1 + current_g)
					{
						nextNode->g = 1 + current_g;
						nextNode->parent_x = curS.first;
						nextNode->parent_y = curS.second;
						compPrior(nextNode,goalposeX,goalposeY);
						if (!(nextNode->closed)){
							OPEN.push(std::make_pair(newX, newY));
						}else{
							INCONS.push(std::make_pair(newX, newY));
						}
					}
				}
			}
			if(nodes[goalposeX][goalposeY] != NULL){
				 threshold = nodes[goalposeX][goalposeY]->prior;
			}
		}
		currentIter++;
		std::chrono::steady_clock::time_point checkpoint = std::chrono::steady_clock::now();
		timeDifference = std::chrono::duration_cast<std::chrono::milliseconds>(checkpoint - begin).count();
	}
	std::pair<int,int> bA = bestAction(robotposeX,robotposeY,goalposeX,goalposeY);
	deleteNodes(&initializedNodes);
	return bA;
}

std::pair< std::vector<double>, std::vector<int> > sssp_macro(int start_i)
{
	std::vector<double> v;
	std::vector<int> parents;
	for(int i = 0; i < 24; i++) v.push_back(INF);
	for(int i = 0; i < 24; i++) parents.push_back(-2);
	bool out[24];
	std::fill(out,out+24,false);
	v[start_i] = 0;
	parents[start_i] = -1;
	int nodes_visited = 0;
	int c_node = start_i;
	double c_max = 0;
	while(nodes_visited < 2*N_NOTABLES && c_max < INF)
	{
		for(int i = 0; i < 2*N_NOTABLES; i++)
		{
			if(i != c_node && !out[i] && notable_adjmat[c_node][i] < INF)
			{
				double obj = v[c_node] + notable_adjmat[c_node][i];
				if(obj < v[i])
				{
					parents[i] = c_node;
					v[i] = obj;
				}
			}
		}
		out[c_node] = true;
		nodes_visited++;
		c_node = -1;
		c_max = INF;
		for(int i = 0; i < 2*N_NOTABLES; i++)
		{
			if(v[i] < c_max && !out[i]){
				c_node = i;
				c_max = v[i];
			}
		}
	}
	return std::make_pair(v, parents);
}

int new_target(int start_i, int fuel)
{
	std::vector<int> fuels;
	int targ;
	for(int i = 0; i < N_NOTABLES; i++)
	{
		if(notables[i]->type == TARGET) targ = i+N_NOTABLES;
		else if(notables[i]->type == FUEL){
			fuels.push_back(i);
			fuels.push_back(i+N_NOTABLES);
		}
	}
	//Estimated enough fuel to get to goal
	if(fuel >= true_adjmat[start_i][targ].first)
	{
		int prev = targ;
		int cur = true_adjmat[start_i][targ].second;
		while(cur != start_i && cur > 0)
		{
			prev = cur;
			cur = true_adjmat[start_i][prev].second;
		}
		return prev;
	}

	//Greedily choose the closest fuel station to the goal that you can reach

	else{
		double cur_min = INF;
		int cur_best_fuel = -1;
		for(int i = 0; i < fuels.size(); i++)
		{
			if(fuel >= true_adjmat[start_i][fuels[i]].first)
			{
				double heur = true_adjmat[fuels[i]][targ].first;
				if(heur < cur_min)
				{
					cur_min = heur;
					cur_best_fuel = fuels[i];
				}
			}
		}
		if(cur_best_fuel < 0)
			return -1;

		int prev = cur_best_fuel;
		int cur = true_adjmat[start_i][cur_best_fuel].second;
		while(cur != start_i && cur > 0)
		{
			prev = cur;
			cur = true_adjmat[start_i][prev].second;
		}
		return prev;
	}
}

static void planner(
		   double* map,
		   int x_size,
 		   int y_size,
           int robotposeX,
            int robotposeY,
            int goalposeX,
            int goalposeY,
						double* payloadPos,
						int payloads,
						double* fuelPos,
						int fuels,
			 int C,
			 int D,
			 int P,
			 int F,
			 int Fcur,
			 int Fmax,
			 char *p_actionX,
			 char *p_actionY
		   )
{
		printf("Entered planner...\n");
    //8-connected grid
		const double decreasingEpsilons[] = {10.0, 8.0, 6.0, 4.0, 3.0, 2.0, 1.5, 1.0};
		const double just20[] = {20.0};

		if(!initEdges)
		{
			printf("Initializing macroedges...\n");
			for(int i = 0; i < 24; i++){
				for(int j = 0; j < 24; j++){
					notable_adjmat[i][j] = INF;
					true_adjmat[i][j] = std::make_pair(INF,-1);
				}
			}

			initEdges = true;

			notables[0] = new notable_t(ROBOT, robotposeX, robotposeY);
			notables[1] = new notable_t(TARGET, goalposeX, goalposeY);
			for(int p = 0; p < payloads; p++)
			{
				notables[2+p] = new notable_t(PAYLOAD, payloadPos[GETMAPINDEX(p,0,payloads,2)], payloadPos[GETMAPINDEX(p,1,payloads,2)]);
				N_NOTABLES++;
			}
			for(int f = 0; f < fuels; f++)
			{
				notables[2+payloads+f] = new notable_t(FUEL, fuelPos[GETMAPINDEX(f,0,fuels,2)], fuelPos[GETMAPINDEX(f,1,fuels,2)]);
				N_NOTABLES++;
			}

			for(int i = 0; i < N_NOTABLES; i++)
			{
				for(int j = i+1; j < N_NOTABLES; j++)
				{
					printf("Initializing macro-edge between %d and %d...\n",i,j);
					int counter = 0;
					int rX = notables[i]->x;
					int rY = notables[i]->y;
					int gX = notables[j]->x;
					int gY = notables[j]->y;
					while(ABS(rX-gX) > 1 || ABS(rY-gY) > 1)
					{
						std::pair<int,int> bA = A_star_micro(map, x_size, y_size, rX, rY, gX, gY, just20, 1, 50);
						rX += bA.first;
						rY += bA.second;
						counter++;
					}
					notable_adjmat[i][j] = counter * (C+D);
					notable_adjmat[j][i] = counter * (C+D);
					notable_adjmat[i+N_NOTABLES][j+N_NOTABLES] = counter * (C+D+P);
					notable_adjmat[j+N_NOTABLES][i+N_NOTABLES] = counter * (C+D+P);
				}
			}
			for(int p = 0; p < payloads; p++)
			{
				notable_adjmat[2+p][2+p+N_NOTABLES] = 0;
			}
			for(int i = 0; i < 2*N_NOTABLES; i++)
			{
				std::pair< std::vector<double>, std::vector<int> > sssp = sssp_macro(i);
				for(int j = 0; j < 2*N_NOTABLES; j++)
				{
					true_adjmat[i][j] = std::make_pair(sssp.first[j], sssp.second[j]);
				}
			}

			CURRENT_TARGET = new_target(0, Fcur);
		}
		if(!payloaded){
			for(int i = 0; i < N_NOTABLES; i++)
			{
				if(notables[i]->type == PAYLOAD && robotposeX == notables[i]->x && robotposeY == notables[i]->y)
				{
					payloaded = true;
				}
				//Currently fueling
				if(notables[i]->type == FUEL && robotposeX == notables[i]->x && robotposeY == notables[i]->y && ABS(Fcur-Fmax) > 1)
				{
					*p_actionX = 0;
					*p_actionY = 0;
			    return;
				}
			}
		}

		//Current target reached, time to set a new one!
		if(notables[CURRENT_TARGET]->x == robotposeX && notables[CURRENT_TARGET]->y == robotposeY && payloaded == (CURRENT_TARGET >= N_NOTABLES ? true : false))
		{
			CURRENT_TARGET = new_target(CURRENT_TARGET, Fcur);
		}

		std::pair<int,int> bA = A_star_micro(map, x_size, y_size, robotposeX, robotposeY, notables[CURRENT_TARGET]->x, notables[CURRENT_TARGET]->y, decreasingEpsilons, 8, 200);

		*p_actionX = bA.first;
		*p_actionY = bA.second;
    return;
}


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )

{

		printf("Entering MEX\n");
    /* Check for proper number of arguments */
    if (nrhs != 11) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Eleven input arguments required.");
    } else if (nlhs != 1) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required.");
    }

    /* get the dimensions of the map and the map matrix itself*/

		printf("Checkpoint\n");
    int x_size = mxGetM(MAP_IN);
    int y_size = mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);

    /* get the dimensions of the robotpose and the robotpose itself*/
				printf("Checkpoint\n");
    int robotpose_M = mxGetM(ROBOT_IN);
    int robotpose_N = mxGetN(ROBOT_IN);
    if(robotpose_M != 1 || robotpose_N != 2){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidrobotpose",
                "robotpose vector should be 1 by 2.");
    }
    double* robotposeV = mxGetPr(ROBOT_IN);
    int robotposeX = (int)robotposeV[0];
    int robotposeY = (int)robotposeV[1];

    /* get the dimensions of the goalpose and the goalpose itself*/
				printf("Checkpoint\n");
    int goalpose_M = mxGetM(GOAL_IN);
    int goalpose_N = mxGetN(GOAL_IN);
    if(goalpose_M != 1 || goalpose_N != 2){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidgoalpose",
                "goalpose vector should be 1 by 2.");
    }
				printf("Checkpoint\n");
    double* goalposeV = mxGetPr(GOAL_IN);
    int goalposeX = (int)goalposeV[0];
    int goalposeY = (int)goalposeV[1];

				printf("Checkpoint\n");
		double* payloadpos = mxGetPr(PAYLOAD_IN);
		int payloadpose_N = mxGetN(PAYLOAD_IN);

				printf("Checkpoint\n");
		double* fuelpos = mxGetPr(FUEL_IN);
		int fuelpose_N = mxGetN(FUEL_IN);
		/*
			#define C_IN			prhs[5]
			#define D_IN			prhs[6]
			#define P_IN			prhs[7]
			#define F_IN			prhs[8]
			#define FCUR_IN		prhs[9]
			#define FMAX_IN		phrs[10]
		*/
				printf("Checkpoint\n");
		int C = (int)(mxGetPr(C_IN)[0]);
		int D = (int)(mxGetPr(D_IN)[0]);
		int P = (int)(mxGetPr(P_IN)[0]);
		int F = (int)(mxGetPr(F_IN)[0]);
		int Fcur = (int)(mxGetPr(FCUR_IN)[0]);
		int Fmax = (int)(mxGetPr(FMAX_IN)[0]);

    /* Create a matrix for the return action */
				printf("Checkpoint\n");
    ACTION_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)2, mxINT8_CLASS, mxREAL);
    char* action_ptr = (char*)  mxGetPr(ACTION_OUT);

    /* Do the actual planning in a subroutine */
				printf("Checkpoint\n");
    planner(map, x_size, y_size, robotposeX, robotposeY, goalposeX, goalposeY, payloadpos, payloadpose_N, fuelpos, fuelpose_N, C, D, P, F, Fcur, Fmax, &action_ptr[0], &action_ptr[1]);
    return;

}
