#include <ilcplex/ilocplex.h>//cplex
#include <iostream>
#include <stdio.h>
#include <cstdlib>
/*MIT License

Copyright(c) 2026 Wei Yang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <time.h>
#include <math.h>
#include <cstring>
#include <algorithm>
#include <vector>
#include <fstream>
#include <unistd.h>
#define MAX 1000
#define MAXPOOL 10
#define MAXdoor 40

using namespace std;
typedef long long ll;
vector<vector<ll>> c;
vector<int> mate_V, mate_U, nhbor, parent;
vector<ll> alpha, beta, slack, min_col;
vector<bool> label_V, label_U;
struct partion
{
	int num_truck;
	int truck[200];
};
partion * X, *Y;

struct SP
{
	int  obj;
	double time;
};
std::vector<SP> SearchProgress;

int M, N, I, J;
long long minf = 9999999999;
long long now_f = 0;
long long best_minf = 9999999999;
int maxiter = 2000000000;
long long objective;

int Flow[MAX][MAX];
int Distance[MAX][MAX];
int C[MAX][MAX];
int Rs[MAX];
int Rr[MAX];
int Ss[MAX];
int Sr[MAX];
int Tot_Ss[MAX];
int Tot_Sr[MAX];
int s[MAX];
int t[MAX];
int best_s[MAX], best_t[MAX], temp_best_s[MAX], temp_best_t[MAX], temp_s[MAX], temp_t[MAX];
int new_flow[MAX][MAX] = { 0 };
int cost_s[MAX][MAX], cost_t[MAX][MAX];
bool shift[MAX] = { 0 };
int oversize_s[MAX], oversize_t[MAX], Tot_OverS, Tot_OverT, previous_Tot_OS, previous_Tot_OT, PTot_OS, PTot_OT;
int ** last_swaped_s, **last_swaped_t, **last_swaped_s2, **last_swaped_t2;
int * determine;
int iteration;
double op, ed, t_f_best = 0.0;//op is start time, ed is end time, t_f_best is the time when best solution was finded;
double r1 = 0.1, r2 = 0.2, r3 = 0.3;//pertrubation strength
double p1 = 1.1;//Tabu tenure p1*M
int TTP = 5; //Threshold for distinguishing Case3 from Case4
long t_all = 1200;//second
double Q1, Q2;
long long ReStartTimes = 10;
int depth_LS = 1;
int depthall = 500; //Maximum iteration in every LA-ITS (LandScapeDetectionTabuSearch())
int pertubation_l, flag;
int gap_iter = 2000;
int weak_depth = 20000, medium_depth = 40000, strong_depth = 100000;
double tollerance = 0.0001;
double relative_gap = 0.0001;
int count_gap = 0;
int success = 0;
int pthreads = 1;
int num_obj = 1;
char   resultfilename[50] = "result/CDAP/";
//char   poolfilename[50] = "pool/";
char   filename[30] = "tuning/Instances/CDAP/";
char * resultpath = "SetB_100x10S15.txt";


int Permutation_local_search_start_from_s_with_depth(int _row, int _col, long long & _current_value, long long & _minf, long long & _best_value, int * _s, int * _t, int * temp_best_s, int * temp_best_t,
	int(&matrix_a)[1000][1000], int(&matrix_b)[1000][1000], int(&g)[1000][1000], double & ed, double t_all, int _depth);

void pertubation(int pertubation_l, int* s, int *t, int M, int I, int N, int J);
void determine_pertubation_length2(int & pertubation_l, int &flag, int & max_gap,
	long long previous_value,//上max_gap次的结果
	long long current_value,//目前的结果
	long long previous_best,//上一次扰动之前搜索到的结果
	long long best_value); //这一次扰动之前搜索到的结果
void determine_pertubation_length3(int & pertubation_l, int &flag, int & max_gap,
	long long previous_value,//上max_gap次的结果
	long long current_value,//目前的结果
	long long previous_best,//上一次扰动之前搜索到的结果
	long long best_value); //这一次扰动之前搜索到的结果

/*Huangarin algorithm*/
bool unlabelled_U(int u) { return  !label_U[u]; }
bool unmatched_V(int v) { return mate_V[v] == -1; }
bool unmatched_U(int u) { return mate_U[u] == -1; }
bool admissible_U(int u) { return slack[u] == 0LL; }


void augment(int v, int exposed_u) {

	int aux = mate_V[v];

	mate_V[v] = exposed_u;
	mate_U[exposed_u] = v;

	if (parent[v] != -1)
		augment(parent[v], aux);

}

void update_slack(int v, int length) {

	// for unlabelled u in U
	int N = length;
	for (int u = 0; u < N; u++)
		if (unlabelled_U(u)) {
			ll bound = c[v][u] - alpha[v] - beta[u];
			if (0LL <= bound && bound < slack[u]) {
				slack[u] = bound;
				nhbor[u] = v;
			}
		}
}

ll update_alpha_beta(int length) {
	int N = length;
	ll theta = numeric_limits<ll>::max();

	// for unlabelled u in U
	for (int u = 0; u < N; u++)
		if (unlabelled_U(u))
			theta = min(theta, slack[u]);

	// skip if theta == 0 (no update needed)
	if (theta > 0LL) {

		// integrality is ensured
		theta /= 2LL;

		for (int i = 0; i < N; i++) {
			if (label_V[i]) alpha[i] += theta;
			else alpha[i] -= theta;
			if (label_U[i]) beta[i] -= theta;
			else beta[i] += theta;
		}
	}

	return theta;
}

int search_augmenting_alternating_path(int length) {
	int N = length;
	while (true) {
		ll theta = update_alpha_beta(N);

		vector<int> admissibles = vector<int>();

		// for unlabelled u in U
		for (int u = 0; u < N; u++)
			if (unlabelled_U(u)) {
				slack[u] -= 2LL * theta;
				if (admissible_U(u)) {
					// unlabelled, admissible and unmatched => path found
					if (unmatched_U(u)) return u;
					else admissibles.push_back(u);
				}
			}

		for (int u : admissibles) {
			label_U[u] = true;
			label_V[mate_U[u]] = true;
			parent[mate_U[u]] = nhbor[u];
			update_slack(mate_U[u], N);
		}
	}

}

void initialize_search(int length) {
	int N = length;
	nhbor = vector<int>(N, -1);
	parent = vector<int>(N, -1);
	slack = vector<ll>(N, numeric_limits<ll>::max());
	label_V = vector<bool>(N, false);
	label_U = vector<bool>(N, false);

}

void initialize_alpha_beta(int length) {
	int N = length;
	mate_V = vector<int>(N, -1);
	mate_U = vector<int>(N, -1);
	alpha = vector<ll>(N, 0LL);
	beta = vector<ll>(N);
	for (int u = 0; u < N; u++)
		beta[u] = min_col[u];

}
void Hungarian_input(int(&Prob)[MAX][MAX], int length) {
	int N = length;
	c = vector<vector<ll>>(N, vector<ll>(N));
	min_col = vector<ll>(N, numeric_limits<ll>::max());
	for (int v = 0; v < N; v++)
		for (int u = 0; u < N; u++) {
			c[v][u] = Prob[v][u];

			// multiply by 2LL to ensure integrality
			c[v][u] *= 2LL;
			min_col[u] = min(min_col[u], c[v][u]);
		}
}
long long Hungarian_algorithm(int(&Prob)[MAX][MAX], int length, int *permutation) {
	Hungarian_input(Prob, length);
	int N = length;
	initialize_alpha_beta(N);

	for (int i = 0; i < N; i++) {
		initialize_search(N);

		// start with unmatched v in V
		for (int v = 0; v < N; v++)
			if (unmatched_V(v)) {
				label_V[v] = true;
				update_slack(v, N);
			}

		int u = search_augmenting_alternating_path(N);

		augment(nhbor[u], u);
	}
	for (int v = 0; v < N; v++)
		permutation[v] = mate_V[v];
	long long opt_cost = 0LL;
	for (int i = 0; i < N; i++)
		opt_cost += alpha[i] + beta[i];
	return opt_cost / 2LL;
}



void allocate_memory()
{
	//last_swaped = new int*[M];
	determine = new int[MAXPOOL];
	last_swaped_s = new int*[M + 1];
	last_swaped_t = new int*[N + 1];
	last_swaped_s2 = new int*[M + 1];
	last_swaped_t2 = new int*[N + 1];
	X = new partion[I];
	Y = new partion[J];
	for (int i = 0; i < M + 1; i++) {
		last_swaped_s[i] = new int[M + 1];
		last_swaped_s2[i] = new int[M + 1];
	}
	for (int i = 0; i < N + 1; i++) {
		last_swaped_t[i] = new int[N + 1];
		last_swaped_t2[i] = new int[N + 1];
	}
	for (int i = 0; i < M + 1; i++) {
		for (int j = 0; j < M + 1; j++) {
			last_swaped_s[i][j] = 0;
			last_swaped_s2[i][j] = 0;
		}
	}
	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < N + 1; j++) {
			last_swaped_t[i][j] = 0;
			last_swaped_t2[i][j] = 0;
		}
	}
	for (int i = 0; i < I; i++) {
		X[i].num_truck = 0;
	}
	for (int i = 0; i < J; i++) {
		Y[i].num_truck = 0;
	}

}
void delete_memory()
{
	delete[] X;
	delete[] Y;
	for (int i = 0; i < M; i++) {
		delete[] last_swaped_s[i];
	}
	for (int i = 0; i < N; i++) {
		delete[] last_swaped_t[i];
	}
	delete[] last_swaped_s;
	delete[] last_swaped_t;
}
int init()
{
	minf = 9999999999;

	for (int i = 0; i < MAX; i++)
	{
		s[i] = 0;//row
		t[i] = 0;//col

	}
	for (int i = 0; i < M + 1; i++) {
		for (int j = 0; j < M + 1; j++) {
			last_swaped_s[i][j] = 0;
			last_swaped_s2[i][j] = 0;
		}
	}
	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < N + 1; j++) {
			last_swaped_t[i][j] = 0;
			last_swaped_t2[i][j] = 0;
		}
	}
	return 0;
}
void ReadInstance(char * filename) {
	ifstream readinstance(filename, ios::in);
	readinstance >> M >> I >> J >> N;
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			readinstance >> Flow[i][j];
		}
	}
	for (int i = 0; i < I; i++) {
		for (int j = 0; j < J; j++) {
			readinstance >> Distance[i][j];
		}

	}
	for (int i = 0; i < M; i++) {
		readinstance >> Ss[i];
	}
	Ss[M] = 0;
	for (int i = 0; i < I; i++) {
		readinstance >> Rs[i];
	}
	for (int i = 0; i < N; i++) {
		readinstance >> Sr[i];
	}
	Sr[N] = 0;
	for (int i = 0; i < J; i++) {
		readinstance >> Rr[i];
	}
	for (int i = 0; i < M; i++) {
		readinstance >> temp_best_t[i];
	}
}

long long compute_value(int(&Flow)[1000][1000], int(&Distance)[1000][1000], int *s, int *t, int M, int N, int I, int J) {
	long long value = 0;
	int ** new_Flow;
	new_Flow = new int*[I];
	for (int i = 0; i < I; i++) {
		new_Flow[i] = new int[J];
	}
	for (int i = 0; i < I; i++) {
		for (int j = 0; j < J; j++) {
			new_Flow[i][j] = 0;
		}
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			new_Flow[s[i]][t[j]] += Flow[i][j];
		}
	}
	for (int i = 0; i < I; i++) {
		for (int j = 0; j < J; j++) {
			value += new_Flow[i][j] * Distance[i][j];
		}
	}
	for (int i = 0; i < I; i++) {
		delete[] new_Flow[i];
	}
	delete[] new_Flow;
	return value;
}
long long compute_value2(int(&Flow)[1000][1000], int(&Distance)[1000][1000], int *s, int *t, int M, int N, int I, int J) {
	long long value = 0;
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			value += Flow[i][j] * Distance[s[i]][t[j]];
		}
	}
	return value;
}

int compute_delta_s(int i1, int i2, int m, int *_s, int *_t, int(&Flow)[1000][1000], int(&Distance)[1000][1000], int N) {
	int delta = 0;
	for (int n = 0; n < N; n++) {
		delta += Distance[i2][_t[n]] * Flow[m][n] - Distance[i1][_t[n]] * Flow[m][n];
	}
	return delta;
}
int compute_delta_t(int j1, int j2, int n, int *_s, int *_t, int(&Flow)[1000][1000], int(&Distance)[1000][1000], int M) {
	int delta = 0;
	for (int m = 0; m < M; m++) {
		delta += Distance[_s[m]][j2] * Flow[m][n] - Distance[_s[m]][j1] * Flow[m][n];
	}
	return delta;
}
int compute_cost_assignmet_one_intruck(int intruck, int indoor, int *_t, int N, int(&Flow)[MAX][MAX], int(&Distance)[MAX][MAX]) {
	int cost = 0;
	for (int i = 0; i < N; i++) {
		cost += Distance[indoor][_t[i]] * Flow[intruck][i];
	}
	return cost;
}
int compute_cost_assignmet_one_outtruck(int outtruck, int outdoor, int *_s, int M, int(&Flow)[MAX][MAX], int(&Distance)[MAX][MAX]) {
	int cost = 0;
	for (int i = 0; i < M; i++) {
		cost += Distance[_s[i]][outdoor] * Flow[i][outtruck];
	}
	return cost;
}
void initial_cost(int *_s, int *_t, int M, int I, int N, int J, int(&Flow)[MAX][MAX], int(&Distance)[MAX][MAX], int(&cost_s)[MAX][MAX], int(&cost_t)[MAX][MAX]) {
	for (int i = 0; i <= M; i++) {
		for (int j = 0; j <= I; j++) {
			cost_s[i][j] = 0;
		}
	}
	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			cost_t[i][j] = 0;
		}
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < I; j++) {
			cost_s[i][j] = compute_cost_assignmet_one_intruck(i, j, _t, N, Flow, Distance);
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < J; j++) {
			cost_t[i][j] = compute_cost_assignmet_one_outtruck(i, j, _s, M, Flow, Distance);
		}
	}
}
void update_S_after_change_T(int M, int I, int(&cost_s)[MAX][MAX], int outtruck, int origindoor, int newdoor) {
	for (int m = 0; m < M; m++) {
		for (int i = 0; i < I; i++) {
			cost_s[m][i] += Flow[m][outtruck] * (Distance[i][newdoor] - Distance[i][origindoor]);
		}
	}
}
void update_T_after_change_S(int N, int J, int(&cost_t)[MAX][MAX], int intruck, int origindoor, int newdoor) {
	for (int n = 0; n < N; n++) {
		for (int j = 0; j < J; j++) {
			cost_t[n][j] += Flow[intruck][n] * (Distance[newdoor][j] - Distance[origindoor][j]);
		}
	}
}
void compute_Tot_capacity(int * _vec, int len_m, int len_d, int*Tot_Ss, int *Ss) {
	for (int i = 0; i < len_d; i++) {
		Tot_Ss[i] = 0;
	}
	for (int i = 0; i < len_m; i++) {
		Tot_Ss[_vec[i]] += Ss[i];
	}
}

bool test_solution(int(&Flow)[1000][1000], int(&Distance)[1000][1000], int *s, int *t, int *Rs, int *Ss, int* Rr, int *Sr, int M, int N, int I, int J) {
	int ** new_Flow;
	int tot = 0;
	new_Flow = new int*[I];
	for (int i = 0; i < I; i++) {
		new_Flow[i] = new int[J];
	}
	for (int i = 0; i < I; i++) {
		for (int j = 0; j < J; j++) {
			new_Flow[i][j] = 0;
		}
	}

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			new_Flow[s[i]][t[j]] += Flow[i][j];
		}
	}
	for (int i = 0; i < I; i++) {
		tot = 0;
		for (int j = 0; j < J; j++) {
			tot += new_Flow[i][j];
		}
		if (tot > Rs[i]) {
			cout << "infeasible" << endl;
			system("pause");
			return 0;
		}
	}
	for (int i = 0; i < J; i++) {
		tot = 0;
		for (int j = 0; j < I; j++) {
			tot += new_Flow[j][i];
		}
		if (tot > Rr[i]) {
			cout << "infeasible" << endl;
			system("pause");
			return 0;
		}
	}
	return 1;
}
void display_result(FILE *fp)
{
	fprintf(fp, "O:%lld and computed:%lld  and if feasible=%d\n", best_minf, compute_value(Flow, Distance, best_s, best_t, M, N, I, J), test_solution(Flow, Distance, best_s, best_t, Rs, Ss, Rr, Sr, M, N, I, J));
	fprintf(fp, "T:%lfs\n", t_f_best);
	fprintf(fp, "s:");
	for (int i = 0; i < M; i++)
		fprintf(fp, "%d ", best_s[i]);
	fprintf(fp, "\n");
	fprintf(fp, "t:");
	for (int i = 0; i < N; i++)
		fprintf(fp, "%d ", best_t[i]);
	fprintf(fp, "\n");
	fprintf(fp, "***********\n");
}
void generate_random_solution(int *sol, int M, int I) {
	for (int m = 0; m < M; m++) {
		sol[m] = rand() % I;
	}
}
void generate_random_solution2(int n, int *permutation)
{
	int i;
	for (i = 0; i < n; i = i + 1) permutation[i] = i;
	for (i = 0; i < n - 1; i = i + 1) swap(permutation[i], permutation[i + rand() % (n - i)]);

}
void generate_random_n_one(bool * zeroone, int n, int num)
{
	for (int i = 0; i < n; i++) {
		zeroone[i] = 0;
	}
	for (int i = 0; i < num; i++)
	{
		int rand_num = rand() % n;
		while (zeroone[rand_num] == 1) {
			rand_num = rand() % n;

		}
		zeroone[rand_num] = 1;
	}
}

int compute_solution_delta_s_and_t(int M, int N, int * _s, int * _t, int * guide_s, int * guide_t)
{
	int delta_s_and_t = 0;
	for (int i = 0; i < M; i++)
	{
		if (_s[i] != guide_s[i])
			delta_s_and_t++;
	}
	for (int i = 0; i < N; i++) {
		if (_t[i] != guide_t[i])
		{
			delta_s_and_t++;
		}
	}
	return delta_s_and_t;
}
void Vector_to_Partion(partion *_X, int * _s, int M, int I) {
	int count[100] = { 0 };
	for (int i = 0; i < M; i++) {
		_X[_s[i]].truck[count[_s[i]]] = i;
		count[_s[i]]++;
	}

	for (int i = 0; i < I; i++) {
		_X[i].num_truck = count[i];
		_X[i].truck[count[i]] = M;
	}
}
void Permutation_To_Vector(int * _p_s, int * _p_t, int * _s, int * _t, int M, int N) {
	for (int i = 0; i < M; i++) {
		_s[i] = _p_s[_s[i]];
	}
	for (int i = 0; i < N; i++) {
		_t[i] = _p_t[_t[i]];
	}
}
void Permutation_To_Vector2(int * _p_s, int * _p_t, int * _s, int * _t, int M, int N) {
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < I; j++) {
			if (_s[i] == _p_s[j]) {
				_s[i] = j;
				break;
			}
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < J; j++) {
			if (_t[i] == _p_t[j]) {
				_t[i] = j;
				break;
			}
		}
	}
}
void Permutation_To_Vector3(int * _p_s, int * _p_t, int * _s, int * _t, int M, int N) {
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < I; j++) {
			if (_s[i] == j) {
				_s[i] = _p_s[j];
				break;
			}
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < J; j++) {
			if (_t[i] == j) {
				_t[i] = _p_t[j];
				break;
			}
		}
	}
}
void generate_new_flow(int *_s, int * _t, int M, int N) {
	for (int i = 0; i < I; i++) {
		for (int j = 0; j < J; j++) {
			new_flow[i][j] = 0;
		}
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			new_flow[_s[i]][_t[j]] += Flow[i][j];
		}
	}
}
void cout_solution(int *s, int *t) {
	for (int i = 0; i < M; i++) {
		cout << s[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < N; i++) {
		cout << t[i] << " ";
	}
	cout << endl;
}

void rand_swap_Vector(int * _vector, partion * _X, int len_m, int len_d, int * _Tot_Ss, int *_Ss, int *_Rs) {
	int indoor1, indoor2, truck1, truck2;
	bool flag = 0;
	for (int i = 0; i < len_d; i++) {
		_Tot_Ss[i] = 0;
	}
	for (int i = 0; i < len_m; i++) {
		_Tot_Ss[_vector[i]] += _Ss[i];
	}
	while (!flag) {
		indoor1 = rand() % len_d;
		indoor2 = rand() % len_d;
		while (indoor2 == indoor1) {
			indoor2 = rand() % len_d;
		}
		truck1 = _X[indoor1].truck[rand() % _X[indoor1].num_truck];
		for (int i = 0; i < _X[indoor2].num_truck; i++) {
			truck2 = _X[indoor2].truck[i];
			if (_Tot_Ss[indoor1] + _Ss[truck2] - _Ss[truck1] <= _Rs[indoor1] && _Tot_Ss[indoor2] + _Ss[truck1] - _Ss[truck2] <= _Rs[indoor2]) {
				flag = 1;
				break;
			}
		}
		if (!flag) {
			truck1 = _X[indoor1].truck[rand() % _X[indoor1].num_truck];
			for (int i = 0; i < _X[indoor2].num_truck; i++) {
				truck2 = _X[indoor2].truck[i];
				if (_Tot_Ss[indoor1] + _Ss[truck2] - _Ss[truck1] <= _Rs[indoor1] && _Tot_Ss[indoor2] + _Ss[truck1] - _Ss[truck2] <= _Rs[indoor2]) {
					flag = 1;
					break;
				}
			}
		}
	}
	_vector[truck1] = indoor2;
	_vector[truck2] = indoor1;
}

void rand_swap_door(int* temp_s, int *temp_t, int M, int I, int N, int J) {
	int rand_door1, rand_door2;
	rand_door1 = rand() % I;
	rand_door2 = rand() % I;
	while (rand_door1 == rand_door2) {
		rand_door2 = rand() % I;
	}
	for (int i = 0; i < M; i++) {
		if (temp_s[i] == rand_door1) {
			temp_s[i] = rand_door2;
		}
		else if (temp_s[i] == rand_door2) {
			temp_s[i] = rand_door1;
		}
	}
	rand_door1 = rand() % J;
	rand_door2 = rand() % J;
	while (rand_door1 == rand_door2) {
		rand_door2 = rand() % J;
	}
	for (int i = 0; i < M; i++) {
		if (temp_t[i] == rand_door1) {
			temp_t[i] = rand_door2;
		}
		else if (temp_t[i] == rand_door2) {
			temp_t[i] = rand_door1;
		}
	}
}
void rand_swap_door_on_side(int * t, int N, int J) {
	int rand_door1, rand_door2;
	rand_door1 = rand() % J;
	rand_door2 = rand() % J;
	while (rand_door1 == rand_door2) {
		rand_door2 = rand() % J;
	}
	for (int i = 0; i < N; i++) {
		if (t[i] == rand_door1) {
			t[i] = rand_door2;
		}
		else if (t[i] == rand_door2) {
			t[i] = rand_door1;
		}
	}
}

long long generate_initial_solution(int *_Rs, int *_Ss, int * sol, int row_num, int I, int t)
{
	long long objectiveValue;
	bool * xx;
	xx = new bool[I];
	generate_random_n_one(xx, I, num_obj);
	IloEnv env;
	try {
		IloModel model(env);

		/*begin of instantiating objective variables where
		variables are presented in two dimensional array
		and each variable such Xij is {0,1}*/
		IloArray<IloBoolVarArray> x(env, row_num);
		for (int i = 0; i < row_num; i++)
		{
			x[i] = IloBoolVarArray(env, I);
		}

		for (int i = 0; i < row_num; i++)
		{
			for (int j = 0; j < I; j++)
			{
				x[i][j] = IloBoolVar(env);
			}
		}//end of Xij
		IloNumVarArray T(env, I);
		for (int i = 0; i < I; i++) {
			T[i] = IloNumVar(env);
		}
		//-----------------------end of instantiating objective variables-------------------------------//

		int m = ceil(I / 2.0);
		//--------------------------begin of instantiating objective function---------------------------//
		IloNumExpr linear_sum(env);

		for (int i = 0; i < I; i++)
		{

			linear_sum += xx[i] * T[i];


		}
		//linear_sum = T[m];
		model.add(IloMaximize(env, linear_sum));
		//linear_sum.end();
		//------------------------------end of objective function----------------------------------------//

		/*each origin/destination must be allocated to one and only one inbound/outbound door*/
		for (int i = 0; i < row_num; i++)
		{
			IloNumExpr constraint1(env);
			for (int j = 0; j < I; j++)
			{
				constraint1 += x[i][j];
			}
			model.add(constraint1 == 1);
			constraint1.end();
		}

		for (int i = 0; i < I; i++)
		{
			IloNumExpr constraint2(env);
			for (int j = 0; j < row_num; j++)
			{
				constraint2 += x[j][i] * _Ss[j];
			}
			model.add(_Rs[i] >= constraint2 + T[i]);
			model.add(T[i] >= 0);
			constraint2.end();
		}
		//-------------------------------------------end of constraints of variable Xij--------------------------//

		//------------------------------------------------end of contriants-------------------------------------//

		//begin of solving 
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, fabs(t));
		cplex.setParam(IloCplex::Param::Threads, pthreads);
		cplex.setParam(IloCplex::Param::MIP::Limits::Solutions, 1);
		//clock_t begin = clock();
		cplex.solve();
		//clock_t end = clock();
		if (cplex.getStatus() == !IloAlgorithm::Optimal)
		{
			cout << "cplex didn't get optimal result" << endl;
		}
		objectiveValue = cplex.getObjValue();
		IloNumArray solution(env);
		for (int i = 0; i < row_num; i++)
		{
			cplex.getValues(solution, x[i]);
			for (int j = 0; j < I; j++)
			{
				if (fabs(solution[j] - 1) < 1e-4)
				{
					sol[i] = j;
					break;
				}
			}
		}

		//de-allocate heap memory
		for (int i = 0; i < row_num; i++)
		{
			x[i].end();
		}
		x.end();
		//cplex.end();
	}
	catch (IloException& e)
	{
		cout << endl;
		cerr << "concert exception caught: " << e << endl;
	}
	catch (...)
	{
		cout << endl;
		cerr << "Unknown exception caught" << endl;
	}

	env.end();
	return objectiveValue;
}

void determine_pertubation_length2(int & pertubation_l, int &flag, int & max_gap,
	long long previous_value,//上max_gap次的结果
	long long  current_value,//目前的结果
	long long previous_best,//上一次扰动之前搜索到的结果
	long long  best_value) //这一次扰动之前搜索到的结果
{
	double gap, best_gap;
	//bool flag = 0;
	gap = (double)(previous_value - current_value) / previous_value;

	if (gap < relative_gap) {
		best_gap = (double)(previous_best - best_value) / previous_best;
		if ((best_gap < tollerance) && (flag >= 4))
		{
			pertubation_l = 1;
			flag++;
			max_gap = strong_depth;
		}
		else if ((best_gap < tollerance) && (flag < 4)) {
			//pertubation_l= rand() % 3 + 1;
			pertubation_l = 1;
			flag++;
			max_gap = medium_depth;
		}
		else {
			//pertubation_l = rand() % 3 + 1;
			pertubation_l = 1;
			flag = 0;
			max_gap = weak_depth;
		}
	}
	else {
		pertubation_l = 0;
		//flag = 0;
	}
	//cout << " relative:" << relative_gap << " tollerance=" << tollerance << " pertubation=" << pertubation_l << " flag=" << flag << endl;
}

void determine_pertubation_length3(int & pertubation_l, int &flag, int & max_gap,
	long long previous_value,//上max_gap次的结果
	long long  current_value,//目前的结果
	long long previous_best,//上一次扰动之前搜索到的结果
	long long  best_value) //这一次扰动之前搜索到的结果
{
	double gap, best_gap;
	//bool flag = 0;
	gap = (double)(previous_value - current_value) / previous_value;

	if (gap < relative_gap) {
		best_gap = (double)(previous_best - best_value) / previous_best;
		if ((best_gap < tollerance) && (flag > TTP))
		{
			pertubation_l = std::max(1, (int)(r3*I));//
			flag++;
			max_gap = strong_depth;
		}
		else if ((best_gap < tollerance) && (flag <= TTP)) {
			//pertubation_l= rand() % 3 + 1;
			pertubation_l = std::max(2, (int)(r2*I));
			flag++;
			max_gap = medium_depth;
		}
		else {
			//pertubation_l = rand() % 3 + 1;
			pertubation_l = std::max(3, (int)(r1*I));
			flag = 0;
			max_gap = weak_depth;
		}
	}
	else {
		pertubation_l = 0;
		//flag = 0;
	}
	//cout << " relative:" << relative_gap << " tollerance=" << tollerance << " pertubation=" << pertubation_l << " flag=" << flag << endl;
}

void compute_over_capacity(int *Tot_Ss, int *_Rs, int *Tot_Sr, int *_Rr) {
	Tot_OverS = 0;
	for (int i = 0; i < I; i++) {

		oversize_s[i] = max(0, Tot_Ss[i] - _Rs[i]);
		Tot_OverS += oversize_s[i];
	}
	Tot_OverT = 0;
	for (int i = 0; i < J; i++) {
		oversize_t[i] = max(0, Tot_Sr[i] - _Rr[i]);
		Tot_OverT += oversize_t[i];
	}
}
bool check_Constraints() {// return 1, the solusion is infeasible.
	for (int i = 0; i < I; i++)
		if (oversize_s[i] > 0) { return 1; }
	for (int i = 0; i < J; i++)
		if (oversize_t[i] > 0) { return 1; }
	return 0;
}

bool check_Constraints(int *oversize, int _I) {// return 1, the solusion is infeasible.
	for (int i = 0; i < _I; i++)
		if (oversize[i] > 0) { return 1; }

	return 0;
}
void pertubation(int pertubation_l, int* s, int *t, int M, int I, int N, int J) {
	for (int i = 0; i < pertubation_l; i++) {
		rand_swap_door(s, t, M, I, N, J);
	}
}


void LandScapeDetectionTabuSearch(int *_vector_s, int *_vector_t, long long& _current_value, long long& _current_best, long long& _previous_best, long long& previous_pertubation_best, long long&_best_value, int *_best_s, int *_best_t, partion * _X, partion * _Y, int M, int I, int N, int J,
	int * _Tot_Ss, int *_Ss, int *_Rs, int * _Tot_Sr, int *_Sr, int *_Rr, int(&cost_s)[MAX][MAX], int(&cost_t)[MAX][MAX],
	int ** last_swaped_s, int ** last_swaped_t, int & iteration, int depth, double &op, double &ed,  int &previous_iter, long t_all) {
	long long _current_minf;
	_current_minf = _current_value;
	int indoor1, indoor2, truck1, truck2, select1, select2, temp_delta, min_delta;
	bool flagY = false;
	Vector_to_Partion(_X, _vector_s, M, I);
	Vector_to_Partion(_Y, _vector_t, N, J);
	for (int iter = iteration; iter < (iteration + depth); iter++) {
		min_delta = 999999999;
		select1 = -1;
		select2 = -1;
		flagY = false;
		for (int i = 0; i < I; i++) {
			for (int j = 0; j <= _X[i].num_truck; j++) {
				truck1 = _X[i].truck[j];
				for (int k = i + 1; k < I; k++) {
					for (int l = 0; l <= _X[k].num_truck; l++) {
						truck2 = _X[k].truck[l];
						if ((truck1 != truck2) && (_Tot_Ss[i] + _Ss[truck2] - _Ss[truck1] <= _Rs[i]) && (_Tot_Ss[k] + _Ss[truck1] - _Ss[truck2] <= _Rs[k])) {
							temp_delta = cost_s[truck1][k] - cost_s[truck1][i] + cost_s[truck2][i] - cost_s[truck2][k];
							if ((temp_delta < min_delta) && ((last_swaped_s[truck1][truck2] + p1 * M < iter) || (temp_delta + _current_value < _current_minf))) {
								min_delta = temp_delta;
								indoor1 = i;
								indoor2 = k;
								select1 = j;
								select2 = l;
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < J; i++) {
			for (int j = 0; j <= _Y[i].num_truck; j++) {
				truck1 = _Y[i].truck[j];
				for (int k = i + 1; k < J; k++) {
					for (int l = 0; l <= _Y[k].num_truck; l++) {
						truck2 = _Y[k].truck[l];
						if ((truck1 != truck2) && (_Tot_Sr[i] + _Sr[truck2] - _Sr[truck1] <= _Rr[i]) && (_Tot_Sr[k] + _Sr[truck1] - _Sr[truck2] <= _Rr[k])) {
							temp_delta = cost_t[truck1][k] - cost_t[truck1][i] + cost_t[truck2][i] - cost_t[truck2][k];
							if ((temp_delta < min_delta) && ((last_swaped_t[truck1][truck2] + p1 * N < iter) || (temp_delta + _current_value < _current_minf))) {
								min_delta = temp_delta;
								indoor1 = i;
								indoor2 = k;
								select1 = j;
								select2 = l;
								flagY = true;
							}
						}
					}
				}
			}
		}
		if (select1 != -1 && select2 != -1) {
			if (!flagY) {
				int temp1, temp2;
				_current_value += min_delta;
				temp1 = _X[indoor1].truck[select1];
				temp2 = _X[indoor2].truck[select2];
				last_swaped_s[temp1][temp2] = iter;
				last_swaped_s[temp2][temp1] = iter;
				//cout << temp1 << " " << temp2 << " s" << endl;
				//cout_solution(_vector_s, _vector_t);
				if (temp1 == M || temp2 == M) {//shift_move
					if (temp1 == M) {//door1插入temp2，indoor2移除temp2
						_X[indoor1].truck[_X[indoor1].num_truck] = temp2;
						_X[indoor1].num_truck += 1;
						_X[indoor1].truck[_X[indoor1].num_truck] = M;
						for (int i = select2; i < _X[indoor2].num_truck; i++) {
							_X[indoor2].truck[i] = _X[indoor2].truck[i + 1];
						}
						_X[indoor2].num_truck += -1;
						_Tot_Ss[indoor1] += _Ss[temp2] - _Ss[temp1];/* 更新容量 */
						_Tot_Ss[indoor2] += _Ss[temp1] - _Ss[temp2];/* 更新容量 */
						_vector_s[temp2] = indoor1;
						update_T_after_change_S(N, J, cost_t, temp2, indoor2, indoor1);
					}
					else {//door2插入temp1，indoor1移除temp1
						_X[indoor2].truck[_X[indoor2].num_truck] = temp1;
						_X[indoor2].num_truck += 1;
						_X[indoor2].truck[_X[indoor2].num_truck] = M;
						for (int i = select1; i < _X[indoor1].num_truck; i++) {
							_X[indoor1].truck[i] = _X[indoor1].truck[i + 1];
						}
						_X[indoor1].num_truck += -1;
						_Tot_Ss[indoor1] += _Ss[temp2] - _Ss[temp1];/* 更新容量 */
						_Tot_Ss[indoor2] += _Ss[temp1] - _Ss[temp2];/* 更新容量 */
						_vector_s[temp1] = indoor2;
						update_T_after_change_S(N, J, cost_t, temp1, indoor1, indoor2);
					}
				}
				else {//swap_move
					_X[indoor1].truck[select1] = temp2;
					_X[indoor2].truck[select2] = temp1;
					_Tot_Ss[indoor1] += _Ss[temp2] - _Ss[temp1];/* 更新容量 */
					_Tot_Ss[indoor2] += _Ss[temp1] - _Ss[temp2];/* 更新容量 */
					_vector_s[temp1] = indoor2;/* 更新解 */
					_vector_s[temp2] = indoor1;
					update_T_after_change_S(N, J, cost_t, temp1, indoor1, indoor2);
					update_T_after_change_S(N, J, cost_t, temp2, indoor2, indoor1);

				}
			}
			else {
				int temp1, temp2;
				_current_value += min_delta;
				temp1 = _Y[indoor1].truck[select1];
				temp2 = _Y[indoor2].truck[select2];
				last_swaped_t[temp1][temp2] = iter;
				last_swaped_t[temp2][temp1] = iter;
				//cout << temp1 << " " << temp2 << " t" << endl;
				//cout_solution(_vector_s, _vector_t);
				if (temp1 == N || temp2 == N) {//shift_move
					if (temp1 == N) {//door1插入temp2，indoor2移除temp2
						_Y[indoor1].truck[_Y[indoor1].num_truck] = temp2;
						_Y[indoor1].num_truck += 1;
						_Y[indoor1].truck[_Y[indoor1].num_truck] = N;
						for (int i = select2; i < _Y[indoor2].num_truck; i++) {
							_Y[indoor2].truck[i] = _Y[indoor2].truck[i + 1];
						}
						_Y[indoor2].num_truck += -1;
						_Tot_Sr[indoor1] += _Sr[temp2] - _Sr[temp1];/* 更新容量 */
						_Tot_Sr[indoor2] += _Sr[temp1] - _Sr[temp2];/* 更新容量 */
						_vector_t[temp2] = indoor1;
						update_S_after_change_T(M, I, cost_s, temp2, indoor2, indoor1);
					}
					else {//door2插入temp1，indoor1移除temp1
						_Y[indoor2].truck[_Y[indoor2].num_truck] = temp1;
						_Y[indoor2].num_truck += 1;
						_Y[indoor2].truck[_Y[indoor2].num_truck] = N;
						for (int i = select1; i < _Y[indoor1].num_truck; i++) {
							_Y[indoor1].truck[i] = _Y[indoor1].truck[i + 1];
						}
						_Y[indoor1].num_truck += -1;
						_Tot_Sr[indoor1] += _Sr[temp2] - _Sr[temp1];/* 更新容量 */
						_Tot_Sr[indoor2] += _Sr[temp1] - _Sr[temp2];/* 更新容量 */
						_vector_t[temp1] = indoor2;
						update_S_after_change_T(M, I, cost_s, temp1, indoor1, indoor2);
					}
				}
				else {//swap_move
					_Y[indoor1].truck[select1] = temp2;
					_Y[indoor2].truck[select2] = temp1;
					_Tot_Sr[indoor1] += _Sr[temp2] - _Sr[temp1];/* 更新容量 */
					_Tot_Sr[indoor2] += _Sr[temp1] - _Sr[temp2];/* 更新容量 */
					_vector_t[temp1] = indoor2;/* 更新解 */
					_vector_t[temp2] = indoor1;
					update_S_after_change_T(M, I, cost_s, temp1, indoor1, indoor2);
					update_S_after_change_T(M, I, cost_s, temp2, indoor2, indoor1);

				}
			}
			//update temp_best solution and best_solution found so far
			if (_current_value < _current_minf) {
				_current_minf = _current_value;
				if (_current_value < _current_best) {
					_current_best = _current_value;
					for (int m = 0; m < M; m++) {
						temp_best_s[m] = _vector_s[m];
					}
					for (int m = 0; m < N; m++) {
						temp_best_t[m] = _vector_t[m];
					}
				}
				if (_current_value < _best_value) {

					_best_value = _current_value;

					ed = clock();

					for (int m = 0; m < M; m++) {
						_best_s[m] = _vector_s[m];
					}
					for (int m = 0; m < N; m++) {
						_best_t[m] = _vector_t[m];
					}
					cout << _best_value <<endl;//<< " F " << compute_value2(Flow, Distance, _best_s, _best_t, M, N, I, J) << 
					SP tmp;
					tmp.obj = _best_value;
					tmp.time = (double)(ed - op) / (double)CLOCKS_PER_SEC;
					SearchProgress.push_back(tmp);
				}
			}

			//triger the landscape detection 
			if (iter == previous_iter + gap_iter) {
				determine_pertubation_length3(pertubation_l, flag, gap_iter, _previous_best, _current_best, previous_pertubation_best, _best_value);
				previous_iter = iter;
				_previous_best = _current_best;
				//_current_best = _current_value;
				//cout << "yess" << endl;
				if (pertubation_l != 0) {
					previous_pertubation_best = _best_value;
					for (int m = 0; m < M; m++) {
						_vector_s[m] = temp_best_s[m];
					}
					for (int m = 0; m < N; m++) {
						_vector_t[m] = temp_best_t[m];
					}
					pertubation(pertubation_l, _vector_s, _vector_t, M, I, N, J);

					_current_value = compute_value2(Flow, Distance, _vector_s, _vector_t, M, N, I, J);
					_current_best = _current_value;
					//cout << "the value after perturbation=" << _current_value << endl;

					//previous_best = _current_value;// this should not be updated
					compute_Tot_capacity(_vector_s, M, I, _Tot_Ss, _Ss);
					compute_Tot_capacity(_vector_t, N, J, _Tot_Sr, _Sr);
					initial_cost(_vector_s, _vector_t, M, I, N, J, Flow, Distance, cost_s, cost_t);
					Vector_to_Partion(_X, _vector_s, M, I);
					Vector_to_Partion(_Y, _vector_t, N, J);
					//cout << "pertubation S" << gap_iter << endl;
					//cout << test_solution(Flow, Distance, _vector_s, _vector_t, _Rs, _Ss, _Rr, _Sr,  M, N, I, J) << " S " << gap_iter<<endl;
				}
			}
		}
		else {
			rand_swap_door_on_side(_vector_s, M, I);
			_current_value = compute_value2(Flow, Distance, _vector_s, _vector_t, M, N, I, J);
			compute_Tot_capacity(_vector_s, M, I, _Tot_Ss, _Ss);
			initial_cost(_vector_s, _vector_t, M, I, N, J, Flow, Distance, cost_s, cost_t);
			Vector_to_Partion(_X, _vector_s, M, I);
			
		}
		if ((long)((clock() - op) / CLOCKS_PER_SEC) > t_all)
			break;
	}
	iteration += depth;
}

void TabuSearch(int *_vector_s, int *_vector_t, long long& _current_value, long long& _current_best, long long& _previous_best, long long& previous_pertubation_best, long long&_best_value, int *_best_s, int *_best_t, partion * _X, partion * _Y, int M, int I, int N, int J,
	int * _Tot_Ss, int *_Ss, int *_Rs, int * _Tot_Sr, int *_Sr, int *_Rr, int(&cost_s)[MAX][MAX], int(&cost_t)[MAX][MAX],
	int ** last_swaped_s, int ** last_swaped_t, int & iteration, int depth, double &op, double &ed, int &previous_iter, long t_all) {
	long long _current_minf;
	_current_minf = _current_value;
	int indoor1, indoor2, truck1, truck2, select1, select2, temp_delta, min_delta;
	bool flagY = false;
	Vector_to_Partion(_X, _vector_s, M, I);
	Vector_to_Partion(_Y, _vector_t, N, J);
	for (int iter = iteration; iter < (iteration + depth); iter++) {
		min_delta = 999999999;
		select1 = -1;
		select2 = -1;
		flagY = false;
		for (int i = 0; i < I; i++) {
			for (int j = 0; j <= _X[i].num_truck; j++) {
				truck1 = _X[i].truck[j];
				for (int k = i + 1; k < I; k++) {
					for (int l = 0; l <= _X[k].num_truck; l++) {
						truck2 = _X[k].truck[l];
						if ((truck1 != truck2) && (_Tot_Ss[i] + _Ss[truck2] - _Ss[truck1] <= _Rs[i]) && (_Tot_Ss[k] + _Ss[truck1] - _Ss[truck2] <= _Rs[k])) {
							temp_delta = cost_s[truck1][k] - cost_s[truck1][i] + cost_s[truck2][i] - cost_s[truck2][k];
							if ((temp_delta < min_delta) && ((last_swaped_s[truck1][truck2] + p1 * M < iter) || (temp_delta + _current_value < _current_minf))) {
								min_delta = temp_delta;
								indoor1 = i;
								indoor2 = k;
								select1 = j;
								select2 = l;
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < J; i++) {
			for (int j = 0; j <= _Y[i].num_truck; j++) {
				truck1 = _Y[i].truck[j];
				for (int k = i + 1; k < J; k++) {
					for (int l = 0; l <= _Y[k].num_truck; l++) {
						truck2 = _Y[k].truck[l];
						if ((truck1 != truck2) && (_Tot_Sr[i] + _Sr[truck2] - _Sr[truck1] <= _Rr[i]) && (_Tot_Sr[k] + _Sr[truck1] - _Sr[truck2] <= _Rr[k])) {
							temp_delta = cost_t[truck1][k] - cost_t[truck1][i] + cost_t[truck2][i] - cost_t[truck2][k];
							if ((temp_delta < min_delta) && ((last_swaped_t[truck1][truck2] + p1 * N < iter) || (temp_delta + _current_value < _current_minf))) {
								min_delta = temp_delta;
								indoor1 = i;
								indoor2 = k;
								select1 = j;
								select2 = l;
								flagY = true;
							}
						}
					}
				}
			}
		}
		if (select1 != -1 && select2 != -1) {
			if (!flagY) {
				int temp1, temp2;
				_current_value += min_delta;
				temp1 = _X[indoor1].truck[select1];
				temp2 = _X[indoor2].truck[select2];
				last_swaped_s[temp1][temp2] = iter;
				last_swaped_s[temp2][temp1] = iter;
				//cout << temp1 << " " << temp2 << " s" << endl;
				//cout_solution(_vector_s, _vector_t);
				if (temp1 == M || temp2 == M) {//shift_move
					if (temp1 == M) {//door1插入temp2，indoor2移除temp2
						_X[indoor1].truck[_X[indoor1].num_truck] = temp2;
						_X[indoor1].num_truck += 1;
						_X[indoor1].truck[_X[indoor1].num_truck] = M;
						for (int i = select2; i < _X[indoor2].num_truck; i++) {
							_X[indoor2].truck[i] = _X[indoor2].truck[i + 1];
						}
						_X[indoor2].num_truck += -1;
						_Tot_Ss[indoor1] += _Ss[temp2] - _Ss[temp1];/* 更新容量 */
						_Tot_Ss[indoor2] += _Ss[temp1] - _Ss[temp2];/* 更新容量 */
						_vector_s[temp2] = indoor1;
						update_T_after_change_S(N, J, cost_t, temp2, indoor2, indoor1);
					}
					else {//door2插入temp1，indoor1移除temp1
						_X[indoor2].truck[_X[indoor2].num_truck] = temp1;
						_X[indoor2].num_truck += 1;
						_X[indoor2].truck[_X[indoor2].num_truck] = M;
						for (int i = select1; i < _X[indoor1].num_truck; i++) {
							_X[indoor1].truck[i] = _X[indoor1].truck[i + 1];
						}
						_X[indoor1].num_truck += -1;
						_Tot_Ss[indoor1] += _Ss[temp2] - _Ss[temp1];/* 更新容量 */
						_Tot_Ss[indoor2] += _Ss[temp1] - _Ss[temp2];/* 更新容量 */
						_vector_s[temp1] = indoor2;
						update_T_after_change_S(N, J, cost_t, temp1, indoor1, indoor2);
					}
				}
				else {//swap_move
					_X[indoor1].truck[select1] = temp2;
					_X[indoor2].truck[select2] = temp1;
					_Tot_Ss[indoor1] += _Ss[temp2] - _Ss[temp1];/* 更新容量 */
					_Tot_Ss[indoor2] += _Ss[temp1] - _Ss[temp2];/* 更新容量 */
					_vector_s[temp1] = indoor2;/* 更新解 */
					_vector_s[temp2] = indoor1;
					update_T_after_change_S(N, J, cost_t, temp1, indoor1, indoor2);
					update_T_after_change_S(N, J, cost_t, temp2, indoor2, indoor1);

				}
			}
			else {
				int temp1, temp2;
				_current_value += min_delta;
				temp1 = _Y[indoor1].truck[select1];
				temp2 = _Y[indoor2].truck[select2];
				last_swaped_t[temp1][temp2] = iter;
				last_swaped_t[temp2][temp1] = iter;
				//cout << temp1 << " " << temp2 << " t" << endl;
				//cout_solution(_vector_s, _vector_t);
				if (temp1 == N || temp2 == N) {//shift_move
					if (temp1 == N) {//door1插入temp2，indoor2移除temp2
						_Y[indoor1].truck[_Y[indoor1].num_truck] = temp2;
						_Y[indoor1].num_truck += 1;
						_Y[indoor1].truck[_Y[indoor1].num_truck] = N;
						for (int i = select2; i < _Y[indoor2].num_truck; i++) {
							_Y[indoor2].truck[i] = _Y[indoor2].truck[i + 1];
						}
						_Y[indoor2].num_truck += -1;
						_Tot_Sr[indoor1] += _Sr[temp2] - _Sr[temp1];/* 更新容量 */
						_Tot_Sr[indoor2] += _Sr[temp1] - _Sr[temp2];/* 更新容量 */
						_vector_t[temp2] = indoor1;
						update_S_after_change_T(M, I, cost_s, temp2, indoor2, indoor1);
					}
					else {//door2插入temp1，indoor1移除temp1
						_Y[indoor2].truck[_Y[indoor2].num_truck] = temp1;
						_Y[indoor2].num_truck += 1;
						_Y[indoor2].truck[_Y[indoor2].num_truck] = N;
						for (int i = select1; i < _Y[indoor1].num_truck; i++) {
							_Y[indoor1].truck[i] = _Y[indoor1].truck[i + 1];
						}
						_Y[indoor1].num_truck += -1;
						_Tot_Sr[indoor1] += _Sr[temp2] - _Sr[temp1];/* 更新容量 */
						_Tot_Sr[indoor2] += _Sr[temp1] - _Sr[temp2];/* 更新容量 */
						_vector_t[temp1] = indoor2;
						update_S_after_change_T(M, I, cost_s, temp1, indoor1, indoor2);
					}
				}
				else {//swap_move
					_Y[indoor1].truck[select1] = temp2;
					_Y[indoor2].truck[select2] = temp1;
					_Tot_Sr[indoor1] += _Sr[temp2] - _Sr[temp1];/* 更新容量 */
					_Tot_Sr[indoor2] += _Sr[temp1] - _Sr[temp2];/* 更新容量 */
					_vector_t[temp1] = indoor2;/* 更新解 */
					_vector_t[temp2] = indoor1;
					update_S_after_change_T(M, I, cost_s, temp1, indoor1, indoor2);
					update_S_after_change_T(M, I, cost_s, temp2, indoor2, indoor1);

				}
			}
			//update temp_best solution and best_solution found so far
			if (_current_value < _current_minf) {
				_current_minf = _current_value;
				if (_current_value < _current_best) {
					_current_best = _current_value;
					for (int m = 0; m < M; m++) {
						temp_best_s[m] = _vector_s[m];
					}
					for (int m = 0; m < N; m++) {
						temp_best_t[m] = _vector_t[m];
					}
				}
				if (_current_value < _best_value) {

					_best_value = _current_value;

					ed = clock();

					for (int m = 0; m < M; m++) {
						_best_s[m] = _vector_s[m];
					}
					for (int m = 0; m < N; m++) {
						_best_t[m] = _vector_t[m];
					}
					//cout << _best_value << " F " << compute_value2(Flow, Distance, _best_s, _best_t, M, N, I, J) << endl;
					cout << _best_value << endl;//<< " F " << compute_value2(Flow, Distance, _best_s, _best_t, M, N, I, J) << 
					SP tmp;
					tmp.obj = _best_value;
					tmp.time = (double)(ed - op) / (double)CLOCKS_PER_SEC;
					SearchProgress.push_back(tmp);
				}
			}
			if (iter == previous_iter + gap_iter) {
				previous_iter =iter;
				if (gap_iter == weak_depth) {
					pertubation_l = std::max(1, (int)(r1*I));
				}
				else if (gap_iter == medium_depth) {
					pertubation_l = std::max(1, (int)(r2*I));
				}
				else
				{
					pertubation_l = std::max(1, (int)(r3*I));
				}
				pertubation(pertubation_l, _vector_s, _vector_t, M, I, N, J);
				_current_value = compute_value2(Flow, Distance, _vector_s, _vector_t, M, N, I, J);
				_current_best = _current_value;
				//cout << "the value after perturbation=" << _current_value << " "<<gap_iter<<endl;

				//previous_best = _current_value;// this should not be updated
				compute_Tot_capacity(_vector_s, M, I, _Tot_Ss, _Ss);
				compute_Tot_capacity(_vector_t, N, J, _Tot_Sr, _Sr);
				initial_cost(_vector_s, _vector_t, M, I, N, J, Flow, Distance, cost_s, cost_t);
				Vector_to_Partion(_X, _vector_s, M, I);
				Vector_to_Partion(_Y, _vector_t, N, J);
				//cout << "pertubation S" << gap_iter << endl;
				//cout << test_solution(Flow, Distance, _vector_s, _vector_t, _Rs, _Ss, _Rr, _Sr,  M, N, I, J) << " S " << gap_iter<<endl;
			}
		}
		else {
			rand_swap_door_on_side(_vector_s, M, I);
			_current_value = compute_value2(Flow, Distance, _vector_s, _vector_t, M, N, I, J);
			compute_Tot_capacity(_vector_s, M, I, _Tot_Ss, _Ss);
			initial_cost(_vector_s, _vector_t, M, I, N, J, Flow, Distance, cost_s, cost_t);
			Vector_to_Partion(_X, _vector_s, M, I);

		}
		if ((long)((clock() - op) / CLOCKS_PER_SEC) > t_all)
			break;
	}
	//iteration += depth;
}

int Permutation_local_search_start_from_s_with_depth(int _row, int _col, long long & _current_value, long long & _minf, long long & _best_value, int * _s, int * _t, int * temp_best_s, int * temp_best_t,
	int(&matrix_a)[1000][1000], int(&matrix_b)[1000][1000], int(&g)[1000][1000], double & ed, double t_all, int _depth)
{
	int t1, t2, count;
	int break_num = 0;
	int _local_optimal = _current_value;
	_minf = _current_value;
	//_best_value= _current_value;
	//cout << "_best_value" << _best_value << endl;
	count = 0;
	for (int i = 0; i < _row; i++) {
		temp_best_s[i] = _s[i];
	}
	for (int i = 0; i < _col; i++) {
		temp_best_t[i] = _t[i];
	}
	int break_num1 = 0, break_num2 = 0;
	while (count < _depth)
	{
		break_num = 0;
		//find best_t, and use the best_s
		for (int j = 0; j < _col; j++)
		{
			for (int f = 0; f < _col; f++)
			{
				g[j][f] = 0;
				for (int i = 0; i < _row; i++)
				{
					g[j][f] += matrix_a[i][j] * matrix_b[_s[i]][f];
				}
			}
		}
		//_current_value = solve_LAP(g, _t, _col, pthreads);
		_current_value = Hungarian_algorithm(g, _col, _t);
		if (_current_value < _local_optimal)
		{
			_local_optimal = _current_value;
			break_num++;
			break_num1 = 1;
			if (_current_value < _minf) {
				_minf = _current_value;
				for (int m = 0; m < _row; m++) {
					temp_best_s[m] = _s[m];
				}
				for (int m = 0; m < _row; m++) {
					temp_best_t[m] = _t[m];
				}
			}
			if (_current_value < _best_value)
			{
				ed = clock();
				_best_value = _current_value;
				cout << _best_value << endl;
			}
		}
		else break_num1 = 0;
		/*if (break_num == 0)
		{
			t1 = rand() % _col;
			t2 = rand() % _col;
			while (t1 == t2) {
				t2 = rand() % _col;
			}
			swap(_t[t1], _t[t2]);
			count++;
			_local_optimal = 999999999;
			break_num = 1;
		}*/

		for (int i = 0; i < _row; i++)
		{
			for (int k = 0; k < _row; k++)
			{
				g[i][k] = 0;
				for (int j = 0; j < _col; j++)
				{
					g[i][k] += matrix_a[i][j] * matrix_b[k][_t[j]];
				}
			}
		}
		//_current_value = solve_LAP(g, _s, _row, pthreads);
		_current_value = Hungarian_algorithm(g, _row, _s);
		if (_current_value < _local_optimal)
		{
			_local_optimal = _current_value;
			break_num2 = 1;
			if (_current_value < _minf) {
				_minf = _current_value;
				for (int m = 0; m < _row; m++) {
					temp_best_s[m] = _s[m];
				}
				for (int m = 0; m < _row; m++) {
					temp_best_t[m] = _t[m];
				}
			}
			if (_current_value < _best_value)
			{
				ed = clock();
				_best_value = _current_value;
				cout << _best_value << endl;
			}

		}
		else break_num2 = 0;

		if (break_num1 == 0 && break_num2 == 0) {
			t1 = rand() % _row;
			t2 = rand() % _row;
			while (t1 == t2) {
				t2 = rand() % _row;
			}
			swap(_s[t1], _s[t2]);
			t1 = rand() % _col;
			t2 = rand() % _col;
			while (t1 == t2) {
				t2 = rand() % _col;
			}
			swap(_t[t1], _t[t2]);
			count++;
			_local_optimal = 999999999;
		}
		/*if (break_num == 1)
		{
			t1 = rand() % _row;
			t2 = rand() % _row;
			while (t1 == t2) {
				t2 = rand() % _row;
			}
			swap(_s[t1], _s[t2]);
			count++;
			_local_optimal = 999999999;

		}*/
		if ((long)((clock() - op) / CLOCKS_PER_SEC) >= t_all)
		{
			return 0;
		}
	}
	return 0;
}

void compute_C_according_t(int M, int N, int I, int J, int(&Flow)[1000][1000], int(&Distance)[1000][1000], int(&C)[1000][1000], int * _t) {
	for (int m = 0; m < M; m++)
	{
		for (int i = 0; i < I; i++)
		{
			C[m][i] = 0;
			for (int n = 0; n < N; n++)
			{
				C[m][i] += Flow[m][n] * Distance[i][_t[n]];
			}
		}
	}
}
void compute_C_according_s(int M, int N, int I, int J, int(&Flow)[1000][1000], int(&Distance)[1000][1000], int(&C)[1000][1000], int * _s) {
	for (int n = 0; n < N; n++)
	{
		for (int j = 0; j < J; j++)
		{
			C[n][j] = 0;
			for (int m = 0; m < M; m++)
			{
				C[n][j] += Distance[_s[m]][j] * Flow[m][n];
			}
		}
	}
}



void tabu3New(int * p_s, int *p_t, int* temp_p_s, int* temp_p_t, int &iteration, int depth) {//不重置previous iter和previous_best
	long long local_optimal, temp_value;
	long long current_best, previous_value, previous_pertubation_best;
	int iter_without_improvement = 0, previous_iter = 0;
	int enlarge = 0;
	//temp_value = 99999999;
	depth = 1999999999;
	/* generate initial solution */
	generate_initial_solution(Rs, Ss, s, M, I, rand() % 2 + 8);
	generate_initial_solution(Rr, Sr, t, N, J, rand() % 2 + 8);

	compute_Tot_capacity(s, M, I, Tot_Ss, Ss);
	compute_Tot_capacity(t, N, J, Tot_Sr, Sr);
	initial_cost(s, t, M, I, N, J, Flow, Distance, cost_s, cost_t);
	Vector_to_Partion(X, s, M, I);
	Vector_to_Partion(Y, t, N, J);
	iteration = 200;
	//depth = 100000;
	previous_value = 999999999;
	previous_pertubation_best = 999999999;
	best_minf = 999999999;

	

	previous_value = now_f;
	previous_pertubation_best = now_f;
	local_optimal = now_f;
	current_best = now_f;
	LandScapeDetectionTabuSearch(s, t, now_f, current_best, previous_value, previous_pertubation_best, best_minf, best_s, best_t, X, Y, M, I, N, J, Tot_Ss, Ss,
		Rs, Tot_Sr, Sr, Rr, cost_s, cost_t, last_swaped_s, last_swaped_t, iteration, depth, op, ed, previous_iter, t_all);

}

void tabu(int * p_s, int *p_t, int* temp_p_s, int* temp_p_t, int &iteration, int depth) {//不重置previous iter和previous_best
	long long local_optimal, temp_value;
	long long current_best, previous_value, previous_pertubation_best;
	int iter_without_improvement = 0, previous_iter = 0;
	int enlarge = 0;
	//temp_value = 99999999;
	depth = 1999999999;
	/* generate initial solution */
	generate_initial_solution(Rs, Ss, s, M, I, rand() % 2 + 8);
	generate_initial_solution(Rr, Sr, t, N, J, rand() % 2 + 8);
	now_f = compute_value2(Flow, Distance, s, t, M, N, I, J);
	/* end */

	compute_Tot_capacity(s, M, I, Tot_Ss, Ss);
	compute_Tot_capacity(t, N, J, Tot_Sr, Sr);
	initial_cost(s, t, M, I, N, J, Flow, Distance, cost_s, cost_t);
	Vector_to_Partion(X, s, M, I);
	Vector_to_Partion(Y, t, N, J);
	iteration = 200;
	previous_value = now_f;
	previous_pertubation_best = now_f;
	local_optimal = now_f;
	current_best = now_f;
	TabuSearch(s, t, now_f, current_best, previous_value, previous_pertubation_best, best_minf, best_s, best_t, X, Y, M, I, N, J, Tot_Ss, Ss,
		Rs, Tot_Sr, Sr, Rr, cost_s, cost_t, last_swaped_s, last_swaped_t, iteration, depth, op, ed, previous_iter, t_all);

}

void hybrid_search3New(int * p_s, int *p_t, int* temp_p_s, int* temp_p_t, int &iteration, int depth) {//不重置previous iter和previous_best
	long long local_optimal, temp_value;
	long long current_best, previous_value, previous_pertubation_best;
	int iter_without_improvement = 0, previous_iter = 0;
	int enlarge = 0;
	//temp_value = 99999999;

	/* generate initial solution */
	generate_initial_solution(Rs, Ss, s, M, I, rand() % 2 + 8);
	generate_initial_solution(Rr, Sr, t, N, J, rand() % 2 + 8);
	now_f = compute_value2(Flow, Distance, s, t, M, N, I, J);
	/* end */

	compute_Tot_capacity(s, M, I, Tot_Ss, Ss);
	compute_Tot_capacity(t, N, J, Tot_Sr, Sr);
	initial_cost(s, t, M, I, N, J, Flow, Distance, cost_s, cost_t);
	Vector_to_Partion(X, s, M, I);
	Vector_to_Partion(Y, t, N, J);
	iteration = 200; 
	//depth = 100000;
	

	previous_value = now_f;
	previous_pertubation_best = now_f;
	local_optimal = now_f;
	current_best = now_f;
	LandScapeDetectionTabuSearch(s, t, now_f, current_best, previous_value, previous_pertubation_best, best_minf, best_s, best_t, X, Y, M, I, N, J, Tot_Ss, Ss,
		Rs, Tot_Sr, Sr, Rr, cost_s, cost_t, last_swaped_s, last_swaped_t, iteration, depth, op, ed, previous_iter, t_all);
	if (now_f < local_optimal) {
		local_optimal = now_f;
	}
	//cout << test_solution(Flow, Distance, s, t, Rs, Ss, Rr, Sr, M, N, I, J) << "tabu1" << endl;
	//cout_solution(_vector_s, _vector_t);
	while (1) {
		generate_new_flow(s, t, M, N);
		for (int l = 0; l < I; l++) {
			p_s[l] = l;
		}
		for (int l = 0; l < J; l++) {
			p_t[l] = l;
		}
		temp_value = now_f;
		Permutation_local_search_start_from_s_with_depth(I, J, temp_value, now_f, best_minf, p_s, p_t, temp_p_s, temp_p_t, new_flow, Distance, C, ed, t_all, depth_LS);
		Permutation_To_Vector(temp_p_s, temp_p_t, s, t, M, N);
		//cout << test_solution(Flow, Distance, _vector_s, _vector_t, _Rs, _Ss, _Rr, _Sr, M, N, I, J) << "permutation" << endl;
		if (now_f < current_best) {
			current_best = now_f;
		}
		if (now_f == best_minf) {
			for (int i = 0; i < M; i++) {
				best_s[i] = s[i];
			}
			for (int i = 0; i < N; i++) {
				best_t[i] = t[i];
			}
		}
		compute_Tot_capacity(s, M, I, Tot_Ss, Ss);
		compute_Tot_capacity(t, N, J, Tot_Sr, Sr);
		initial_cost(s, t, M, I, N, J, Flow, Distance, cost_s, cost_t);
		//previous_value = now_f;
		//previous_pertubation_best = best_minf;
		Vector_to_Partion(X, s, M, I);
		Vector_to_Partion(Y, t, N, J);
		LandScapeDetectionTabuSearch(s, t, now_f, current_best, previous_value, previous_pertubation_best, best_minf, best_s, best_t, X, Y, M, I, N, J, Tot_Ss, Ss,
			Rs, Tot_Sr, Sr, Rr, cost_s, cost_t, last_swaped_s, last_swaped_t, iteration, depth, op, ed,  previous_iter, t_all);
		//cout << test_solution(Flow, Distance, s, t, Rs, Ss, Rr, Sr, M, N, I, J) << "tabu2" << endl;

		if ((long)((clock() - op) / CLOCKS_PER_SEC) > t_all)
			break;
	}
}
void exportSearchProgress(char * filename)
{
	std::ofstream myfile(filename, std::ios::app);
	for (SP state : SearchProgress)
		myfile << state.obj << ";" << state.time << std::endl;
	myfile << "****" << std::endl;
}
void testGetOpt(int argc, char *argv[]) {
	int opt;  // getopt() 的返回值
	const char *optstring = "i:a:b:c:d:e:"; // 

	while ((opt = getopt(argc, argv, optstring)) != -1)
	{

		switch (opt)
		{
		case 'i':
			resultpath = optarg;
			cout << "testfile= " << resultpath << endl;
			break;
		case 'a':
			sscanf(optarg, "%d", &depthall);
			cout << "depth_out= " << depthall << endl;
			break;
		case 'b':
			sscanf(optarg, "%d", &weak_depth);
			cout << "weak_depth= " << weak_depth << endl;
			break;
		case 'c':
			sscanf(optarg, "%d", &medium_depth);
			cout << "medium_depth= " << medium_depth << endl;
			break;
		case 'd':
			sscanf(optarg, "%d", &strong_depth);
			cout << "strong_depth= " << strong_depth << endl;
			break;
		case 'e':
			sscanf(optarg, "%d", &TTP);
			cout << "pmax= " << TTP << endl;
			break;
		}
	}
}
int main(int argc, char ** argv)
{
	//load_parameter(parameter);
	testGetOpt(argc, argv);
	double totaltime = 0.0, argtime = 0.0;
	long long totalvalue = 0;

	iteration = 0;
	int p_s[MAX], p_t[MAX], temp_p_s[MAX], temp_p_t[MAX];
	bool subpool_index[MAXPOOL];//determine if solution[i] was used to guide the pathrelink, 1 means using, 0 means  not using.

	char test[20];
	//snprintf(test, sizeof(test), "%d", select_a);
	strcat(resultfilename, test);
	strcat(resultfilename, resultpath);
	cout << resultfilename << endl;
	strcat(filename, resultpath);
	//strcat(poolfilename, test);
	//strcat(poolfilename, resultpath);
	cout << filename << endl;
	printf("CLOCKS_PER_SEC:%d\n", CLOCKS_PER_SEC);
	printf("ReStartTimes:%d\n", ReStartTimes);
	printf("t_all:%ds\n", t_all);
	ReadInstance(filename);
	srand((int)time(NULL));
	allocate_memory();
	op = clock();
	ed = clock();
	
	int enlarge = 0, temp_value;
	temp_value = 99999999;
	
	/* generate initial solution */
	for (int iter = 0; iter < ReStartTimes; iter++)//ReStart 10 times
	{
		best_minf = 9999999999;
		op = clock();
		ed = clock();
		init();
		iteration = 100;
		gap_iter = weak_depth;

	
		hybrid_search3New(p_s, p_t, temp_p_s, temp_p_t, iteration, depthall);
		//tabu3New(p_s, p_t, temp_p_s, temp_p_t, iteration, depthall);
		//tabu(p_s, p_t, temp_p_s, temp_p_t, iteration, depthall);
		t_f_best = (ed - op) / CLOCKS_PER_SEC;
		totaltime += t_f_best;
		totalvalue += best_minf;
		//exportSearchProgress(resultfilename);
		FILE *fp = fopen(resultfilename, "a+");
		display_result(fp);
		fclose(fp);
		printf("\n%d Times ReStart\n", (iter + 1));
		printf("time:%lfs\n", t_f_best);
		printf("best_minf:%lld\n", best_minf);

	}
	FILE *fp = fopen(resultfilename, "a+");
	fprintf(fp, "ArgT:%lfs argtime:%lfs\n", totaltime, totaltime / ReStartTimes);
	fprintf(fp, "ArgO:%lld argvalue:%lf\n", totalvalue, (double)(totalvalue*1.0 / ReStartTimes));
	fclose(fp);
	delete_memory();
	return 0;
}

