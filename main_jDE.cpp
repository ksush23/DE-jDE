/*
  CEC20 Test Function Suite for Single Objective Bound Constrained Numerical Optimization
*/
#include <WINDOWS.H>    
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <time.h>
#include <tuple>
#include <iomanip> 
#include <fstream>
#include <string>


using namespace std;


void cec20_test_func(double*, double*, int, int, int);
double* OShift, * M, * y, * z, * x_bound;
int ini_flag = 0, n_flag, func_flag, * SS;

const int NP = 100;
const int D = 10;
const int Gmax = 100000;
const int boundary_min = -100;
const int boundary_max = 100;
const double func_1 = 100.;
const double func_2 = 1100.;
const double func_3 = 700.;
const double func_4 = 1900.;
const double func_5 = 1700.;
const double func_6 = 1600.;
const double func_7 = 2100.;
const double func_8 = 2200.;
const double func_9 = 2400.;
const double func_10 = 2500.;
const double eps = 10e-8;
static int func_n;

std::tuple<vector<double>, double> mutation(vector<vector<double>> P_matrix, int index, double f) {
	double probability = 0;
	double new_f = 0;
	while (probability == 0 || probability == 1) {
		probability = rand() / (RAND_MAX + 1.);
	}
	if (probability <= 0.1) {
		double p = (double)rand() / RAND_MAX;
		f = 0.1 + p * 0.8;
	}
	else {
		new_f = f;
	}

	int j = index;
	while (j == index) {
		j = rand() % NP;
	}
	int k = index;
	while (k == index || k == j) {
		k = rand() % NP;
	}
	int m = index;
	while (m == index || m == k || m == j) {
		m = rand() % NP;
	}
	vector<double> new_vector(D);
	for (int i = 0; i < D; i++) {
		new_vector[i] = P_matrix[k][i] - P_matrix[m][i];
		new_vector[i] *= new_f;
		new_vector[i] += P_matrix[j][i];
		if (new_vector[i] <= boundary_min || new_vector[i] >= boundary_max) {
			double f = (double)rand() / RAND_MAX;
			new_vector[i] = boundary_min + f * (boundary_max - boundary_min);
		}
	}
	return std::make_tuple(new_vector, new_f);
}

std::tuple<vector<double>, double> crossover(vector<double> x, vector<double> v, double cr) {
	double probability = 0;
	double new_cr = 0;
	while (probability == 0 || probability == 1) {
		probability = rand() / (RAND_MAX + 1.);
	}
	if (probability <= 0.1) {
		new_cr = rand() / (RAND_MAX + 1.);
	}
	else {
		new_cr = cr;
	}
	vector<double> crossovered_vector(D);
	for (int i = 0; i < D; i++) {
		probability = 0;
		while (probability == 0 || probability == 1) {
			probability = rand() / (RAND_MAX + 1.);
		}
		if (probability <= new_cr) {
			crossovered_vector[i] = v[i];
		}
		else {
			crossovered_vector[i] = x[i];
		}
	}
	return std::make_tuple(crossovered_vector, new_cr);
}

double function(vector<double> x) {
	double sum = 0;
	for (int i = 0; i < D; i++) {
		//sum += x[i] * x[i];
		//sum += 100 * ((x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1])) + (1 - x[i]) * (1 - x[i]);
		sum += (0 - (x[i] * sin(sqrt(abs(x[i])))));
	}
	return sum;
}

double best_f(vector<vector<double>> P_matrix) {
	double current_best = function(P_matrix[0]);
	double ind[D];
	for (int j = 0; j < D; j++) {
		ind[j] = P_matrix[0][j];
	}
	cec20_test_func(ind, &current_best, D, 1, func_n);
	for (int i = 1; i < NP; i++) {
		double new_best;
		for (int j = 0; j < D; j++) {
			ind[j] = P_matrix[i][j];
		}
		cec20_test_func(ind, &new_best, D, 1, func_n);
		if (new_best < current_best) {
			current_best = new_best;
		}
	}
	return current_best;
}

vector<double> best_x(vector<vector<double>> P_matrix) {
double current_best[D];
double current_best_f;
for (int i = 0; i < D; i++) {
	current_best[i] = P_matrix[0][i];
}
cec20_test_func(current_best, &current_best_f, D, 1, func_n);
for (int i = 1; i < NP; i++) {
	double new_best_f;
	double new_best[D];
	for (int j = 0; j < D; j++) {
		new_best[j] = P_matrix[i][j];
	}
	cec20_test_func(new_best, &new_best_f, D, 1, func_n);

	if (new_best_f < current_best_f) {
		for (int j = 0; j < D; j++) {
			current_best[j] = new_best[j];
		}
	}
}
vector<double> best(D);
for (int i = 0; i < D; i++) {
	best[i] = current_best[i];
}
return best;
}

void main()
{
	for (int kk = 9; kk <= 9; kk++) {
		func_n = kk;
		/*if (kk == 6 || kk == 7) {
			continue;
		}*/
		ofstream outData, outX;
		string file1 = "outfile_biig" + to_string(kk) + ".csv";
		string file2 = "statistics_biig" + to_string(kk) + ".csv";
		outData.open(file1, ios::app);
		outX.open(file2, ios::app);
		for (int kkk = 0; kkk < 30; kkk++) {
			srand(time(NULL));
			double func;
			switch (func_n) {
			case 1:
				func = func_1;
				break;
			case 2:
				func = func_2;
				break;
			case 3:
				func = func_3;
				break;
			case 4:
				func = func_4;
				break;
			case 5:
				func = func_5;
				break;
			case 6:
				func = func_6;
				break;
			case 7:
				func = func_7;
				break;
			case 8:
				func = func_8;
				break;
			case 9:
				func = func_9;
				break;
			case 10:
				func = func_10;
				break;
			}
			vector<vector<double>> P(NP);
			for (int i = 0; i < NP; i++)
				P[i].resize(D);
			for (int i = 0; i < NP; i++) {
				for (int j = 0; j < D; j++) {
					double f = (double)rand() / RAND_MAX;
					P[i][j] = boundary_min + f * (boundary_max - boundary_min);
				}
			}
			int G = 0;
			vector<vector<double>> P_new(NP);
			for (int i = 0; i < NP; i++)
				P_new[i].resize(D);

			vector<double> f_best(Gmax);
			f_best[G] = best_f(P);
			G++;

			vector<double> cr_array(NP);
			for (int i = 0; i < NP; i++) {
				cr_array[i] = 0.5;
			}
			vector<double> f_array(NP);
			for (int i = 0; i < NP; i++) {
				f_array[i] = 0.9;
			}

			while (G < Gmax) {
				for (int i = 0; i < NP; i++) {
					vector<double> x = P[i];
					vector<double> v;
					double f;
					tie(v, f) = mutation(P, i, cr_array[i]);
					vector<double> u;
					double cr;
					tie(u, cr) = crossover(x, v, f_array[i]);

					double f_x = 0, f_u = 0;
					//f_x = function(x);
					//f_u = function(u);
					double ind[D];
					for (int i = 0; i < D; i++) {
						ind[i] = x[i];
					}
					cec20_test_func(ind, &f_x, D, 1, func_n);
					double new_ind[D];
					for (int i = 0; i < D; i++) {
						new_ind[i] = v[i];
					}
					cec20_test_func(new_ind, &f_u, D, 1, func_n);

					if (f_u <= f_x) {
						f_array[i] = f;
						cr_array[i] = cr;
						for (int j = 0; j < D; j++) {
							P_new[i][j] = u[j];
						}
					}
					else {
						for (int j = 0; j < D; j++) {
							P_new[i][j] = x[j];
						}
					}
				}
				for (int i = 0; i < NP; i++) {
					for (int j = 0; j < D; j++) {
						P[i][j] = P_new[i][j];
					}
				}
				f_best[G] = best_f(P);
				/*if (abs(f_best[G] - func) < eps) {
					break;
				}*/
				if (G == 1000 || G == 2000 || G == 3000 || G == 4000 || G == 5000 || G == 6000 || G == 7000 || G == 8000 || G == 9000 || G == 10000 || G == 11000 || G == 20000) {
					cout << f_best[G] << endl;
				}
				G++;
			}
			std::cout << "best x: " << "(";
			vector<double> best_ind = best_x(P);
			for (int i = 0; i < D; i++) {
				std::cout << setprecision(8) << best_ind[i] << ", ";
				outX << setprecision(8) << best_ind[i] << " ";
			}
			std::cout << ")" << endl;
			std::cout << f_best[Gmax - 1] << endl;
			for (int i = 0; i < Gmax; i++) {
				if (f_best[i] == 0) {
					break;
				}
				//std::cout << setprecision(15) << f_best[i] << "  ";
				outData << setprecision(8) << f_best[i] << " ";
			}
			outData << endl;
			outX << endl;
		}
	}
}


