#include <stdlib.h>
#include <iostream>
#include <math.h> 
#include<ctime>
#include <fstream>

// use modulus for indexing for periodicity
using namespace std;

void initialize_surface(int surf[], int L);
void reduce_surface(int surf[], int L);
void print_surface(int surf[], int L);
void st_deposit(int surf[], int L, int t);
int st_check(int surf[], int L, int i);
bool is_kink(int surf[], int L, int i);
double global_width(int surf[], int L);
double surface_mean(int surf[], int L);
void measure_width(int t0, int tau, double lambda, int rep, int L);
double global_width_avg(int surf[], int L);

int main() {
	srand(time(0));

	int L = 8192, t0 = 1e2, tau = 1e5, rep = 25e1;
	double lambda = 1.01;
	
	measure_width(t0, tau, lambda, rep, L);

	return 0;
}

void initialize_surface(int surf[], int L) {
	// initializes flat surface - all zeros
	for (int i = 0; i < L; i++)
		surf[i] = 0;
}

void print_surface(int surf[], int L) {
	for (int i = 0; i < L; i++)
		cout << surf[i] << " ";
	cout << endl;
}

void reduce_surface(int surf[], int L) {
	// simplifies surface such that the minimum height is zero.
	// relative heights are unaffected.	
	int min = surf[0];
	// find minimum
	for (int i = 1; i < L; i++)
		if (surf[i] < min) min = surf[i];
	// subtract minimum
	for (int i = 0; i < L; i++)
		surf[i] -= min;
}

void st_deposit(int surf[], int L, int t) {
	// t is the number of monolayers (L particles) to be deposited.
	for (int k = 0; k < L*t; k++) {
		// k dummy index for total number of L*t depositions.
		// perform ST rule checks and deposit
		surf[st_check(surf, L, rand() % L)] += 1;
	}
}

int st_check(int surf[], int L, int i) {
	// given a surface and an index, performs rule checks to
	// output a deposition index - either i + 1 or i - 1.
	bool l, r;
	i = i + L; // for periodicity on left hand side, will take mod L later anyways
	// deposit at i if i is kink or valley.
	if (is_kink(surf, L, i%L)) return i%L;
	// perform left check
	if (surf[(i-1) % L] < surf[i%L]) l = true;
	else if (surf[(i-1)%L] == surf[i%L] && is_kink(surf, L, (i-1)%L)) l = true;
	else l = false;
	// perform right check
	if (surf[(i+1) % L] < surf[i%L]) r = true;
	else if (surf[(i+1)%L] == surf[i%L] && is_kink(surf, L, (i+1)%L)) r = true;
	else r = false;
	if (l && r) {
		if (rand()%2 == 0) return (i-1)%L;
		else return (i+1)%L;	
	}
	else if (l) return (i-1)%L;
	else if (r) return (i+1)%L;
	else return i%L;
} 

bool is_kink(int surf[], int L, int i) {
	// checks whether the position i is a kink (or valley)
	i = i + L;
	if ((surf[(i+1)%L] > surf[i%L]) || (surf[(i-1)%L] > surf[i%L])) return true;
	else return false;
}

double global_width(int surf[], int L) {
	// returns the global width of the surface at some time
	// reduce surface before calculations
	reduce_surface(surf, L);
	double mean = surface_mean(surf, L);
	double sum = 0;
	for (int i = 0; i < L; i++) sum += (surf[i] - mean)*(surf[i] - mean);
	return sum;
}

double surface_mean(int surf[], int L) {
	// returns the mean height of surface
	double sum = 0;
	for (int i = 0; i < L; i++) sum += surf[i];
	return sum / L;
}

void measure_width(int t0, int tau, double lambda, int rep, int L) {
	// t0 is initial time
	// tau is the final time
	// rep is the number of repetitions - for width measurements
	// lambda is the proportionality constant between consecutive times
	// L is surface length

	int surf[L];
	int t, dt;
	int N = int(log(tau/t0)/log(lambda)) + 1;
	double width[N];

	cout << "number of data points: " << N << endl;

	for (int i = 0; i < N; i++) {
		width[i] = 0;
	}
	ofstream data;
	data.open("data.txt");

	for (int n = 0; n < rep; n++) {
		initialize_surface(surf, L);
		t = t0;
		st_deposit(surf, L, t0);
		for (int m = 0; m < N; m++) {
			width[m] += global_width(surf, L);
			dt = (lambda-1)*t;
			st_deposit(surf, L, dt);
			t += dt;
		}
		cout << "rep: " << n+1 << endl;
	}
	for (int i = 0; i < N; i++) {
		data << t0 * pow(lambda, i) << ", " << sqrt(width[i]/(L*rep)) << endl;
	}
	data.close();
}














