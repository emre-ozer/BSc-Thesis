#include <iostream>
#include <armadillo>
#include<fstream>
#include <random>
#include <algorithm>

using namespace std;
using namespace arma;


vec u_t(vec u);
vec theta(vec x, vec y);
vec itheta(vec x, vec y);
vec delta(vec x, vec y);
vec w1(vec u);
vec w2(vec u);
vec w3(vec u);
vec k_1(vec u);
vec maxx(vec x, vec y);

vec evolve(vec u_0, double t_end);
vec noise(int N);
void measure_width(double t0, double t_end, double r, double rep);

int N = 32;
double dt = 0.001;
double t_end = 1e5;
double r = 1.01;
int rep = 1e2;
double t0 = 2;
double a = 1.0, tau = 1.0;
double reg = 1.0;

mt19937 rng(time(NULL));
normal_distribution<double> distribution(0.0,1.0);

int main() {

	arma_rng::set_seed_random();
	measure_width(t0, t_end, r, rep);

	return 0;
}

vec k_1(vec u) {
	return a / tau * (w1(u) + w2(u) + w3(u));
}

vec w1(vec u) {
	vec hp1 = shift(u,-1);
	vec hp2 = shift(u,-2);
	vec hm1 = shift(u,1);
	vec hm2 = shift(u,2);
	vec w = zeros<vec>(N);
	w += theta(hm1,hm2) % theta(hp1,hp2) % delta(hm1,u) % delta(u,hp1);
	w += itheta(u,hm1) % theta(u,hp1);
	w += theta(u,hm1) % itheta(u,hp1);
	w += itheta(u,hm1) % itheta(u,hp1);
	return w;
}

vec w2(vec u) {
	vec hp1 = shift(u,-1);
	vec hp2 = shift(u,-2);
	vec hp3 = shift(u,-3);
	vec hm1 = shift(u,1);
	vec w = zeros<vec>(N);
	w += delta(u,hp1) % delta(hp1,hp2) % itheta(u,hm1) % (theta(hp2,hp3) + 0.5*itheta(hp2,hp3));
	w += delta(hp1,hp2) % itheta(u,hp1) % (theta(hp2,hp3) + 0.5*itheta(hp2,hp3));
	w += 0.5*delta(hp1,u) % itheta(u,hm1) % itheta(hp2,hp1);
	w += 0.5*itheta(u,hp1) % itheta(hp2,hp1);
	return w;
}

vec w3(vec u) {
	vec hp1 = shift(u,-1);
	vec hm1 = shift(u,1);
	vec hm2 = shift(u,2);
	vec hm3 = shift(u,3);
	vec w = zeros<vec>(N);
	w += delta(hm2,hm1) % delta(hm1,u) % itheta(u,hp1) % (theta(hm2,hm3) + 0.5*itheta(hm2,hm3));
	w += delta(hm2,hm1) % itheta(u,hm1) % (theta(hm2,hm3) + 0.5*itheta(hm2,hm3));
	w += 0.5*delta(hm1,u) % itheta(hm2, hm1) % itheta(u,hp1);
	w += 0.5*itheta(hm2,hm1) % itheta(u,hm1);
	return w;
}

vec theta2(vec x, vec y) {
	vec t = zeros<vec>(N);
	for (int i = 0; i < N; i++){
		if (x[i] >= y[i]) t[i] = 1.0;
		else t[i] = 0.0;
	}
	return t;
}

vec theta(vec x, vec y) {
	vec w = x - y;
	vec t = zeros<vec>(N);
	for (int i = 0; i < N; i++) {
		if (w[i] < 0.0 && w[i] > -reg) t[i] = w[i]/reg + 1.0;
		else if (w[i] >= 0.0) t[i] = 1.0;
	}
	return t;
}

vec itheta(vec x, vec y) {
	return 1 - theta(x, y);
}

vec delta(vec x, vec y) {
	return theta(x,y) + theta(y,x) - 1;
}

// time derivative of field
vec u_t(vec u) {
	vec k1 = k_1(u);
	return k1 + sqrt(a*k1/dt) % noise(N);
}

vec evolve(vec u0, double t_end) {
	vec u = u0;
	vec k1, k2, k3, k4;
	for (double t = 0; t <= t_end; t += dt) {
		u += dt * u_t(u);
	}
	return u;
}

void measure_width(double t0, double t_end, double r, double rep) {
	double delta_t, t;
	vec u;
	int M = int(log(t_end/t0)/log(r)) + 1;
	double width[M];

	cout << "number of data points: " << M << endl;

	for (int i = 0; i < M; i++) {
		width[i] = 0;
	}
	ofstream data;
	data.open("data.txt");

	for (int n = 0; n < rep; n++) {
		u = zeros<vec>(N);
		t = t0;
		u = evolve(u, t0);
		for (int m = 0; m < M; m++) {
			width[m] += var(u);
			delta_t = (r-1)*t;
			u = evolve(u, delta_t);
			t += delta_t;
		}
		cout << "rep: " << n+1 << endl;
	}
	for (int i = 0; i < M; i++) {
		data << t0 * pow(r, i) << ", " << sqrt(width[i]/(rep)) << endl;
	}
	data.close();
}

vec noise(int N) {
	vec noise(N);
	for (int i = 0; i < N; i++){
		noise[i] = distribution(rng);
	}
	return noise;
}












