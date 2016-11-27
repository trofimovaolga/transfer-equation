#include <cstdio> 
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;
//mt19937 r;
void halt(int s) {
	//r();
	system("pause");
	exit(s);
}

const double k1 = 1.54, k2 = 1;				//коэффициенты преломления стекла и воздуха
const long double z0 = 0, z1 = 0.7, z2 = 3.5;	//границы
double eps = 1e-8;

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

bool is_eq(long double a, long double b) {
	return abs(a - b) < eps;
}

int L(long double z) {
	if ((z0 <= z) && (z < z1)) {
		return 1;
	}
	if ((z1 <= z) && (z <= z2)) {
		return 2;
	}
	cout << "mistake in L" << endl;
	halt(0);
	return 0;
}

double k(long double v) {
	return v > 0 ? k2 / k1 : k1 / k2;
}

double psi(long double v) {
	long double p = 1 - k(v)*k(v)*(1 - v*v);
	return (p >= 0) ? sgn(v)*sqrt(p) : 0;
}

long double R(long double v) {
	long double R1 = (k(v)*psi(v) - v) / (k(v)*psi(v) + v);
	long double R2 = (psi(v) - k(v)*v) / (psi(v) + k(v)*v);

	return (R1*R1 + R2*R2) / 2;
}

double ksi(long double z, long double v) {
	int Lvar = L(z);
	if (Lvar == 1){
		if (v > 0) { 
			return z0; 
		} 
		else { 
			return z1;
		}
	}
	if (Lvar == 2) {
		if (z == z1) {
			if (v > 0) {
				return z0;
			}
			else {
				return z2;
			}
		}
		if (v > 0) { 
			return z1; 
		}
		else { 
			return z2; 
		}
	}
	cout << "mistake in ksi" << endl;
	halt(0);
	return 0;
}

long double alpha() {
	double k = (double)rand() + 10, d = ((double)RAND_MAX + 20);
	long double a = k / d;
	return a;
}

long double l(long double z) {
	int Lvar = L(z);
	if (Lvar == 1) { 
		return abs(z1 - z0); 
	}
	if (Lvar == 2) { 
		return abs(z2 - z1); 
	}
	cout << "mistake in l" << endl;
	halt(0);
	return 0;
}

long double h(long double z, long double v) {
	long double k = ksi(z, v);
	if (k == z0) {
		return 1;
	}
	else if (k == z2) {
		return 0;
	}
	else { return 0;  }
	cout << "mistake in h" << endl;
	halt(0);
}

long double lambda(int L) {
	if (L == 1) { 
		return 0.22/0.226; 
	}
	if (L == 2) { 
		return 0.156/0.162; 
	}
	cout << "mistake in lambda" << endl;
	halt(0);
	return 0;
}

long double exponenta(long double z, long double v) {
	long double expon = 0.0;
	long double d = 0.0;
	int Lvar = L(z);

	if (Lvar == 1) {
		if (abs(v) > eps) {
			if (v > 0) { 
				d = z / v; 
			}
			else { 
				d = abs((z - z1) / v); 
			}
			expon = expl(-d);
		}
	}
	if (Lvar == 2) {
		if (abs(v) > eps) {
			if (v > 0) {
				if (is_eq(z, z1)) { 
					d = z1 / v; 
				}
				else { 
					d = (z - z1) / v; 
				}
			}
			else { 
				d = (z - z2) / v; 
			}
			expon = expl(-d);
		}
	}
	if (isinf(expon)) {
		cout << "exponenta is inf, d = " << d << ", e = " << expon << endl;
		halt(0);
	}
	return expon;
}

long double f_0(long double z, long double v) {
	return h(z, v) * exponenta(z, v);
}

long double func(long double z, long double v, int n) {
	if (n == 0) { 
		return f_0(z, v); 
	}
	long double z_new, v_new, e = exponenta(z, v);
	long double log_arg = 1 - alpha()*(1 - e);
	if (log_arg <= 0) {
		cout << "bad log argument" << endl;
		halt(0);
	}
	long double t = -logl(log_arg);
	
	z_new = z - t*v; 
	if (abs(z) >= z2) { 
		cout << "bad z variable" << endl;
		halt(0);
	}
	v_new = 2 * alpha() - 1;

	long double AS = (1 - e)*lambda(L(z));
	if (is_eq(ksi(z, v), z1)) {
		long double r = R(v), p = psi(v);
		long double B = r*func(z1, -v, n - 1) + (1 - r)*func(z1, p, n - 1);
		return B*e + AS*func(z_new, v_new, n - 1) + f_0(z, v);
	}
	return AS*func(z_new, v_new, n - 1) + f_0(z, v);
}

int main() {
	srand((unsigned)time(0));  
	ofstream out;
	out.open("result.txt");
	long double AS = 0, z_new = 0, v_new = 0.5;
	int N = 8000;							//число траекторий
	int n = 8;								//номер точки на траектории

	const int m = 20;
	long double s[m];
	long double step = z2 / m;
	s[0] = z_new;

	for (int k = 1; k < m; k++) {
		s[k] = s[k - 1] + step;
		//cout << s[k] << endl;
	}

	for (int q = 1; q < m; q++) {	//z_new = s[q];
		for (int i = 1; i <= N; i++) {
			AS = (AS * (i - 1) + func(s[q], v_new, n)) / i;
		}
		out << AS << "\n";
	}
	system("pause");
	out.close();
	return 0;
}

//макс 10 слоев сделать массивы