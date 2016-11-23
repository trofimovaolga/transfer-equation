#include <cstdio> 
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;				
const long double z0 = 0, z1 = 0.5, z2 = 1;	//границы

long double alpha() {
	long double a = ((long double)rand() / RAND_MAX);
	return a;
}

int L(long double z_new) {

	if ((z0 <= z_new) && (z_new < z1)) { 
		return 1; 
	}
	if ((z1 <= z_new) && (z_new <= z2)) { 
		return 2; 
	}
	return 0;
}

long double l(long double z_new) {
	if (L(z_new) == 1) { 
		return abs(z1 - z0); 
	}
	if (L(z_new) == 2) { 
		return abs(z2 - z1); 
	}
	return 0;
}

long double h(long double v, int L) {
	if (L = 1) {
		if (v > 0) { return 1; }
		else { return 0; }
	}
	if (L = 2) {
		if (v > 0) { return 0; }
		else { return 0.5; }
	}
	return 0;
}

long double lambda(int L) {
	if (L == 1) { return 0.9; }
	if (L == 2) { return 0.92; }
	return 0;
}

long double exponenta(long double z, long double v, int L) {
	long double expon = 0;
	long double d = 0;
	
	if (L == 1) {
		if (abs(v) > 0.00000000001) {
			if (v > 0) { d = z / v; }
			else { d = abs((z - l(z)) / v); }
			return expon = expl(-d);
		}
	}
	if (L == 2) {
		if (abs(v) > 0.00000000001) {
			if (v > 0) { 
				if (z == l(z)) { d = l(z) / v; }
				else { d = (z - l(z)) / v; }
			}
			else { d = abs((z - l(z)) / v); }
			return expon = expl(-d);
		}
	}
	return 0;
}

long double f_0(long double z, long double v, int L) {
	return h(v, L) * exponenta(z, v, L);
}

long double R() {
	long double R1, R2;
	R1 = 1; 
	R2 = 1;

	return (R1*R1 + R2*R2) / 2;
}

long double func(long double z_new, long double v_new, int n, int L) {
	if (n == 0) { 
		return f_0(z_new, v_new, L); 
	}
	long double z_old = z_new;
	long double v_old = v_new;
	long double z, v, e = exponenta(z_old, v_old, L);
	long double t = -logl(1 - alpha()*(1 - e));
	z = z_old - t*v_old;
	v = 2 * alpha() - 1;

	long double AS = (1 - e)*lambda(L);
	long double result = AS*func(z, v, n-1, L) + f_0(z_old, v_old, L); //+Bf*exp
	return result;
}

int main() {
	ofstream out;
	out.open("result.txt");
	long double AS = 0, z_new = 0, v_new = 0.5;
	int N = 10000;							//число траекторий
	int n = 20;								//номер точки на траектории

	const int m = 10;
	long double s[m];
	long double step = z2 / m;
	s[0] = z_new;

	for (int k = 1; k < m; k++) {
		s[k] = s[k - 1] + step;
		//cout << s[k] << endl;
	}

	for (int q = 1; q < m; q++) {
		srand(1);
		z_new = s[q];
	    cout << z_new << endl;
		int Lvar = L(s[q]); 
		cout << Lvar << endl;
		for (int i = 1; i <= N; i++) {
			long double f_var = func(s[q], v_new, n, Lvar);
			AS = (AS * (i-1) + f_var) / i;
		}
		long double f = AS;
		out << f << "\n";
	}
	system("pause");
	out.close();
	return 0;
}