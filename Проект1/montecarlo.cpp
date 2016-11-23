#include <cstdio> 
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;	

const double k1 = 1.5, k2 = 1;				//ęîýôôčöčĺíňű ďđĺëîěëĺíč˙ ńňĺęëŕ č âîçäóőŕ
const long double z0 = 0, z1 = 0.5, z2 = 1;	//ăđŕíčöű

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

int L(long double z_new) {

	if ((z0 <= z_new) && (z_new < z1)) {
		return 1;
	}
	if ((z1 <= z_new) && (z_new <= z2)) {
		return 2;
	}
	//cout << "mistake in L" << endl;
	return 0;
}

double k(long double v) {
	if ((0 < v) && (v <= 1)) { return k2 / k1; }
	else 
		if ((-1 <= v) && (v < 0)) { return k1 / k2; }
	//cout << "mistake in k" << endl;
	return 0;
}

double psi(long double v) {
	long double p = 1 - k(v)*k(v)*(1 - v*v);
	if (p >= 0) { return sgn(v)*sqrt(p); }
	return 0; 
}

long double R(long double v) {
	long double R1 = (k(v)*psi(v) - v) / (k(v)*psi(v) + v);
	long double R2 = (psi(v) - k(v)*v) / (psi(v) + k(v)*v);

	return (R1*R1 + R2*R2) / 2;
}

double ksi(long double z, long double v) {
	if (L(z) == 1) {
		if (v > 0) { return z0; } 
		else { return z1;}
	}
	if (L(z) == 2) {
		if (v > 0) { return z1; }
		else { return z2; }
	}
	//cout << "mistake in ksi" << endl;
	return 0;
}

long double alpha() {
	long double a = ((long double)rand() / RAND_MAX);
	return a;
}

long double l(long double z_new) {
	if (L(z_new) == 1) { 
		return abs(z1 - z0); 
	}
	if (L(z_new) == 2) { 
		return abs(z2 - z1); 
	}
	//cout << "mistake in l" << endl;
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
	long double d = 0;
	
	if (L == 1) {
		if (abs(v) > 0.00000000001) {
			if (v > 0) { d = z / v; }
			else { d = abs((z - l(z)) / v); }
			return expl(-d);
		}
	}
	if (L == 2) {
		if (abs(v) > 0.00000000001) {
			if (v > 0) { 
				if (z == l(z)) { d = l(z) / v; }
				else { d = (z - l(z)) / v; }
			}
			else { d = abs((z - l(z)) / v); }
			return expl(-d);
		}
	}
	//cout << "mistake in exp" << endl;
	return 0;
}

long double f_0(long double z, long double v, int L) {
	return h(v, L) * exponenta(z, v, L);
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
	
	if (ksi(z, v) == z1) {
		long double r = R(v);
		long double B = r*func(ksi(z, v), -v, n-1, L) + (1 - r)*func(ksi(z, v), psi(v), n-1, L);
		return B*e + AS*func(z, v, n-1, L) + f_0(z_old, v_old, L); 
	}
	return AS*func(z, v, n-1, L) + f_0(z_old, v_old, L);
}

int main() {
	ofstream out;
	out.open("result.txt");
	long double AS = 0, z_new = 0, v_new = 0.5;
	int N = 10000;							//÷čńëî ňđŕĺęňîđčé
	int n = 20;								//íîěĺđ ňî÷ęč íŕ ňđŕĺęňîđčč

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
		int Lvar = L(s[q]); 
		for (int i = 1; i <= N; i++) {
			long double f_var = func(s[q], v_new, n, Lvar);
			AS = (AS * (i-1) + f_var) / i;
		}
		long double f = AS;
		out << f << "\n";
	}
	//system("pause");
	out.close();
	return 0;
}
