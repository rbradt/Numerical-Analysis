#include "Gauss.hpp"

int main() {
	int n = 5;
	
	vector<double> a(n-1, -1);
	vector<double> b(n, 40);
	vector<double> c(n-1, -1);
	vector<double> d(n, 4);
	vector<double> x(n, 0);
	
	*--b.end() = -20;
	*b.begin() = -20;
	
	Tri(a,b,c,d,x);
	
	for(int i = 0; i < n; i++)
		cout << "x(" << i+1 << ") = " << x[i] << endl;
	
	vector<vector<double>> A = {{0.4096, 0.1234, 0.3678, 0.2943},
							   {0.2246, 0.3872, 0.4015, 0.1129},
							   {0.3645, 0.1920, 0.3781, 0.0643},
							   {0.1784, 0.4002, 0.2786, 0.3927}};
							   
	vector<vector<double>> A2 = {{1,1,1},{1,0,1},{1,1,0}};
	vector<double> b2 = {1,2,3};
	vector<double> x2(3, 0);
							   
	vector<double> x1(4, 0);
	vector<double> b1 = {0.4043, 0.1550, 0.4240, 0.2557};
	
	SppForwardGauss(A,x1,b1);
	
	for(int i = 0; i < 4; i++)
		cout << "x(" << i+1 << ") = " << x1[i] << endl;
}