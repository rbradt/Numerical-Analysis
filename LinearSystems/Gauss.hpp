#ifndef __GAUSS_HPP__
#define __GAUSS_HPP__

#include<vector>
#include<iostream>

using namespace std;

double max(const vector<double>& x) {
	double max = abs(x[0]);
	for(int i = 1; i < x.size(); ++i)
		if(abs(x[i]) > max)
			max = abs(x[i]);
		
	return max;
}

/********************************************************************************************
 *             Scaled Partial Pivoting Forward Gaussian Elimination Method
 ********************************************************************************************
 * Takes a nxn matrix A, and nx1 vectors x and b such that Ax=b. Solves for x.
 * - Solves linear systems of equations.
 ********************************************************************************************/
void SppForwardGauss(vector<vector<double>>& A, vector<double>& x, vector<double> b) {
	vector<int> rank(A.size(), 0);		// Index vector
	vector<double> scale(A.size(), 0);	// Scale vector
	
	// Initialize the index and scale vectors
	for(int i = 0; i < A.size(); ++i) {
		rank[i] = i;
		scale[i] = max(A[i]);
	}
	
	// iterate through all pivot points (except the last)
	for(int n = 0; n < A.size()-1; ++n) {	
		/*****************  Determine rank of rows using the scale vector  *****************/
		double rmax = abs(A[rank[n]][n])/scale[rank[n]]; // maximum scaled pivot
		int k = n;									// index of maximum scaled pivot
		
		// determine the maximum scaled pivot
		for(int i = n+1; i < A.size(); ++i) {
			double r = abs(A[rank[i]][n])/scale[rank[i]];
			if(r > rmax){
				rmax = r;
				k = i;
			}
		}
		
		// swap the index of the max scaled pivot and the index of the current first element
		int temp = rank[n];
		rank[n] = rank[k];
		rank[k] = temp;
		
		/*****************  Gaussian Elimination  *****************/
		// iterate through non-pivot rows
		for(int i = n+1; i < A.size(); ++i)	{	
			// calculate multiplier for this row
			double mult = A[rank[i]][n]/A[rank[n]][n];
			A[rank[i]][n] = mult;
			// eliminate the row
			// iterate through non-pivot columns
			for(int j = n+1; j < A.size(); ++j) {
				// eliminate column element
				A[rank[i]][j] -= mult*A[rank[n]][j];
			}
		}
	}
	
	cout << "Rank, Scale: " << endl;
		for(int i = 0; i<A.size(); ++i) {
			cout << rank[i] << ", ";
			cout << scale[i] << endl;
		}
	
	cout << "Matrix A:" << endl;
	for(int i =0; i<A.size(); ++i) {
		for(int j = 0; j<A.size(); ++j) {
				cout << A[i][j] << ", ";
		}
		cout << b[i] << endl;
	}
	
	/*****************  Backwards Substitution  *****************/
	// calculate values of the b vector
	for(int n=0; n < A.size()-1; ++n)
		for(int m = n+1; m < A.size(); ++m)
			b[rank[m]] -= A[rank[m]][n]*b[rank[n]];
	
	// backwards substitute
	x[A.size()-1] = b[rank[A.size()-1]]/A[A.size()-1][A.size()-1];
	// iterate through pivots in reverse skipping the last one
	for(int i = A.size()-2; i>-1; --i) { 
		// subtract values to the right of the pivot from b
		double s = b[rank[i]];
		for(int j = i + 1; j < A.size(); ++j)
			s -= A[rank[i]][j]*x[j];
		
		// calculate x by dividing the sum by the pivot
		x[i] = s/A[rank[i]][i]; 
	}
}

void Tri(const vector<double>& a, vector<double> b, vector<double> c, vector<double> d, vector<double>& x) {
	for(int i = 1; i < d.size(); ++i) {
		// calculate multiplier
		double m = a[i-1]/d[i-1];
		// eliminate the row past the pivot row
		d[i] -= m*c[i-1];
		// update b vector
		b[i] -= m*b[i-1];
	}
	
	// backwards substitution to solve for x
	x[d.size()-1] = b[d.size()-1]/d[d.size()-1];
	for(int i = d.size()-2; i > -1; --i)
		x[i] = (b[i] - c[i]*x[i+1])/d[i];
}

#endif