#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
vector <vector<double>> mmultiply (vector <vector<double>> A, vector <vector<double>> B){ //matrix multiplication
    if (A[0].size() == B.size()){
        vector <vector<double>> C (A.size());
        for (int i = 0; i < C.size(); i++){
            C[i].resize(B[0].size());
        }
        for (int i = 0; i < A.size(); i++){
            for (int j = 0; j < B[0].size(); j++){
                for (int k = 0; k < A[0].size(); k++){
                    C[i][j] += A[i][k]*B[k][j];
                }
            }
        }
        return C;
    }
}
pair<vector<vector<double>>,vector<vector<double>>> LUP (vector <vector<double>> A){ //returns a matrix that stores L,U and a matrix P (transposition matrix)
    vector <vector<double>> P(A.size());
    for(int i = 0; i < P.size(); i++){
        P[i].resize(A.size(),0);
        P[i][i] = 1;
    }
    for(int i = 0; i<A.size()-1; i++){
        double lead = INT64_MIN;
        double nlead = -1;
        for (int j = i; j < A.size(); j++){
            if (abs(A[j][i]) > lead){
                lead = abs(A[j][i]);
                nlead = j;
            }
        }
        swap(A[i],A[nlead]);
        swap(P[i],P[nlead]);
        for (int j = i+1; j < A.size(); j++){
            A[j][i] = A[j][i]/A[i][i];
            for (int k = i+1; k<A.size(); k++){
                A[j][k] = A[j][k]-A[j][i]*A[i][k];
            }
        }
    }
    return make_pair(A,P);
}
vector <double> LUPsolve(vector <vector<double>> A, vector<vector<double>> b){ //solves the equation system by using the results of LUP function
    pair<vector<vector<double>>,vector<vector<double>>> LpUaP = LUP (A);
    vector<vector<double>> LU = LpUaP.first;
    b = mmultiply(LpUaP.second,b);
    vector<double> y(b.size());
    for(int i = 0; i<b.size(); i++){
        y[i] = b[i][0];
    }
    for (int i = 0; i < A.size(); i++){
        for (int k = 0; k<i;k++){
            y[i]-=LU[i][k]*y[k];
        }
    }
    vector<double> x(b.size());
    for(int i = b.size()-1; i>=0; i--){
        x[i] = y[i];
        for (int k = i+1; k<b.size(); k++){
            x[i] -= LU[i][k]*x[k];
        }
        x[i] = x[i]/LU[i][i];
    }
    return x;
}
double f(double x){
    return pow(x,2) + asin(x-0.2);
}
vector <vector<double>> points_eq(double (*f)(double), double a, double b, int n){
    vector <vector<double>> result(2);
    double step = (b-a)/(n);
    for (int i = 0; i <= n; i++){
        result[0].push_back(f(a));
        result[1].push_back(a),
        a+=step;
    }
    return result;
}
vector<vector<double>> points_op(double (*f)(double), double a, double b, int n){
    vector <vector<double>> result(2);
    double xi;
    for (int i = n; i >= 0; i--){
        xi = 0.5*((b-a)*cos(M_PI*(2*i+1)/(2*(n+1)))+b+a);
        result[0].push_back(f(xi));
        result[1].push_back(xi);
    }
    return result;
}
double lm(double x, int n, const vector<double>& table) {
    double result = 1;
    for (int i = 0; i < table.size(); i++) {
        if (i != n) {
            result *= (x - table[i]) / (table[n] - table[i]);
        }
    }
    return result;
}
double lp(double x, const vector<vector <double>> &table){
    int n = table[0].size();
    double result = 0;
    for (int i = 0; i < n; i++){
        result+=lm(x,i,table[1])*table[0][i];
    }
    return result;
}
vector<double> div_dif(const vector<vector<double>>& table) {
    int n = table[1].size()-1;
    vector<vector<double>> result(n);
    vector<double> coeffs = {table[0][0]}; // Initialize with function value at first point

    for (int i = 0; i < n; i++) {
        result[0].push_back((table[0][i+1]-table[0][i])/(table[1][i+1]-table[1][i]));
    }
    coeffs.push_back(result[0][0]);

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n-i; j++) {
            result[i].push_back((result[i-1][j+1] - result[i-1][j])/(table[1][j+i+1] - table[1][j]));
        }
        coeffs.push_back(result[i][0]); // Append new coefficient to coeffs vector
    }

    return coeffs; // Return vector of coefficients
}
// Calculates Newton's polynomial at point x for the given table
double np(double x, const vector<vector<double>>& table) {
    vector<double> coeffs = div_dif(table);
    double result = coeffs[0];

    for (int i = 1; i < coeffs.size(); i++) {
        double term = coeffs[i];
        for (int j = 0; j < i; j++) {
            term *= (x - table[1][j]);
        }
        result += term;
    }

    return result;
}
vector<vector<double>> lin_coeffs(const vector<vector<double>>& table){
    // Calculate number of data points
    int n = table[1].size();

    // Create vector b for solving the linear equation system
    vector<vector<double>> b ((n-1)*2);
    b[0].push_back(table[0][0]); // First point
    for (int i = 1; i <= (n-2); i++){
        b[2*i-1].push_back(table[0][i]); // Midpoints (y values)
        b[2*i].push_back(table[0][i]); // Midpoints (y values)
    }
    b[2*n-3].push_back(table[0][n-1]); // Last point

    // Create matrix A for solving the linear equation system
    vector<vector<double>> A ((n-1)*2);
    A[0].resize((n-1)*2,0); // Initialize first row of A with zeros
    A[0][0] = table[1][0]; // First point x value
    A[0][1] = 1; // Coefficient of constant term
    for (int i = 1; i <= (n-2); i++){
        A[2*i-1].resize((n-1)*2,0); // Initialize row for ith midpoint with zeros
        A[2*i-1][2*i-2] = table[1][i]; // Coefficient of linear term for ith midpoint
        A[2*i-1][2*i-1] = 1; // Coefficient of constant term for ith midpoint
        A[2*i].resize((n-1)*2,0); // Initialize row for (i+1)th midpoint with zeros
        A[2*i][2*i] = table[1][i]; // Coefficient of linear term for (i+1)th midpoint
        A[2*i][2*i+1] = 1; // Coefficient of constant term for (i+1)th midpoint
    }
    A[2*n-3].resize((n-1)*2,0); // Initialize last row of A with zeros
    A[2*n-3][2*n-4] = table[1][n-1]; // Last point x value
    A[2*n-3][2*n-3] = 1; // Coefficient of constant term for last point

    // Solve the linear equation system using LUP decomposition and store the coefficients in a vector
    vector <double> coeffs = LUPsolve(A,b);

    // Create a vector of vectors to store the polynomial coefficients for each interval
    vector<vector<double>> result(n-1);
    for (int i = 0; i < n-1; i++){
        result[i].push_back(coeffs[2*i]); // Coefficient of linear term for ith interval
        result[i].push_back(coeffs[2*i+1]); // Coefficient of constant term for ith interval
    }
    return result;
}
double lin_interp(double x,const vector<vector<double>>& table) {
    int n = table[0].size();
    if (x < table[1][0] || x > table[1][n-1]) {
        // Return an error code if x is outside the range of the table
        return NAN;
    }

    // Find the interval that x belongs to
    int i;
    for (i = 0; i < n-1; i++) {
        if (x >= table[1][i] && x <= table[1][i+1]) {
            break;
        }
    }

    // Calculate the coefficients for the linear spline
    vector<vector<double>> coeffs = lin_coeffs(table);

    // Evaluate the linear spline at x
    double y = coeffs[i][1] + coeffs[i][0] * x;

    return y;
}
vector<vector<double>> quad_coeffs(const vector<vector<double>>& table) {
    // Calculate number of data points
    int n = table[1].size();

    // Create vector b for solving the linear equation system
    vector<vector<double>> b((n - 1) * 3, vector<double>(1));
    b[0][0] = table[0][0]; // First point
    for (int i = 1; i <= (n - 2); i++) {
        b[3 * i - 2][0] = table[0][i]; // Midpoints (y values)
        b[3 * i - 1][0] = 0; // Midpoints (y values)
        b[3 * i][0] = table[0][i]; // Midpoints (y values)
    }
    b[3 * n - 5][0] = table[0][n - 1]; // Last point

// Create matrix A for solving the linear equation system
    vector<vector<double>> A((n - 1) * 3, vector<double>((n - 1) * 3, 0));
    A[0][0] = table[1][0] * table[1][0]; // Coefficient of quadratic term for 1st midpoint
    A[0][1] = table[1][0]; // Coefficient of linear term
    A[0][2] = 1; // Coefficient of constant term
    for (int i = 1; i <= (n - 2); i++) {
        A[3 * i - 2][3 * i - 3] = table[1][i] * table[1][i]; // Coefficient of quadratic term for ith midpoint
        A[3 * i - 2][3 * i - 2] = table[1][i]; // Coefficient of linear term for ith midpoint
        A[3 * i - 2][3 * i - 1] = 1; // Coefficient of constant term for ith midpoint
        A[3 * i - 1][3 * i - 3] = table[1][i] * 2; //Coefficients of equation for derivatives
        A[3 * i - 1][3 * i - 2] = 1;
        A[3 * i - 1][3 * i - 1] = 0;
        A[3 * i - 1][3 * i] = -table[1][i]*2;
        A[3 * i - 1][3 * i + 1] = -1;
        A[3 * i - 1][3 * i + 2] = 0;
        A[3 * i][3 * i] = table[1][i] * table[1][i]; // Coefficient of quadratic term for (i+1)th midpoint
        A[3 * i][3 * i + 1] = table[1][i]; // Coefficient of linear term for (i+1)th midpoint
        A[3 * i][3 * i + 2] = 1; // Coefficient of constant term for (i+1)th midpoint
    }
    A[3*n - 5][3*n - 6] = table[1][n - 1] * table[1][n - 1]; // Second to last midpoint x value squared
    A[3*n - 5][3*n - 5] = table[1][n - 1]; // Coefficient of linear term
    A[3*n - 5][3*n - 4] = 1; // Coefficient of constant term
    A[3*n - 4][3*n - 6] = 2 * table[1][n - 1]; //Coefficients of equation for last derivative
    A[3*n - 4][3*n - 5] = 1;
    // Solve the linear equation system Ax = b for x using LUP decomposition
   vector<double> x = LUPsolve(A, b);

// Reshape x into vector of vector of coefficients
    vector<vector<double>> coeffs(n - 1, vector<double>(3));
    for (int i = 0; i <= (n - 2); i++) {
        coeffs[i][0] = x[3 * i];
        coeffs[i][1] = x[3 * i + 1];
        coeffs[i][2] = x[3 * i + 2];
    }
    return coeffs;
}

double quad_interp(double x, const vector<vector<double>>& table) {
    int n = table[0].size();
    if (x < table[1][0] || x > table[1][n-1]) {
        // Return an error code if x is outside the range of the table
        return NAN;
    }

    // Find the interval that x belongs to
    int i;
    for (i = 0; i < n-1; i++) {
        if (x >= table[1][i] && x <= table[1][i+1]) {
            break;
        }
    }

    // Calculate the coefficients for the quadratic spline
    vector<vector<double>> coeffs = quad_coeffs(table);

    // Evaluate the quadratic spline at x
    double y = coeffs[i][2] + coeffs[i][1] * x + coeffs[i][0] * pow(x, 2);

    return y;
}

double cubic_interp ( double x, const vector<vector<double>>& table){
    // Calculate number of data points
    int n = table[1].size();
    vector<vector<double>> b(n-2,vector<double>(1));
    vector<double> h(n-1);
    for(int i = 0; i < n-1; i++){
        h[i] = table[1][i+1] - table[1][i];
    }
    for(int i = 1; i < n-1; i++){
        b[i-1][0] = 6*((table[0][i+1]-table[0][i])/h[i] - (table[0][i]-table[0][i-1])/h[i-1]);
    }
    vector<vector<double>> H(n-2, vector<double>(n-2, 0));
    H[0][0] = 2 * (h[0]+h[1]);
    H[0][1] = h[1];
    for (int i = 1; i < n-3; i++){
        H[i][i-1] = h[i];
        H[i][i] = 2 * (h[i]+h[i+1]);
        H[i][i+1] = h[i+1];
    }
    H[n-3][n-4] = h[n-3];
    H[n-3][n-3] = 2 * (h[n-3]+h[n-2]);
    vector <double> y_der2(n, 0), y_der1(n,0), solved = LUPsolve(H, b);
    for(int i = 1; i < n-1; i++){
        y_der2[i] = solved[i-1];
    }
    for(int i = 0; i < n-1; i++){
        y_der1[i] = (table[0][i+1] - table[0][i])/h[i] - y_der2[i+1]*h[i]/6 - y_der2[i]*h[i]/3;
    }
    if (x < table[1][0] || x > table[1][n-1]) {
        // Return an error code if x is outside the range of the table
        return NAN;
    }

    // Find the interval that x belongs to
    int i;
    for (i = 0; i < n-1; i++) {
        if (x >= table[1][i] && x <= table[1][i+1]) {
            break;
        }
    }
    //Evaluate the cubic spline at x
    double y = table[0][i] + y_der1[i]*(x-table[1][i]) + y_der2[i] * pow(x - table[1][i],2)/2 + (y_der2[i+1] - y_der2[i]) * pow(x - table[1][i],3)/(6 * h[i]);
    return y;
}


int main() {
    int m = 300, n;
    double a = -0.799999, b = 1.199999, interpolated, error, mx_eq = -1, mx_op = -1;
    vector<int> points_amount = {3,10,25,50};
    for (int k = 0; k < 4; k++) {
        n = points_amount[k];
        cout << endl << "Amount of points: " << n << endl;
        cout.precision(8);
        double step = (b - a) / (m);
        double x = a;
        vector<vector<double>> table = points_eq(&f, a, b, n);
        for (int i = 0; i <= m; i++) {
            interpolated = cubic_interp(x, table);
            error = abs(f(x) - interpolated);
            mx_eq = max(mx_eq, error);
            if (!isnan(interpolated)) {
                cout << fixed << x << ' ' << interpolated << endl;
                x += step;
            }
        }
        cout << endl << endl;
        x = a;
        table = points_op(&f, a, b, n);
        for (int i = 0; i <= m; i++) {
            interpolated = cubic_interp(x, table);
            error = abs(f(x) - interpolated);
            mx_op = max(mx_op, error);
            if (!isnan(interpolated)) {
                cout << fixed << x << ' ' << interpolated << endl;
            }
            x += step;
        }
        cout << endl << n << ' ' << m << ' ' << mx_eq << ' ' << mx_op << endl;
        mx_eq = -1;
        mx_op = -1;
    }
    return 0;
}
