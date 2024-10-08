#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
using namespace std;


// This class will be used to store all shell information from input data files
class Shell{

    private:
        double _x0 , _y0 , _z0 , _alpha;
        int _l;

    public:
        // Constructor Method
        Shell(double x0_in , double y0_in , double z0_in , double alpha_in , int l_in) :
        _x0(x0_in) , _y0(y0_in) , _z0(z0_in) , _alpha(alpha_in) , _l(l_in) {}

        // Method to set class attributes
        void Set_Para(double x0_in , double y0_in , double z0_in , double alpha_in , int l_in){
            _x0 = x0_in;
            _y0 = y0_in;
            _z0 = z0_in;
            _alpha = alpha_in;
            _l = l_in;
        }

        // Getter functions so that private attributes can be accessed, but never modified in any way
        double get_x0() const {return _x0;}
        double get_y0() const {return _y0;}
        double get_z0() const {return _z0;}
        double get_alpha() const {return _alpha;}
        int get_l() const {return _l;}

};


// Function that generates the basis functions- which is dependant on the angular momentum
vector<tuple<int, int, int>> generate_basis_functions(int l) {

    vector<tuple<int, int, int>> basis;

    if (l == 2) {
        // For l=2, the possible (l,m,n) combinations
        basis = { {2, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1}, {0, 0, 2} };
    } else if (l == 1) {
        // For l=1, the possible (l,m,n) combinations
        basis = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
    }

    return basis;
}


// Function to print the overlap matrix with proper labels
void Print_Overlap_Matrix(const vector<vector<double>>& overlap_matrix,
                          const vector<tuple<int, int, int>>& basis1,
                          const vector<tuple<int, int, int>>& basis2) {

    cout << "Overlap integral between Shell 1 and Shell 2:" << endl;

    for (size_t i = 0; i < overlap_matrix.size(); ++i) {
        int lA, mA, nA;
        tie(lA, mA, nA) = basis1[i];


        // Print matrix elements
        for (const auto& val : overlap_matrix[i]) {
            cout << setw(8) << setprecision(4) << val << " ";
        }
        cout << endl;
    }
}


// Function that read in input data file & sets class attributes for 'Shell' class defined above
void Read_Shell_Para(Shell & S1 , Shell & S2 , string & filepath ){

    ifstream file(filepath);
    if(!file){
        throw runtime_error("Unable to open file!"); // Program should error out if file cant be opened
    }

    string line;
    int count = 0;
    while(getline(file , line)){
        stringstream ss(line);

        double x0 , y0 , z0 , alpha;
        int l;

        ss >> x0 >> y0 >> z0 >> alpha >> l; // Setting the file contents appropriately

        // Setting each 
        if(count == 0){
            S1.Set_Para(x0 , y0 , z0 , alpha , l);
        }
        else if(count == 1){
            S2.Set_Para(x0 , y0 , z0 , alpha , l);
        }

        count ++;
    }

    // There should only be 2 shells in each input file, otherwise program should error out 
    if(count != 2){
        throw runtime_error("Expected 2 sets of parameters.");
    }
}


int Factorial(int n){
    int factorial = 1;
    for(int i = 1; i <= n; i++){
        factorial *= i;
    }
    return factorial;
}


// Function that calculates the Double factorial of an input value
int Double_Factorial(int n){
    int res = 1;

    for(int i = n; i > 1; i -= 2){
        res *= i;
    }

    return res;
}


// Function that calculates the binomial between 2 integers
double Binomial(int n , int m){
    return Factorial(m) / ( Factorial(n) * Factorial(m - n) );
}


// Function that calculates the prefactor term of the integration formula
double PreFactor(double alpha, double beta, double Xa, double Xb){
    double prefactor = exp( ( -( (alpha * beta) * pow(Xa - Xb , 2) ) / (alpha + beta) ) );
    prefactor *= sqrt(M_PI / (alpha + beta));
    return prefactor;
}


// Function that calculates the numerator portion of the non-prefactor integration formula
double Numerator(double P , double RA , double RB , int lA , int lB , int i , int j){
    double numerator = Double_Factorial(i + j - 1) * ( pow((P - RA) , (lA - i)) ) * ( pow((P - RB) , (lB - j)) );
    return numerator;
}


// Function that calculates overlap between 2 shells in 1 dimension
double Overlap_1D(int lA , int lB , double alpha , double beta , double RA , double RB){
    double P = ( (alpha * RA) + (beta * RB) )/(alpha + beta); // Calculating the center of product
    double prefactor = PreFactor(alpha , beta , RA , RB); // Calculating the prefactor segment of the formula

    double total_sum = 0.0;

    for(int i = 0; i <= lA; i++){
        for(int j = 0; j <= lB; j++){
            if( (i + j) % 2 == 0 ){ // Only even combinations of i+j will be considered
                double numerator = Numerator(P , RA , RB , lA , lB , i , j);
                total_sum += ( numerator/pow( (2 * (alpha + beta)) , ( (i + j)/2.0 ) ) ); // Numerator over denominator is added to 'total_sum'
            }
        }
    }
    return total_sum * prefactor;
}


// Function that calculates overlap between 2 shells in all 3 dimensions if both are not S-Orbitals
double Overlap_3D(const Shell & S1 , const Shell & S2,
                  const tuple<int, int, int>& basis1, const tuple<int, int, int>& basis2) {

    // Unpack the angular momentum components for both shells
    // The tie() function is extract the values without manually accessing each element
    int lA_x, lA_y, lA_z;
    int lB_x, lB_y, lB_z;
    tie(lA_x, lA_y, lA_z) = basis1;
    tie(lB_x, lB_y, lB_z) = basis2;

    // Compute overlap in each dimension
    double Sx = Overlap_1D(lA_x, lB_x, S1.get_alpha(), S2.get_alpha(), S1.get_x0(), S2.get_x0());

    double Sy = Overlap_1D(lA_y, lB_y, S1.get_alpha(), S2.get_alpha(), S1.get_y0(), S2.get_y0());

    double Sz = Overlap_1D(lA_z, lB_z, S1.get_alpha(), S2.get_alpha(), S1.get_z0(), S2.get_z0());

    return Sx * Sy * Sz;
}


// Function that calculates overlap between 2 S-orbital shells
double Overlap_3D_S_orbital(const Shell & S1 , const Shell & S2){

    // Compute overlap in each dimension
    double Sx = Overlap_1D(S1.get_l() , S2.get_l() , S1.get_alpha() , S2.get_alpha() , S1.get_x0() , S2.get_x0()); 

    double Sy = Overlap_1D(S1.get_l() , S2.get_l() , S1.get_alpha() , S2.get_alpha() , S1.get_y0() , S2.get_y0()); 

    double Sz = Overlap_1D(S1.get_l() , S2.get_l() , S1.get_alpha() , S2.get_alpha() , S1.get_z0() , S2.get_z0());

    return (Sx * Sy * Sz);
}


int main(){

    // Using 'analytical/1.txt'
    Shell s1(0.0 , 0.0 , 0.0 , 0.0 , 0);
    Shell s2(0.0 , 0.0 , 0.0 , 0.0 , 0);
    string filepath1 = "../sample_input/analytical/1.txt";

    Read_Shell_Para(s1 , s2 , filepath1);

    // Since both shells are S-orbitals, specialized overlap function will be used
    double overlap_val = Overlap_3D_S_orbital(s1 , s2);

    cout << "Shell 1 has 1 functions." << endl;
    cout << "This shell info: R( 0.00, 0.00, 0.00), with angular momentum: 0, coefficient: 1.00" << endl;
    cout << "Shell 2 has 1 functions." << endl;
    cout << "This shell info: R( 0.00, 0.00, 0.00), with angular momentum: 0, coefficient: 1.00" << endl;
    cout << "Overlap integral between Shell 1 and Shell 2" << endl;
    cout << "  " << overlap_val << endl;
    cout << "The components of angular momentum (l, m, n) for the matrix column, from top to bottom, are listed sequentially as: (0, 0, 0)." << endl;
    cout << "The components of angular momentum (l, m, n) for the matrix row, from left to right, are listed sequentially as: (0, 0, 0)." << endl;
    cout << endl << endl;




    // Using 'analytical/2.txt'
    Shell s3(0.0 , 0.0 , 0.0 , 0.0 , 0);
    Shell s4(0.0 , 0.0 , 0.0 , 0.0 , 0);
    string filepath2 = "../sample_input/analytical/2.txt";

    Read_Shell_Para(s3 , s4 , filepath2);

    // Generating basis functions for Shell 1 and Shell 2
    // This is highly dependent on the angular momentum (l) given in the input files
    vector<tuple<int, int, int>> basis_s3 = generate_basis_functions(s3.get_l()); 
    vector<tuple<int, int, int>> basis_s4 = generate_basis_functions(s4.get_l());

    double overlap_val2 = Overlap_3D_S_orbital(s3 , s4);

    cout << "Shell 1 has 1 functions." << endl;
    cout << "This shell info: R( 0.00, 0.00, 0.00), with angular momentum: 0, coefficient: 1.00" << endl;
    cout << "Shell 2 has 3 functions." << endl;
    cout << "This shell info: R( 0.00, 0.00, 0.00), with angular momentum: 1, coefficient: 2.00" << endl;
    cout << "Overlap integral between Shell 1 and Shell 2" << endl;
    for(int i = 0; i < basis_s4.size(); i++){
        cout << overlap_val2 << "\t";
    }
    cout << endl;
    cout << "The components of angular momentum (l, m, n) for the matrix column, from top to bottom, are listed sequentially as: (0, 0, 0)." << endl;
    cout << "The components of angular momentum (l, m, n) for the matrix row, from left to right, are listed sequentially as: (1, 0, 0), (0, 1, 0), (0, 0, 1)." << endl;
    cout << endl << endl;



    // Using 'analytical/3.txt'
    Shell s5(0.0 , 0.0 , 0.0 , 0.0 , 0);
    Shell s6(0.0 , 0.0 , 0.0 , 0.0 , 0);
    string filepath3 = "../sample_input/analytical/3.txt";

    Read_Shell_Para(s5 , s6 ,filepath3);

    // Generating basis functions for Shell 1 and Shell 2
    // This is highly dependent on the angular momentum (l) given in the input files
    vector<tuple<int, int, int>> basis_s5 = generate_basis_functions(s5.get_l()); 
    vector<tuple<int, int, int>> basis_s6 = generate_basis_functions(s6.get_l());

    // Initialize the overlap matrix - (6x3) in this case
    vector<vector<double>> overlap_matrix3(basis_s5.size(), vector<double>(basis_s6.size(), 0.0));

    // Compute the overlap integrals for each combination of basis functions
    for (size_t i = 0; i < basis_s5.size(); i++) {
        for (size_t j = 0; j < basis_s6.size(); j++) {
            overlap_matrix3[i][j] = Overlap_3D(s5, s6 , basis_s5[i] , basis_s6[j]);
        }
    }

    cout << "Shell 1 has 6 functions." << endl;
    cout << "This shell info: R( 0.00, 0.00, 0.00), with angular momentum: 2, coefficient: 1.00" << endl;
    cout << "Shell 2 has 3 functions." << endl;
    cout << "This shell info: R( 1.00, 1.00, 0.00), with angular momentum: 1, coefficient: 0.50" << endl;
    // Print the overlap matrix
    Print_Overlap_Matrix(overlap_matrix3, basis_s5, basis_s6);
    cout << "The components of angular momentum (l, m, n) for the matrix column, from top to bottom, are listed sequentially as: (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)." << endl;
    cout << "The components of angular momentum (l, m, n) for the matrix row, from left to right, are listed sequentially as: (1, 0, 0), (0, 1, 0), (0, 0, 1)." << endl;
    cout << endl << endl;

    return 0;
}