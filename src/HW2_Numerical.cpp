#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;


class Integrand{

    public:
        virtual double eval(double x) = 0;

        // Setting getter functions since these attributes will be private below
        virtual double get_x0() const = 0;
        virtual double get_alpha() const = 0;
        virtual double get_l() const = 0;

};

class Gaussian : public Integrand{

    private:
        double _x0 , _alpha;
        int _l;

    public:

        // Constructor method
        Gaussian(double x0_in , double alpha_in , int l_in)
            : _x0(x0_in) , _alpha(alpha_in) , _l(l_in){}

        // Method where Parameters can be set
        void Set_Para(double x0_in , double alpha_in , int l_in){
            _x0 = x0_in;
            _alpha = alpha_in;
            _l = l_in;
        }

        double eval(double x) override{
            double GA = pow((x - _x0) , _l) * exp(-_alpha * pow((x - _x0) , 2));
            return GA;
        }

        // Method that returns the value of x0
        double get_x0() const override {return _x0;}  

        // Method that returns the value of alpha
        double get_alpha() const override {return _alpha;}

        // Method that returns the angular momentum value
        double get_l() const override {return _l;}      

};

void Read_Gaussian_Para(Gaussian & G1 , Gaussian & G2 , string & filepath){

    ifstream file(filepath);
    if(!file){
        throw runtime_error("Unable to open file!");
    }

    string line;
    int count = 0;
    while(getline(file , line)){
        stringstream ss(line);

        double x0 , alpha;
        int l;

        ss >> x0 >> alpha >> l;

        if(count == 0){
            G1.Set_Para(x0 , alpha , l); // Setting values for first Gaussian obj
        }
        else if(count == 1){
            G2.Set_Para(x0 , alpha , l); // Setting values for second Gaussian obj
        }

        count ++;
    }

    if(count != 2){
        throw runtime_error("Expected 2 sets of parameters."); // There should only be 2 items in each file. Otherwise program should error out 
    }

}

void Trapo(Integrand & F1 , Integrand & F2 , const int num_points , const double stepsize){
    // Implement trapazoid rule
    int a = F1.get_x0() - 5; // Setting lower bound
    int b = F2.get_x0() + 5; // Setting upper bound

    int h = static_cast<int>(b - a) / num_points; // Calculating h

    double integral_val = 0.0;

    for(int i = 0; i < num_points; i++){
        double x = (i * stepsize) + a; // Current x val

        double eval_point = F1.eval(x) * F2.eval(x);

        if(i == 0 || i == h){
            integral_val += (eval_point / 2.0); // First & last point get halved
        }
        else{
            integral_val += eval_point;
        }
    }

    double integral = integral_val * stepsize;

    cout << "1d numerical overlap integral between Gaussian functions is " << scientific << setprecision(17) << integral << endl;

    cout.unsetf(ios::scientific);
}


int main(){

    // Using 'numerical/1.txt'
    Gaussian i1(0.0 , 0.0 , 0.0);
    Gaussian i2(0.0 , 0.0 , 0.0);
    string file_path1 = "../sample_input/numerical/1.txt";

    Read_Gaussian_Para(i1 , i2 , file_path1);

    Trapo(i1 , i2 , 1000000 , 0.5);

    cout << endl << endl;



    // Using 'numerical/2.txt'
    Gaussian i3(0.0 , 0.0 , 0.0);
    Gaussian i4(0.0 , 0.0 , 0.0);
    string file_path2 = "../sample_input/numerical/2.txt";

    Read_Gaussian_Para(i3 , i4 , file_path2);

    Trapo(i3 , i4 , 1000000 , 0.5);

    cout << endl << endl;



    // Using 'numerical/3.txt'
    Gaussian i5(0.0 , 0.0 , 0.0);
    Gaussian i6(0.0 , 0.0 , 0.0);
    string file_path3 = "../sample_input/numerical/3.txt";

    Read_Gaussian_Para(i5 , i6 , file_path3);

    Trapo(i5 , i6 ,  1000000 , 0.5);

    return 0;
}