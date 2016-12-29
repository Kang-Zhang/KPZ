#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include <array>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace std;

//define parameters
const int N = 64;
double Dx = 1;
double Dy = 1;
double lx = 0;
double ly = 0;

random_device rd;
mt19937 gen(rd());

int periodic(int i, int N){
    return (i % N + N) % N;
}
typedef array <array<double, N> ,N> matrix; //shortcut for calling a data type

matrix change_lattice(matrix& L_new, matrix& L, double c_L, int N, double Dx,
                      double Dy, double lx, double ly, double dt){ //reference to avoid wasting memory
    uniform_real_distribution<> dis(-0.5,0.5);
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            //top = [periodic(i-1,N)][j];
            //right = [i][periodic(j+1,N)];
            //bottom = [periodic(i+1,N)][j];
            //left = [i][periodic(j-1,N)];

            double diffX = Dx*(sin(L[i][j] - L[i][periodic(j+1,N)]) + sin(L[i][j] - L[i][periodic(j-1,N)]));
            double diffY = Dy*(sin(L[i][j] - L[periodic(i-1,N)][j]) + sin(L[i][j] - L[periodic(i+1,N)][j]));
            double nonlinX = lx*( cos(L[i][j] - L[i][periodic(j+1,N)]) + cos(L[i][j] - L[i][periodic(j-1,N)]) - 1);
            double nonlinY = ly*(cos(L[i][j] - L[periodic(i-1,N)][j]) + cos(L[i][j] - L[periodic(i+1,N)][j]) - 1);
            double noise = 2*M_PI*c_L*dis(gen);

            //Euler update for KPZ equation
            double L_step = dt*(diffX + diffY + nonlinX + nonlinY + noise);
            L_new[i][j] = fmod((L[i][j] - L_step),2.0*M_PI);
            if (L_new[i][j] < 0){
                L_new[i][j] += 2.0*M_PI;}
        }
    }
    return L_new;
}

double energy (matrix& L_new, int N, double Dx, double Dy){
    double E = 0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
             E += -Dx*(cos(L_new[i][j] - L_new[i][periodic(j+1,N)]) + cos(L_new[i][j] - L_new[i][periodic(j-1,N)]))
                - Dy*(cos(L_new[i][j] - L_new[periodic(i-1,N)][j]) + cos(L_new[i][j] - L_new[periodic(i+1,N)][j]));
        }}
    return E;
}


double vortices(matrix& L_new, int N){
    double num = 0;
    double d[4];
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            //top = [periodic(i-1,N)][j]
            //top_left = [periodic(i-1,N)][periodic(j-1,N)]
            //left = [i][periodic(j-1,N)]

            d[0] = L_new[periodic(i-1,N)][j] - L_new[i][j];
            d[1] = L_new[periodic(i-1,N)][periodic(j-1,N)] - L_new[periodic(i-1,N)][j];
            d[2] = L_new[i][periodic(j-1,N)] - L_new[periodic(i-1,N)][periodic(j-1,N)];
            d[3] = L_new[i][j] - L_new[i][periodic(j-1,N)];

            for(int k=0; k < 4; k++){
                if(d[k] > M_PI){
                    d[k] -= 2*M_PI;}
                else if(d[k] < -M_PI){
                    d[k] += 2*M_PI;}}

            num += fabs(d[0]+d[1]+d[2]+d[3]);
        }}
    return num/2.0*M_PI;
}

void writeEV(int max_iters, int R, double c_L, int N, double Dx, double Dy, double lx, double ly,
             double dt, bool writeE, bool writeV, bool converge_test, double toleranceE)
{
    stringstream cL;
    cL << fixed << setprecision(2) << c_L;
    stringstream l_x;
    l_x << fixed << setprecision(2) << lx;
    stringstream l_y;
    l_y << fixed << setprecision(2) << ly;
    stringstream Lsize;
    Lsize << N;

    ofstream write_outputE;
    ofstream write_outputV;
    ofstream write_outputE_t;
    ofstream write_outputV_t;
    if (writeE){
        if(converge_test){write_outputE_t.open("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "E_t.txt");}
        else{write_outputE.open("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "E.txt");}
    }
    if (writeV){
        if(converge_test){write_outputV_t.open("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "V_t.txt");}
        else{write_outputV.open("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "V.txt");}
    }

    for (int r = 0; r<R;r++){
        //initialise random matrix for starting lattice configuration
        uniform_real_distribution<> dis(0,2*M_PI);
        matrix L;
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                    L[i][j] = dis(gen);}}

        matrix L_new; //initialise random matrix for storing new lattice
        //cout << r << "\n";

        //start of time evolution of lattice
        //add convergence test for energy and number of vortices
        int t = 0;
        double oldEdens = 0;
        double oldV = 0;

        while (t < max_iters){
            t += 1; //increase time step (iterations) by 1
            L = change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt);

            if (writeE){double E = energy(L, N, Dx, Dy);
                if (converge_test){
                    if (fabs((E/(N*N)) - oldEdens) < toleranceE){
                        write_outputE_t << E << " "<< t << "\n";
                        break;}
                    else{oldEdens = E/(N*N);}}
                else{write_outputE << E/(N*N) << " ";}}

            if (writeV){double V = vortices(L, N);
                if (converge_test){
                    if (fabs(V - oldV) == 0){
                        int t_i = t;
                        //need to compare next 30 V values to see if V is still the same
                        double Vtot = 0;
                        for (t = t_i; t < (t_i+30);t++){
                            Vtot += V;
                            oldV = V;}
                        if (Vtot/30 == V){
                            write_outputV_t << V << " "<< t << "\n";
                            break;}}
                    else{oldV = V;}}
                else{write_outputV << V << " ";}}
        }
        //end of one realisation so need new line after each realisation
        if (writeE && converge_test == false){write_outputE << "\n";}
        if (writeV && converge_test == false){write_outputV << "\n";}
        }
    //end of all realisations, need to close all files
    if (writeE){
        if (converge_test){write_outputE_t.close();}
        else{write_outputE.close();}}
    if (writeV){
        if (converge_test){write_outputV_t.close();}
        else{write_outputV.close();}}
}

int main(){
    //initial configuration of NxN lattice, each site with initial phase between 0 and 2pi chosen using random number generator

    //define parameters
    double c_L = 0;
    double dt = 0.05;
    const int max_iters = 3000;
    int R = 20; //number of realisations
    bool writeE = true;
    bool writeV = false;
    bool converge_test = true;
    double toleranceE = 1e-4;

    for (double c_L=0; c_L<=7; c_L+=0.5){
    writeEV(max_iters, R, c_L, N, Dx, Dy, lx, ly, dt, writeE, writeV, converge_test, toleranceE);
    }


}
