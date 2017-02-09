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
const int N = 20;
double Dx = 1;
double Dy = 1;

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
            double nonlinX = (lx*0.5)*( cos(L[i][j] - L[i][periodic(j+1,N)]) + cos(L[i][j] - L[i][periodic(j-1,N)]) - 1);
            double nonlinY = (ly*0.5)*(cos(L[i][j] - L[periodic(i-1,N)][j]) + cos(L[i][j] - L[periodic(i+1,N)][j]) - 1);
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
        cout << c_L << "\n";
}

void writequench(double cL_init, double cL_final, int init_iter, int final_iter, int Gamma_Q, int N, double Dx,
                 double Dy, double lx, double ly, double dt)
{
    stringstream cL_i;
    cL_i << fixed << setprecision(2) << cL_init;
    stringstream cL_f;
    cL_f << fixed << setprecision(2) << cL_final;
    stringstream l_x;
    l_x << fixed << setprecision(2) << lx;
    stringstream l_y;
    l_y << fixed << setprecision(2) << ly;
    stringstream Lsize;
    Lsize << N;
    stringstream GammaQ;
    GammaQ << Gamma_Q;

    ofstream write_quenchV;
    write_quenchV.open("KPZProjectGraphs/" + Lsize.str() + "-" + cL_i.str() + "-" + cL_f.str() + "-" + GammaQ.str() + "-" + l_x.str() + "-" + l_y.str() + "V.txt");

    //obtain initial matrix and initial number of vortices for initial cL value
    matrix L_init; //to store initial lattice in quenching
    matrix L; //to initialise random matrix for configuring initial lattice
    uniform_real_distribution<> dis(0,2*M_PI);
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
                L[i][j] = dis(gen);}}
    matrix L_new; //initialise random matrix for storing new lattice
    for(int i = 0; i < init_iter; i++){
        L = change_lattice(L_new, L, cL_init, N, Dx, Dy, lx, ly, dt);}
    //steady state reached
    L_init = L;
    double V_init = vortices(L_init, N); //initial number of vortices at t = 0
    write_quenchV << V_init << " ";

    //finite quenching:
    matrix L_current; //to store current lattice in time (current iteration)
    matrix L_temp; //initialise temp matrix for storing lattice in time
    L_current = L_init; //at t = 0, current lattice = initial lattice
    double cL_quench = cL_init;

    for(int t = 1; t <= final_iter; t++){
        if(t <= Gamma_Q){
            double dcL = (cL_final - cL_init)/double(Gamma_Q);
            cL_quench += dcL;
            L_current = change_lattice(L_temp, L_current, cL_quench, N, Dx, Dy, lx, ly, dt);
            double V_current = vortices(L_current, N);
            write_quenchV << V_current << " ";
            //check cL values
            cout << cL_quench << "\n";}
        //after quench, lattice has noise corresponding to cL_final value
        else{
        L_current = change_lattice(L_temp, L_current, cL_final, N, Dx, Dy, lx, ly, dt);
        double V_current = vortices(L_current, N);
        write_quenchV << V_current << " ";}
    }
    write_quenchV.close();
}

int main(){
    //initial configuration of NxN lattice, each site with initial phase between 0 and 2pi chosen using random number generator

    //define parameters
    double dt = 0.05;
    double toleranceE = 1e-4;
    double cL[29] = {0, 0.5, 1, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.05, 3.15, 3.25, 3.35, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.85, 3.95, 4, 4.5, 5, 5.5, 6, 6.5, 7};
    double Et[29] = {800,800,800,800,800,800,800,800,900,950,100,1100,1200,1300,1400,1500,1600,1600,1600,1100,1100,1000,900,500,400,400,300,200,160};
    double ER[29] = {10,10,15,15,20,30,30,30,40,40,40,50,50,60,60,70,70,70,70,70,80,90,100,100,100,150,170,200,250};
    double Vt[29] = {700,700,700,700,700,800,800,800,800,800,800,900,900,900,1100,1200,1300,1500,1500,1100,1100,1000,900,500,400,400,300,200,160};
    double VR[29] = {25,30,20,40,20,20,30,40,30,40,40,50,40,50,50,40,50,50,50,60,50,60,60,100,120,140,180,200,250};

    for (int i = 0; i < 29; i++)
    {writeEV(Et[i], ER[i], cL[i], N, Dx, Dy, 0, 0, dt, true, false,false, toleranceE);
    writeEV(Vt[i], VR[i], cL[i], N, Dx, Dy, 0, 0, dt, false, true,false, toleranceE);}


/*
    //define parameters
    double dt = 0.05;
    double toleranceE = 1e-4;
    //writeEV(1400, 70, 3.75, N, Dx, Dy, 0.25, -0.25, dt, true, false, false, toleranceE);
    //writeEV(1100, 70, 3.85, N, Dx, Dy, 0.25, -0.25, dt, true, false, false, toleranceE);
    //writeEV(1400, 70, 3.5, N, Dx, Dy, 0.5, -0.5, dt, true, false, false, toleranceE);

    double lx = 0;
    double ly = 0;
    double cL_final = 2;
    double cL_init = 7;
    int init_iter = 200; //iterations to reach steady state for initial cL - known from previous data
    int final_iter = 700; //iterations to reach steady state for final cL - known from previous data

    //int Gamma_Qfinal = floor(log2(final_iter/2.0));
    for (int i = 1; i <= 8; i++){ // max power always 8 if considering final_iter < 1000
        int Gamma_Q = pow(2,i);
        writequench(cL_init, cL_final, init_iter, final_iter, Gamma_Q, N, Dx, Dy, lx, ly, dt);}
        */

}
