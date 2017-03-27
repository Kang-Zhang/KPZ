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
const int N = 512; //lattice size
double Dx = 1; //diffusion constant
double Dy = 1; //diffusion constant
random_device rd;
mt19937 gen(rd());

//Set periodic boundary conditions:
int periodic(int i, int N){return (i % N + N) % N;}

typedef array <array<double, N> ,N> matrix; //shortcut for calling a data type

/* Function to simulate KPZ equation */
matrix change_lattice(matrix& L_new, matrix& L, double c_L, int N, double Dx, double Dy, double lx, double ly, double dt){ //reference to avoid wasting memory
    uniform_real_distribution<> dis(-0.5,0.5);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            //defining the terms in equation
            //top = [periodic(i-1,N)][j]; right = [i][periodic(j+1,N)]; bottom = [periodic(i+1,N)][j]; left = [i][periodic(j-1,N)]
            double diffX = Dx*(sin(L[i][j] - L[i][periodic(j+1,N)]) + sin(L[i][j] - L[i][periodic(j-1,N)]));
            double diffY = Dy*(sin(L[i][j] - L[periodic(i-1,N)][j]) + sin(L[i][j] - L[periodic(i+1,N)][j]));
            double nonlinX = (lx*0.5)*( cos(L[i][j] - L[i][periodic(j+1,N)]) + cos(L[i][j] - L[i][periodic(j-1,N)]) - 1);
            double nonlinY = (ly*0.5)*(cos(L[i][j] - L[periodic(i-1,N)][j]) + cos(L[i][j] - L[periodic(i+1,N)][j]) - 1);
            double noise = 2*M_PI*c_L*dis(gen);
            //Euler update for KPZ equation
            double L_step = dt*(diffX + diffY + nonlinX + nonlinY + noise);
            L_new[i][j] = fmod((L[i][j] - L_step),2.0*M_PI);
            if (L_new[i][j] < 0){L_new[i][j] += 2.0*M_PI;}}}
    return L_new;}

/*Function to calculate energy*/
double energy (matrix& L_new, int N, double Dx, double Dy){
    double E = 0; //keep track of energy
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
             E += -Dx*(cos(L_new[i][j] - L_new[i][periodic(j+1,N)]) + cos(L_new[i][j] - L_new[i][periodic(j-1,N)]))
                - Dy*(cos(L_new[i][j] - L_new[periodic(i-1,N)][j]) + cos(L_new[i][j] - L_new[periodic(i+1,N)][j]));}}
    return E;}

/*Function to calculate number of vortices*/
double vortices(matrix& L_new, int N){
    double num = 0; //keep track of number of vortices
    double d[4];
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            //top = [periodic(i-1,N)][j]; top_left = [periodic(i-1,N)][periodic(j-1,N)]; left = [i][periodic(j-1,N)]
            d[0] = L_new[periodic(i-1,N)][j] - L_new[i][j];
            d[1] = L_new[periodic(i-1,N)][periodic(j-1,N)] - L_new[periodic(i-1,N)][j];
            d[2] = L_new[i][periodic(j-1,N)] - L_new[periodic(i-1,N)][periodic(j-1,N)];
            d[3] = L_new[i][j] - L_new[i][periodic(j-1,N)];
            //ensure difference is in [-pi, pi]
            for(int k = 0; k < 4; k++){
                if(d[k] > M_PI){d[k] -= 2*M_PI;}
                else if(d[k] < -M_PI){d[k] += 2*M_PI;}}
            num += fabs(d[0]+d[1]+d[2]+d[3]);}} //keep track of total number of vortices
    return num/(2.0*M_PI);}

/*Function to calculate energy and number of vortices for different N and nonlinear constants lx and ly with a convergence test*/
void writeEV(int max_iters, int R, double c_L, int N, double Dx, double Dy, double lx, double ly, double dt, bool writeE, bool writeV, bool converge_test, double toleranceE)
{   //Write data and/or convergence test data to file:
    stringstream cL;    cL << fixed << setprecision(2) << c_L;
    stringstream l_x;   l_x << fixed << setprecision(2) << lx;
    stringstream l_y;   l_y << fixed << setprecision(2) << ly;
    stringstream Lsize; Lsize << N;
    ofstream write_outputE; ofstream write_outputV; ofstream write_outputE_t;   ofstream write_outputV_t;
    if (writeE){
        if(converge_test){write_outputE_t.open("KPZProjectData/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "E_t.txt");}
        else{write_outputE.open("KPZProjectData/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "E.txt");}}
    if (writeV){
        if(converge_test){write_outputV_t.open("KPZProjectData/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "V_t.txt");}
        else{write_outputV.open("KPZProjectData/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "V.txt");}}
    //Repeat for number of realisations R:
    for (int r = 0; r < R;r++){
        //Initialise random matrix for starting lattice configuration:
        uniform_real_distribution<> dis(0,2*M_PI);
        matrix L;
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                    L[i][j] = dis(gen);}}
        matrix L_new; //initialise random matrix for storing new lattice

        //Start of time evolution of lattice:
        int t = 0;
        double oldEdens = 0; //for convergence test of energy
        double oldV = 0; //for convergence test of vortices
        while (t < max_iters){
            t += 1; //increase time step (iterations) by 1
            L = change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt);
            if (writeE){double E = energy(L, N, Dx, Dy);
                if (converge_test){
                    if (fabs((E/(N*N)) - oldEdens) < toleranceE){write_outputE_t << E << " "<< t << "\n";
                        break;}
                    else{oldEdens = E/(N*N);}}
                else{write_outputE << E/(N*N) << " ";}}

            if (writeV){double V = vortices(L, N);
                if (converge_test){
                    if (fabs(V - oldV) == 0){
                        int t_i = t;
                        double Vtot = 0;
                        for (t = t_i; t < (t_i + 30);t++){ //need to compare next 30 V values to see if V is still the same
                            Vtot += V;
                            oldV = V;}
                        if (Vtot/30 == V){write_outputV_t << V << " "<< t << "\n";
                            break;}}
                    else{oldV = V;}}
                else{write_outputV << V << " ";}}
        } //end of one realisation, add new line to create a data matrix in text file
        if (writeE && converge_test == false){write_outputE << "\n";}
        if (writeV && converge_test == false){write_outputV << "\n";}
        } //end of all realisations, need to close all files
    if (writeE){
        if (converge_test){write_outputE_t.close();}
        else{write_outputE.close();}}
    if (writeV){
        if (converge_test){write_outputV_t.close();}
        else{write_outputV.close();}}}

/*Function to simulate finite rate quenching and outputs data to text file*/
void writequench(double cL_init, double cL_final, int init_iter, int Tau_Q, int N, double Dx,double Dy, double lx, double ly, double dt, int R){
    //Write data to file:
    stringstream cL_i;  cL_i << fixed << setprecision(2) << cL_init;
    stringstream cL_f;  cL_f << fixed << setprecision(2) << cL_final;
    stringstream l_x;   l_x << fixed << setprecision(2) << lx;
    stringstream l_y;   l_y << fixed << setprecision(2) << ly;
    stringstream Lsize; Lsize << N;
    stringstream TauQ;  TauQ << Tau_Q;
    stringstream realis;    realis << R;
    ofstream write_quenchV;
    write_quenchV.open("KPZProjectData/" + Lsize.str() + "-" + cL_i.str() + "-" + cL_f.str() + "-" + TauQ.str() + "-" + l_x.str() + "-" + l_y.str() + "-" + realis.str() + "V.txt");

    for (int r = 0; r < R; r++) {
        //Obtain initial matrix and initial number of vortices for initial cL value:
        matrix L; //to initialise random matrix for configuring initial lattice
        uniform_real_distribution<> dis(0,2*M_PI);
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                L[i][j] = dis(gen);}}
        matrix L_new; //initialise random matrix for storing new lattice
        //Obtain initial disordered config. for initial cL value
        for(int i = 0; i < init_iter; i++){
                L = change_lattice(L_new, L, cL_init, N, Dx, Dy, lx, ly, dt);} //steady state reached
        double V_init = vortices(L, N); //initial number of vortices at t = 0
        write_quenchV << V_init << " "; //record value at t = 0

        //Finite quenching:
        double cL_quench = cL_init;
        for(int t = 1; t <= Tau_Q; t++) {
            double dcL = (cL_final - cL_init)/double(Tau_Q);
            cL_quench += dcL; //decrease cL according to quench rate Tau_Q
            L = change_lattice(L_new, L, cL_init, N, Dx, Dy, lx, ly, dt);
            double V_current = vortices(L, N);
            write_quenchV << V_current << " "; // write value ^= current cL
        } //end of realisation & quench, lattice has noise ^= final cL
        write_quenchV << "\n";
        } //end of realisations, need to close files
        write_quenchV.close();}

int main(){
    //define parameters
    double dt = 0.05;
    double toleranceE = 1e-4;
    double cL_final = 2;
    double cL_init = 7;
    int init_iter = 200; //iterations to reach steady state for initial cL; known from previous data

    /*Write steady state data*/
    double cL[29] = {0, 0.5, 1, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.05, 3.15, 3.25, 3.35, 3.45, 3.5, 3.55,
    3.6, 3.65, 3.7, 3.75, 3.85, 3.95, 4, 4.5, 5, 5.5, 6, 6.5, 7}; //noise parameter
    double Et[29] = {600,600,700,700,700,800,800,800,900,950,1000,1100,1200,1300,1400,1500,
    1600,1600,1600,1100,1100,1000,900,500,400,400,300,200,160}; //iterations for energy
    double ER[29] = {10,10,10,15,15,20,20,30,30,30,40,40,40,50,50,50,60,60,60,60,50,70,80,
    80,90,100,150,170,200}; //realisations for energy
    double Vt[29] = {700,700,700,700,700,800,800,800,800,800,800,900,900,900,1100,1200,1300,
    1500,1500,1100,1100,1000,900,500,400,400,300,200,160}; //iterations for vortices
    double VR[29] = {20,30,30,30,30,40,40,40,50,50,50,50,60,60,60,60,70,70,70,80,80,
    80,80,90,110,140,160,180,200}; //realisations for vortices

    //Write energy density and vortices data for all cL values, l_x = 0 and l_y = 0
    for (int i = 0; i < 29; i++){
        writeEV(Et[i], ER[i], cL[i], N, Dx, Dy, 0, 0, dt, true, false,false, toleranceE);
        writeEV(Vt[i], VR[i], cL[i], N, Dx, Dy, 0, 0, dt, false, true,false, toleranceE);}

    /*Write finite rate quenching data*/
    for (int i = 9; i <= 13; i++){
        int Tau_Q = pow(2,i);
        writequench(cL_init, cL_final, init_iter, Tau_Q, N, Dx, Dy, 0, 0, dt, 60);}
}
