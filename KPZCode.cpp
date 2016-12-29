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
    //int periodic_index = (i % N + N) % N;
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
            if (L_new[i][j] < 0)
            {
                L_new[i][j] += 2.0*M_PI;
            }
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

            for(int k=0; k < 4; k++)
            {
                if(d[k] > M_PI)
                {
                    d[k] -= 2*M_PI;
                }
                else if(d[k] < -M_PI)
                {
                    d[k] += 2*M_PI;
                }

            }

            num += fabs(d[0]+d[1]+d[2]+d[3]);
        }}
    return num/2.0*M_PI;
}

void outputV(int max_iters, int R, double c_L, int N, double Dx,
                      double Dy, double lx, double ly, double dt)
{
    stringstream cL;
    cL << fixed << setprecision(2) << c_L;
    stringstream l_x;
    l_x << fixed << setprecision(2) << lx;
    stringstream l_y;
    l_y << fixed << setprecision(2) << ly;
    stringstream Lsize;
    Lsize << N;
    ofstream write_outputV ("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "V.txt");

    for (int r = 0; r<R;r++)
    {
        uniform_real_distribution<> dis(0,2*M_PI);
        matrix L;
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                    L[i][j] = dis(gen);
            }}
        matrix L_new;

        cout << r << "\n";
            for (int i = 0; i < max_iters; i++)
            {
                L = change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt);

                double V = vortices(L, N);
                write_outputV << V << " ";

            }
            write_outputV << "\n";
            }

    write_outputV.close();

}

void outputE(int max_iters, int R, double c_L, int N, double Dx,
                      double Dy, double lx, double ly, double dt)
{
    stringstream cL;
    cL << fixed << setprecision(2) << c_L;
    stringstream l_x;
    l_x << fixed << setprecision(2) << lx;
    stringstream l_y;
    l_y << fixed << setprecision(2) << ly;
    stringstream Lsize;
    Lsize << N;
    ofstream write_outputE ("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "E.txt");

    for (int r = 0; r<R;r++)
    {
        uniform_real_distribution<> dis(0,2*M_PI);
        matrix L;
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                    L[i][j] = dis(gen);
            }}
        matrix L_new;

        cout << r << "\n";
            for (int i = 0; i < max_iters; i++)
            {
                L = change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt);

                double E = energy(L, N, Dx, Dy);
                write_outputE << E/(N*N) << " ";
            }
            write_outputE << "\n"; //new line after each realisation
            }

    write_outputE.close();

}

void outputEV(int max_iters, int R, double c_L, int N, double Dx,
                      double Dy, double lx, double ly, double dt)
{
    stringstream cL;
    cL << fixed << setprecision(2) << c_L;
    //cL << fixed << setprecision(1) << c_L+0.01;
    stringstream l_x;
    l_x << fixed << setprecision(2) << lx;
    stringstream l_y;
    l_y << fixed << setprecision(2) << ly;
    stringstream Lsize;
    Lsize << N;
    ofstream write_outputE ("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "E.txt");
    ofstream write_outputV ("KPZProjectGraphs/" + Lsize.str() + "-" + cL.str() + "-" + l_x.str() + "-" + l_y.str() + "V.txt");

    //array<double, max_iters> meanEdensity;
    //array<double, max_iters> meanVdensity; //this seems to give nan in E file

    for (int r = 0; r<R;r++)
    {
        uniform_real_distribution<> dis(0,2*M_PI);
        matrix L;
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                    L[i][j] = dis(gen);
            }}
        matrix L_new;

        cout << r << "\n";
            for (int i = 0; i < max_iters; i++)
            {
                L = change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt);

                double E = energy(L, N, Dx, Dy);
                write_outputE << E/(N*N) << " ";

                //cout << E/(N*N) << "\n";
                // meanEdensity[i] += E/(N*N);
                //cout << meanEdensity[i] << "\n";

                double V = vortices(L, N);
                write_outputV << V << " ";
                //meanVdensity[i] += V/(N*N);
            }
            write_outputE << "\n"; //new line after each realisation
            write_outputV << "\n";
            }

    //for (int i = 0; i < max_iters; i++ ){
        //meanEdensity[i] /= R;
        //meanVdensity[i] /= R;
        //write_outputE << meanEdensity[i] << " ";
        //cout << meanEdensity[i] << "\n";
        //write_outputV << meanVdensity[i] << " ";}

    write_outputE.close();
    write_outputV.close();


}

int main(){
    //initial configuration of NxN lattice, each site with initial phase between 0 and 2pi chosen using random number generator

    //define parameters
    double c_L = 3.35;
    double dt = 0.05;
    const int max_iters = 100;
    int R = 10; //number of realisations

    //for (double c_L=0; c_L<=7; c_L+=0.5){
    outputEV(max_iters, R, c_L, N, Dx, Dy, lx, ly, dt);
    //outputE(max_iters, R, c_L, N, Dx, Dy, lx, ly, dt);
    //outputV(max_iters, R, c_L, N, Dx, Dy, lx, ly, dt);

    //}


}
