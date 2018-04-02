/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <iostream>

#include "LUDoolittle.hpp"
#include "Matrix.hpp"

using namespace std;


template <typename T>
void printMatrix(const anpi::Matrix<T>& A){

    int const Acols = A.cols();
    int const Arows = A.rows();

    for(int i=0; i<Arows;i++){
        string print=" ";
        for(int j=0; j<Acols;j++){
            if(j==0){
                //print+= "[ ";
                T t = A[i][j];
                print+= "[ " + to_string(t) + " ";
                //print+=" ";
            }else if(j==Acols-1){
                T t = A[i][j];
                print+= to_string(t) + " ]";
                //print+=" ]";

            }else{
                T t = A[i][j];
                print+= to_string(t) + " ";
                //print+= " ";

            }

        }
        cout<<print<<endl;
    }
}

int main() {

    // Some example code
    /*anpi::Matrix<double> A = { {-1,-2,1,2},
                              { 2, 0,1,2},
                              {-1,-1,0,1},
                              { 1, 1,1,1} };*/

    anpi::Matrix<float> A = { { 2, 0,1,2},{-1,-2,1,2},{ 1, 1,1,1},{-1,-1,0,1} };

    anpi::Matrix<double> E = {{4 ,-2, 1},
                           {20,-7,12},
                           {-8,13,17}};

    anpi::Matrix<int> F = {{3 , 2, 1},
                           {4 , 1, 2},
                           {1 , 3, 0}};

    //printMatrix<float>(A);

    std::vector<size_t> b;
    anpi::Matrix<float> LU;
    anpi::Matrix<float> LU1;
    anpi::Matrix<float> LU2;


    //printMatrix(A);


    /*anpi::Matrix<float> B = { {2,3,4},
                              {6,7,8},
                              {1,2,3}  };*/

    /*anpi::Matrix<float> C = { {1,2},
                              {2,1} };*/
    //int a = A.rows()*A.cols();

    /*vector<float> vector1 = {1,2,3,4};
    vector<float> vector2;







    anpi::Matrix<int> D;


    vector2 = A * vector1;


    anpi::Matrix<int> b = {{1,2},{1,2},{1,2}};
    anpi::Matrix<int> r1 = {{6,12},{6,12},{6,12}};


    D = b * r1;*/
    //D = B * C;

    /*anpi::luDoolittle<float>(A,LU,b);
    anpi::luDoolittle<double>(E,LU1,b);*/
    //anpi::luDoolittle<int>(F,LU2,b);

    double eps =  std::numeric_limits<float>::epsilon();

    anpi::luDoolittle<float>(A,LU,b);
    anpi::unpackDoolittle<float>(LU,LU1,LU2);
    anpi::Matrix<float> Ar = LU1 * LU2;

    for (size_t i=0;i<Ar.rows();++i) {
        for (size_t j=0;j<Ar.cols();++j) {
            if(abs(Ar(i,j)-A(i,j)) < eps){
                cout<<true<<endl;
            }
        }
    }
    


    //A[1][1] = 3;









    //anpi::Matrix<float> LU;

    //std::vector<size_t> p;
    //anpi::luDoolittle(A,LU,p);

    //return EXIT_FAILURE;
}
  
