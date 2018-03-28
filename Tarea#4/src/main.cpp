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

int main() {

    // Some example code
    anpi::Matrix<float> A = { {-1,-2,1,2},
                              { 2, 0,1,2},
                              {-1,-1,0,1},
                              { 1, 1,1,1} };


    anpi::Matrix<float> B = { {2,3,4},
                              {6,7,8},
                              {1,2,3}  };

    anpi::Matrix<float> C = { {1,2},
                              {2,1},
            /*{1,2} */};
    //int a = A.rows()*A.cols();

    vector<float> vector1 = {1,2,3,4};
    vector<float> vector2;







    anpi::Matrix<int> D;


    vector2 = A * vector1;


    anpi::Matrix<int> b = {{1,2},{1,2},{1,2}};
    anpi::Matrix<int> r1 = {{6,12},{6,12},{6,12}};


    D = b * r1;
    //D = B * C;



    cout<<vector2[0]<<endl;


    //anpi::Matrix<float> LU;

    //std::vector<size_t> p;
    //anpi::luDoolittle(A,LU,p);

    //return EXIT_FAILURE;
}
  
