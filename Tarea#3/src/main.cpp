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
#include <math.h>
#include "../include/RootInterpolation.hpp"
#include "../include/RootBisection.hpp"
#include "../include/RootSecant.hpp"
#include "../include/RootBrent.hpp"
#include "../include/RootNewtonRaphson.hpp"
#include <../test/testRootFinders.cpp>
#include "../include/PlotPy.hpp"

using namespace std;

template<typename T>
T abs_x(T x){
    return abs(x)-exp(-x);
}


template<typename T>
T test(T x){

    static const T pi = 3.14159265358979323846264338327950288;
    return 0.5*exp(-x)-5.0*cos(pi*x);
}

int main() {
    //float a = anpi::rootNewtonRaphson<float>( abs_x<float> /*abs_x<float>*/, 1.0);


    // Put your main code in here


  //double a = anpi::rootInterpolation<double>( abs_x<double> /*abs_x<float>*/, 0.0 , 1.0 );

   // double a = anpi::rootBrent<double>(abs_x<double>,-3.1,-1.1);

  //double a = anpi::rootSecant<double>(abs_x<double> , 0.0, 1.0);

  //double a = anpi::rootBisection<double>(abs_x<double> , 0.0, 1.0);
  //std::cout <<  a << std::endl; // REMOVE-ME!*/




  return 0;
  //return EXIT_FAILURE;
}
  
