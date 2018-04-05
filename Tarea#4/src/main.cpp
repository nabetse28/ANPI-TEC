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
#include <QR.hpp>

#include "LUDoolittle.hpp"
#include "LUCrout.hpp"
#include "Matrix.hpp"
#include "QR.cpp"
using namespace std;




int main() {

    // Some example code
    anpi::Matrix<float> A = { {-1,-2,1,2},
                              { 2, 0,1,2},
                              {-1,-1,0,1},
                              { 1, 1,1,1} };

    //anpi::Matrix<float> A = { { 2, 0,1,2},{-1,-2,1,2},{ 1, 1,1,1},{-1,-1,0,1} };

    /*anpi::Matrix<float> H = {{2,0,1,2},{-0.5,-2,1.5,3},{-0.5,0.5,1.25,1.25},{0.5,-0.5,-0.2,0.8}};
    anpi::Matrix<float> E = {{4 ,-2, 1},
                           {20,-7,12},
                           {-8,13,17}};

    anpi::Matrix<int> F = {{3 , 2, 1},
                           {4 , 1, 2},
                           {1 , 3, 0}};

    anpi::Matrix<float> G = {{1, 0, 1}, {0, 2, 1}, {1, 1, 0}};

     */




    std::vector<size_t> bd;
    anpi::Matrix<float> LUD;
    anpi::Matrix<float> LD;
    anpi::Matrix<float> UD;
    anpi::Matrix<float> AD;



    std::vector<size_t> bc;
    anpi::Matrix<float> LUC;
    anpi::Matrix<float> LC;
    anpi::Matrix<float> UC;
    anpi::Matrix<float> AC;



    cout<<"-----------------------------LUDOOLITTLE-----------------------------------"<<endl;

    cout<<"Esta es la Matriz A que se le ingresa al metodo DooLittle: "<<endl;
    anpi::printMatrixDoolittle(A);
    anpi::luDoolittle<float>(A,LUD,bd);
    cout<<"Esta es la Matriz LU resultante en metodo DooLittle con su respectivo vector: "<<endl;
    anpi::printMatrixDoolittle(LUD);
    anpi::printVectorDoolittle(bd);
    cout<<endl;
    cout<<"Esta es la Matriz LU que se le ingresa al metodo UnpackDooLittle: "<<endl;
    anpi::printMatrixDoolittle(LUD);
    anpi::unpackDoolittle<float>(LUD,LD,UD);
    cout<<"Esta L obtenida despues del UnpackDooLittle: "<<endl;
    anpi::printMatrixDoolittle(LD);
    cout<<"Esta U obtenida despues del UnpackDooLittle: "<<endl;
    anpi::printMatrixDoolittle(UD);
    cout<<"Esta es la multiplicacion de L y U obtenida despues del UnpackDooLittle: "<<endl;
    AD = LD * UD;
    anpi::printMatrixDoolittle(AD);
    cout<<endl;

    cout<<"-----------------------------LUCROUT-----------------------------------"<<endl;

    cout<<"Esta es la Matriz A que se le ingresa al metodo Crout: "<<endl;
    anpi::printMatrixDoolittle(A);
    anpi::luCrout<float>(A,LUC,bc);
    cout<<"Esta es la Matriz LU resultante en metodo Crout con su respectivo vector: "<<endl;
    anpi::printMatrixDoolittle(LUC);
    anpi::printVectorDoolittle(bc);
    cout<<endl;
    cout<<"Esta es la Matriz LU que se le ingresa al metodo UnpackCrout: "<<endl;
    anpi::printMatrixDoolittle(LUC);
    anpi::unpackCrout<float>(LUC,LC,UC);
    cout<<"Esta L obtenida despues del UnpackCrout: "<<endl;
    anpi::printMatrixDoolittle(LC);
    cout<<"Esta U obtenida despues del UnpackCrout: "<<endl;
    anpi::printMatrixDoolittle(UC);
    cout<<"Esta es la multiplicacion de L y U obtenida despues del UnpackCrout: "<<endl;
    AC = LC * UC;
    anpi::printMatrixDoolittle(AC);
    cout<<endl;






    /// OBJETOS PARA LAS PRUEBAS
    QR<double> qrDouble;
    QR<float> qrFloat;


    /// MATRICES PARA PRUEBAS
    anpi::Matrix<double> ad = {{1, 2, 4}, {4, 2, 6}, {7, 6, 9}};
    anpi::Matrix<double> a(ad.rows(), ad.cols(), double(0));
    std::vector<double> vad = {1, 3, 0};
    std::vector<double> vadx = {0,0,0};



    anpi::Matrix<double> ad1 = {{3, 7, 1,6}, {12,5, 1, 9}, {8, 4,21, 6},{12,32,21,25}};
    anpi::Matrix<double> a1(ad1.rows(), ad1.cols(), double(0));
    std::vector<double> vad1 = {14, 23, 10,5};
    std::vector<double> vadx1 = {0,0,0,0};

    anpi::Matrix<double> ad2 = {{3, 7, 1,6,3}, {12,5, 1,8, 9}, {3,6,12,3,2},{6,12,32,7,8},{1,5,3,2,1}};
    anpi::Matrix<double> a2(ad1.rows(), ad1.cols(), double(0));
    std::vector<double> vad2 = {14, 23, 10,5,4};
    std::vector<double> vadx2 = {0,0,0,0,0};


    anpi::Matrix<float> bda = {{1, 2, 4}, {4, 2, 6}, {7, 6, 9}};
    anpi::Matrix<float> b(ad.rows(), ad.cols(), float(0));
    std::vector<float> vbd = {1, 3, 0};
    std::vector<float> vbdx = {0,0,0};


    anpi::Matrix<float> bd1 = {{3, 7, 1,6}, {12,5, 1, 9}, {8, 4,21, 6},{12,32,21,25}};
    anpi::Matrix<float> b1(ad1.rows(), ad1.cols(), float(0));
    std::vector<float> vbd1 = {14, 23, 10,5};
    std::vector<float> vbdx1 = {0,0,0,0};

    anpi::Matrix<float> bd2 = {{3, 7, 1,6,3}, {12,5, 1,8, 9}, {3,6,12,3,2},{6,12,32,7,8},{1,5,3,2,1}};
    anpi::Matrix<float> b2(ad1.rows(), ad1.cols(), float(0));
    std::vector<float> vbd2 = {14, 23, 10,5,4};
    std::vector<float> vbdx2 = {0,0,0,0,0};


    cout<<"-----------------------------QR-----------------------------------"<<endl;
    //-----------------------------------------QR
    std::cout << "---------------QR PRUEBAS CON TIPO DOUBLE" << std::endl;
    qrDouble.solveQR(ad,vadx,vad);
    std::cout << "SOLVE QR (Ax = b)" << std::endl;
    std::cout << "Matriz A" << std::endl;
    qrDouble.imprimirMatrix(ad);
    std::cout << "Vector b" << std::endl;
    qrDouble.imprimirVector(vad);
    std::cout << "Resultado vector x" << std::endl;
    qrDouble.imprimirVector(vadx);
    std::cout << std::endl;

    qrDouble.solveQR(ad1,vadx1,vad1);
    std::cout << "SOLVE QR (Ax = b)" << std::endl;
    std::cout << "Matriz A" << std::endl;
    qrDouble.imprimirMatrix(ad1);
    std::cout << "Vector b" << std::endl;
    qrDouble.imprimirVector(vad1);
    std::cout << "Resultado vector x" << std::endl;
    qrDouble.imprimirVector(vadx1);
    std::cout << std::endl;

    qrDouble.solveQR(ad2,vadx2,vad2);
    std::cout << "SOLVE QR (Ax = b)" << std::endl;
    std::cout << "Matriz A" << std::endl;
    qrDouble.imprimirMatrix(ad2);
    std::cout << "Vector b" << std::endl;
    qrDouble.imprimirVector(vad2);
    std::cout << "Resultado vector x" << std::endl;
    qrDouble.imprimirVector(vadx2);
    std::cout << std::endl;

    std::cout << "---------------QR PRUEBAS CON TIPO FLOAT" << std::endl;

    qrFloat.solveQR(bda,vbdx,vbd);
    std::cout << "SOLVE QR (Ax = b)" << std::endl;
    std::cout << "Matriz A" << std::endl;
    qrFloat.imprimirMatrix(bda);
    std::cout << "Vector b" << std::endl;
    qrFloat.imprimirVector(vbd);
    std::cout << "Resultado vector x" << std::endl;
    qrFloat.imprimirVector(vbdx);
    std::cout << std::endl;

    qrFloat.solveQR(bd1,vbdx1,vbd1);
    std::cout << "SOLVE QR (Ax = b)" << std::endl;
    std::cout << "Matriz A" << std::endl;
    qrFloat.imprimirMatrix(bd1);
    std::cout << "Vector b" << std::endl;
    qrFloat.imprimirVector(vbd1);
    std::cout << "Resultado vector x" << std::endl;
    qrFloat.imprimirVector(vbdx1);
    std::cout << std::endl;

    qrFloat.solveQR(bd2,vbdx2,vbd2);
    std::cout << "SOLVE QR (Ax = b)" << std::endl;
    std::cout << "Matriz A" << std::endl;
    qrFloat.imprimirMatrix(bd2);
    std::cout << "Vector b" << std::endl;
    qrFloat.imprimirVector(vbd2);
    std::cout << "Resultado vector x" << std::endl;
    qrFloat.imprimirVector(vbdx2);
    std::cout << std::endl;




}
  
