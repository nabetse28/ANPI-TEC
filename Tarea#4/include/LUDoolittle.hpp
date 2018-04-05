/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi {


    /**
     * It fills the vector with 0...n where n is the number of
     * columns in the matrix.
     *
     * @param[in] Vector
     * @param[in] n (columns number)
     *
     */
    void fillVector(std::vector<size_t>& v, size_t n){
        for(size_t i=0;i<n; i++){
            v[i] = i;
        }
    }

    /**
     * It prints the vector like a vector in console.
     *
     * @param[in] Vector
     *
     */
    void printVectorDoolittle(std::vector<size_t>& v){
        string prsize_t = "Vector: ";
        size_t const v_size = v.size();
        for(size_t i=0;i<v_size;i++){
            if(i == 0){
                prsize_t += "[ " + to_string(v[i]) + " ";
            }else if(i==v_size-1){
                prsize_t += to_string(v[i]) +" ]";
            }else{
                prsize_t += to_string(v[i]) +" ";
            }


        }
        cout<<prsize_t<<endl;

    }

    /**
     * It checks if the number a is in the vector
     *
     * @param[in] Vector
     * @param[in] a (number)
     * @param[out] return a bool (false or true)
     *
     */
    template <typename T>
    bool checkPermut(T a,std::vector<T>& v){
        bool res = false;
        for(size_t i= 0; i<v.size();i++){
            if(v[i]==a){
                res = true;
            }
        }
        return res;
    }

    /**
     * Prints the Matrix A like a matrix in console.
     *
     * @param[in] A a square matrix
     */
    template <typename T>
    void printMatrixDoolittle(const anpi::Matrix<T>& A){

        size_t const Acols = A.cols();
        size_t const Arows = A.rows();
        for(size_t i=0; i<Arows;i++){
            string prsize_t=" ";
            for(size_t j=0; j<Acols;j++){
                if(j==0){

                    T t = A[i][j];
                    prsize_t+= "[ " + to_string(t) + " ";

                }else if(j==Acols-1){
                    T t = A[i][j];
                    prsize_t+= to_string(t) + " ]";

                }else{
                    T t = A[i][j];
                    prsize_t+= to_string(t) + " ";

                }

            }
            cout<<prsize_t<<endl;
        }
        cout<<endl;
    }

    /**
     * Auxiliary method used to debug LU decomposition.
     *
     * It separates a packed LU matrix size_to the lower triangular matrix
     * L and the upper triangular matrix U, such that the diagonal of L
     * is composed by 1's.
     */
    template<typename T>
    void unpackDoolittle(const Matrix<T>& LU,
                         Matrix<T>& L,
                         Matrix<T>& U) {


        size_t const Acols = LU.cols();
        size_t const Arows = LU.rows();

        anpi::Matrix<T> newA = LU;


/*
        cout<<"Esta es la matriz LU: "<<endl;
        printMatrixDoolittle(newA);
*/

        // Se compara si es una matriz cuadrada.
        if(Acols!=Arows){
            throw anpi::Exception("No es una matriz rectangular");
        }

        // Se inicializan las matrices L y U con la cantidad de
        // filas y columnas de LU.
        L.allocate(newA.rows(),newA.cols());
        U.allocate(newA.rows(),newA.cols());

        // Se realiza el desempaquetamiento de la matriz LU.
        for(size_t k=0;k<Acols;k++){
            for(size_t r=0;r<Arows;r++){
                if(k==r){
                    L[k][r] = T(1);
                    U[k][r] = newA[k][r];
                }else if(k<r){
                    T factor = newA[r][k]/*((newA[r][k])/(newA[k][k]))*/;
                    L[r][k] = factor;
                    U[k][r] = newA[k][r];
                    /*for(size_t c=0;c<Arows;c++){
                        newA[r][c] = ((newA[r][c])-((factor)*(newA[k][c])));
                        U = newA;
                    }*/
                }
            }
        }

/*
        cout<<"Esta es la matriz L: "<<endl;
        printMatrixDoolittle(L);
        cout<<"Esta es la matriz U: "<<endl;
        printMatrixDoolittle(U);
        anpi::Matrix<T> LU1 = L * U;
        cout<<"Esta es la matriz LU multiplicada"<<endl;
        printMatrixDoolittle(LU1);
*/

    }

    /**
     * Decompose the matrix A size_to a lower triangular matrix L and an
     * upper triangular matrix U.  The matrices L and U are packed size_to
     * a single matrix LU.
     *
     * The L matrix will have in the Doolittle's LU decomposition a
     * diagonal of 1's
     *
     * @param[in] A a square matrix
     * @param[out] LU matrix encoding the L and U matrices
     * @param[out] permut permutation vector, holding the indices of the
     *             original matrix falling size_to the corresponding element.
     *             For example if permut[5]==3 holds, then the fifth row
     *             of the LU decomposition in fact is dealing with the third
     *             row of the original matrix.
     *
     * @throws anpi::Exception if matrix cannot be decomposed, or input
     *         matrix is not square.
     */
    template<typename T>
    void luDoolittle(const Matrix<T>& A,
                     Matrix<T>& LU,
                     std::vector<size_t>& permut) {

        // Se asignan variables para trabajar la funcion.
        anpi::Matrix<T> L;
        anpi::Matrix<T> U;
        anpi::Matrix<T> newA;
        vector<size_t> v1(A.cols());


        newA = A;

        // Se inicializan las matrices con la cantidad correspondiente de
        // filas y columnas de A.
        L.allocate(A.rows(),A.cols());
        U.allocate(A.rows(),A.cols());
        LU.allocate(A.rows(),A.cols());
        newA.allocate(A.rows(),A.cols());

        size_t const Acols = A.cols();
        size_t const Arows = A.rows();

        // Se compara si es una matriz cuadrada.
        if((newA.rows()!=newA.cols())){
            throw anpi::Exception("No es una matriz rectangular");
        }

        // Se rellena el vector con la funcion fill.
        fillVector(v1,Acols);

/*
        cout<<"Esta es la matriz A ingresada: "<<endl;
        printMatrixDoolittle(newA);

*/

        // Se ordena la matriz conforme el pivoteo, donde se pone el mayor de cada
        // columna de primero hasta ordenarla.
        for(size_t i= 0; i<Acols;i++){
            size_t pos = T(0);
            T colm = newA[i][i];
            for(size_t j = i; j<Arows;j++){
                if(abs(newA[j][i])>=abs(colm)) {
                    pos = j;
                    colm = newA[j][i];
                    size_t temp = v1[i];
                    v1[i] = v1[j];
                    v1[j] = temp;
                }
            }
            for(size_t m = 0; m<Arows;m++){
                for(size_t n = 0; n<Acols;n++){
                    if(m==pos){
                        T temp = newA[i][n];
                        newA[i][n] = newA[m][n];
                        newA[m][n] = temp;
                    }
                }
            }

            //printMatrixDoolittle(newA);
        }


        // Se le asignan los valores obtenidos al vector permut.
        permut = v1;


        // Se realiza la funcion de Doolittle
        for(size_t k=0;k<Acols;k++){
            for(size_t r=0;r<Arows;r++){
                if(k==r){
                    L[k][r] = T(1);
                }else if(k<r){
                    T factor = ((newA[r][k])/(newA[k][k])); // Se saca el factor y se agrega a la Matriz L.
                    L[r][k] = factor;
                    for(size_t c=0;c<Arows;c++){
                        newA[r][c] = ((newA[r][c])-((factor)*(newA[k][c]))); // Se obtiene el resultado de la multiplicacion del factor * los datos de la fila correspondiente.

                    }
                }
            }
        }



        // Se le asigna a LU su valor real.
        for(int i=0;i<Acols;i++){
            for(int j=0;j<Arows;j++){
                if(i==j){
                    LU[i][j] = newA[i][j];
                }else{
                    LU[i][j] = newA[i][j] + L[i][j];
                }

            }
        }
        //LU = L + U;

        // Se compara LU con A a ver si tienen la misma cantidad de filas y columnas.
        if((LU.cols()!=newA.cols())||(LU.rows()!=newA.rows())){
            throw anpi::Exception("La matriz LU no coincide con A");
        }


/*
        cout<<"Esta es la matriz L: "<<endl;
        printMatrixDoolittle(L);
        cout<<"Esta es la matriz U: "<<endl;
        printMatrixDoolittle(newA);
        cout<<"Esta es la matriz LU con su respectivo vector de permutacion: "<<endl;
        printMatrixDoolittle(LU);
        printVectorDoolittle(v1);
*/

    }

}

#endif

