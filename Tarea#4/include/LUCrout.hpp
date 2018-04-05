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
#include <iostream>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_CROUT_HPP
#define ANPI_LU_CROUT_HPP

namespace anpi {

    template <typename T>
    void printMatrixCrout(const anpi::Matrix<T>& A){

        size_t const Acols = A.cols();
        size_t const Arows = A.rows();

        for(size_t i=0; i<Arows;i++){
            std::string prsize_t=" ";
            for(size_t j=0; j<Acols;j++){
                if(j==0){

                    T t = A[i][j];
                    prsize_t+= "[ " + std::to_string(t) + " ";

                }else if(j==Acols-1){
                    T t = A[i][j];
                    prsize_t+= std::to_string(t) + " ]";

                }else{
                    T t = A[i][j];
                    prsize_t+= std::to_string(t) + " ";

                }

            }
            std::cout<<prsize_t<<std::endl;
        }
    }

    /**
     * Auxiliary method used to debug LU decomposition.
     *
     * It separates a packed LU matrix into the lower triangular matrix
     * L and the upper triangular mastrix U, such that the diagonal of U
     * is composed by 1's.
     */
    template<typename T>
    void unpackCrout(const Matrix<T>& LU,
                     Matrix<T>& L,
                     Matrix<T>& U) {

        size_t const N = LU.rows();
        anpi::Matrix<T> newLU = LU;

        if(LU.cols()!=LU.rows()){
            throw anpi::Exception("No es una matriz rectangular");
        }

        L.allocate(newLU.rows(), newLU.cols());
        U.allocate(newLU.rows(), newLU.cols());

        for (int i = 0 ; i < N ; i++){      //It iterates by rows
            int j = 0;
            while (j != i){                //Iterates until the diagonal.
                L(i,j) = newLU(i,j);       //The elements from below must be in the L matrix.
                j++;
            }
            /*L(i,j) = newLU(i,j);           //In the diagonal, the elements must be in the L matrix.
            U(i,j) = 1;                    //The U matrix has 1s in its diagonal.
             */
            U(i,j) = newLU(i,j);
            L(i,j) = 1;
            j++;
            while (j != N){               //Iterates from the diagonal until the end.
                U(i,j) = newLU(i,j);      //The elements in top of the diagonal must be in the U matrix.
                j++;
            }
        }
        /*
        std::cout << "Matriz L"<<std::endl;
        printMatrixDoolittle(L);
        std::cout << "Matriz U"<<std::endl;
        printMatrixDoolittle(U);
        std::cout<< "Matriz LU"<<endl;
        anpi::Matrix<T> LU1 = L*U;
        printMatrixCrout(LU1);
         */
    }

    /**
     * Decompose the matrix A into a lower triangular matrix L and an
     * upper triangular matrix U.  The matrices L and U are packed into
     * a single matrix LU.
     *
     * Crout's way of packing assumes a diagonal of
     * 1's in the U matrix.
     *
     * @param[in] A a square matrix
     * @param[out] LU matrix encoding the L and U matrices
     * @param[out] permut permutation vector, holding the indices of the
     *             original matrix falling into the corresponding element.
     *             For example if permut[5]==3 holds, then the fifth row
     *             of the LU decomposition in fact is dealing with the third
     *             row of the original matrix.
     *
     * @throws anpi::Exception if matrix cannot be decomposed, or input
     *         matrix is not square.
     */
    template<typename T>
    void luCrout(const Matrix<T>& A,
                 Matrix<T>& LU,
                 std::vector<size_t>& permut) {

        int const N = A.rows();                      //Definition of N. The matrix is a square matrix.
        anpi::Matrix<T> newA;
        newA = A;
        LU.allocate(A.rows(),A.cols());

        /*
        cout<<"La matriz ingresada es: "<<endl;
        printMatrixCrout(A);

         */

        int const Acols = A.cols();
        int const Arows = A.rows();

        if(A.cols()!=A.rows()){
            throw anpi::Exception("No es una matriz rectangular");
        }


        std::vector<size_t> v1(Acols);

        for (int i = 0 ; i <= N-1 ; i++){
            v1[i] = i;            //Construction of the permutation vector.
        }


        for(int i= 0; i<Acols;i++){
            int pos = T(0);
            T colm = newA[i][i];
            for(int j = i; j<Arows;j++){
                if(abs(newA[j][i])>=abs(colm)) {
                    pos = j;
                    colm = newA[j][i];
                    size_t temp = v1[i];
                    v1[i] = v1[j];
                    v1[j] = temp;
                }
            }
            for(int m = 0; m<Arows;m++){
                for(int n = 0; n<Acols;n++){
                    if(m==pos){
                        T temp = newA[i][n];
                        newA[i][n] = newA[m][n];
                        newA[m][n] = temp;
                    }
                }
            }

        }

        permut = v1;



        for (int j = 0; j < N; j++){       //Applies for both iterations.

            ////////////// U iteration //////////////

            for (int i = 0 ; i <= j ; i++){
                T result; //Result of the addition that has to be done.
                if (i == 0) //If the upper limit is 0, the addition is 0 too.
                    result = T(0);
                else{
                    for (int k = 0 ; k <= i ; k++){
                        result+= LU(i,k)*LU(k,j); //Addition function.
                    }
                }
                LU(i,j) = A(i,j)-result;
            }

            ///////////// L iteration ///////////////

            for (int i = j+1 ; i < N; i++){
                T result; //Result of the addition that has to be done.
                if (j == 0) //If the upper limit is 0, the addition is 0 too.
                    result = T(0);
                else{
                    for (int k = 0 ; k <= j ; k++){
                        result+= LU(i,k)*LU(k,j); //Addition function.
                    }
                }
                LU(i,j) = (A(i,j) - result)/LU(j,j);
            }
        }

        /*
        cout<<"LU: "<<endl;
        printMatrixCrout(LU);
        */
        if((LU.cols()!=A.cols())||(LU.rows()!=A.rows())){
            throw anpi::Exception("La matriz LU no coincide con A");
        }
    }
}
#endif
