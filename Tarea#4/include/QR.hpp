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

#ifndef ANPI_QR_HPP
#define ANPI_QR_HPP

namespace anpi {


    /**
     *  Funcion usada para calcular el determinante de una matriz
     * @param anpi::Matrix<T> A Matriz a calcular determinante
     * @return Retorna el determinante de A
     * */


    template <typename T>
    T cofactor(anpi::Matrix<T> A, T fila, T columna)
    {
        anpi::Matrix<T> subMatriz(A.rows() - 1, A.cols() - 1, anpi::DoNotInitialize);

        int x = 0;
        int y = 0;
        for (int i = 0; i < A.rows(); ++i)
        {
            for (int j = 0; j < A.cols(); ++j)
            {
                if (i != fila && j != columna)
                {
                    subMatriz[x][y] = A[i][j];
                    y = y + 1;
                    if (y >= subMatriz.cols())
                    {
                        x = x + 1;
                        y = 0;
                    }
                }
            }
        }

        T tem = getDeterminante(subMatriz);
        return (pow(-1, fila + columna)) * tem;
    }

    template <typename T>
    T getDeterminante(anpi::Matrix<T> A)
    {
        T resultado = 0;
        if (A.cols() == 1)
        {
            resultado = A[0][0];
        }
        else
        {
            for (int j = 0; j < A.cols(); ++j)
            {
                T tem = cofactor<T>(A, 0, j);
                resultado = resultado + (A[0][j] * tem);
            }
        }
        return resultado;
    }

    template <typename T>
    T getModule(anpi::Matrix<T> A,int fila){
        int x = A.cols();
        T resultado = 0;

        for(int j=0; j<x; ++j){
            resultado = resultado + (A[fila][j]*A[fila][j]);
        }
        resultado = sqrt(resultado);
        return resultado;
    }

    /**
     * Funcion encargada de obtener la transpuesta de una matriz
     * @param A Matriz a la cual se le determina la transpuesta
     * @return La matriz A transpuesta
     */

    template <typename T>
    anpi::Matrix<T> getTranspuesta(anpi::Matrix<T>& A){
        anpi::Matrix<T> resultado;
        resultado.allocate(A.cols(),A.rows());
        for(int i=0;i<A.rows();i++){
            for(int j=0;j<A.cols();j++){
                resultado[j][i]=A[i][j];
            }
        }
        return  resultado;

    }

    /**
     * Funcion que obtiene un fila especifica de una matriz
     * @param A Matriz a la cual se le obtiene la fila
     * @param numFila El numero de fila a obtener
     * @return La fila numero numFila en la matriz A
     */

    template <typename T>
    anpi::Matrix<T> getFila(anpi::Matrix<T>& A,int numFila){

        int size_Y = A.cols();
        anpi::Matrix<T> resultado;
        resultado.allocate(1,size_Y);
        for(int j=0;j<size_Y;++j){
            resultado[0][j]=A[numFila][j];
        }
        return resultado;
    }

    /**
     * Funcion que construye una matriz identidad de un tamaño especifico
     * @param size El tamaño de la matriz identidad deseada
     * @return La matriz identidad de tamaño size
     */

    template <typename T>
    anpi::Matrix<T> getIdentidad(int size){
        anpi::Matrix<T> resultado;
        resultado.allocate(size,size);
        for(int i=0;i<size;++i){
            for(int j=0;j<size;++j){
                if(i==j){resultado[i][j]=1;}
                else{resultado[i][j]=0;}
            }}
        return resultado;
    }

    /**
     * Funcion encargada de multiplicar una matriz con un determinado dato
     * @param dato Valor o valores a multiplacar con la matriz
     * @param A Matriz a multiplicar
     * @return el resultado de la multiplicacion de mat con data
     */
    template <typename T>
    anpi::Matrix<T> getMultiplicacion(T dato,anpi::Matrix<T> A){
        anpi::Matrix<T> resultado;
        resultado.allocate(A.rows(),A.cols());
        for(int i=0;i<A.rows(); ++i){
            for(int j=0; j<A.cols();++j){
                resultado[i][j]= A[i][j]*dato;
            }}
        return resultado;
    }

    /**
     * Funcion encargada de reducir la matriz una unidad
     * @param mat Matriz a reducir
     * @return La matriz reducida
     */

    template <typename T>
    anpi::Matrix<T> getReducida(anpi::Matrix<T> A){
        anpi::Matrix<T> resultado;
        resultado.allocate(A.rows()-1,A.cols()-1);
        for(int i=1;i<A.rows(); ++i){
            for(int j=1; j<A.cols();++j){
                resultado[i-1][j-1]= A[i][j];
            }}
        return resultado;
    }

    template <typename T>
    anpi::Matrix<T> getCompleta(anpi::Matrix<T> A,anpi::Matrix<T> iden){
        anpi::Matrix<T> result;
        result.allocate(iden.rows(),iden.cols());
        result = iden;
        int dif = iden.rows()-A.rows();

        for(int i=dif;i<result.rows(); ++i){
            for(int j=dif; j<result.cols();++j) {
                result[i][j] = A[i-dif][j-dif]; }}
        return result;
    }

    /**
     * Funcion principal de la descomposicion QR de una matrix
     * @param A Matriz
     * @param Q Matriz
     * @param R Matriz
     */
    template <typename T>
    void qr(const anpi::Matrix<T> &A, anpi::Matrix<T> &Q, anpi::Matrix<T> &R){
        if (A.rows() != A.cols())
        {
            throw std::invalid_argument("Para la descomposicion QR la matriz debe de ser cuadrada");
        }
        if (getDeterminante(A) == 0)
        {
            throw std::invalid_argument("Para la descomposicion QR la matriz NO debe de ser singular");
        }
        int size_x = A.rows();
        anpi::Matrix<T> Qn;
        anpi::Matrix<T> temp;
        anpi::Matrix<T> Vector;
        Qn = A;
        anpi::Matrix<T> roots;

        for(int i=0;i<size_x-1;++i){
            if(i>0){Qn = getReducida(Qn);}
            temp = Qn;
            Qn = getTranspuesta(Qn);
            anpi::Matrix<T> iden;
            iden = getIdentidad<T>(Qn.rows());
            Vector = getFila(Qn,0);
            T mod1 = getModule(Vector,0);
            anpi::Matrix<T> e =  getMultiplicacion(mod1,getFila(iden,0));
            anpi::Matrix<T> u = Vector-e;
            mod1 = getModule(u,0);
            T data;
            data= (T(2)/(mod1*mod1));
            anpi::Matrix<T> v = getTranspuesta(u)*u;
            v = getMultiplicacion(data,v);
            v = iden - v;
            u = getCompleta<T>(v, getIdentidad<T>(size_x));
            roots[i] = u;
            Qn = v*temp;
        }
        anpi::Matrix<T> res;
        res = roots[0];
        for(int j=1;j<size_x-1;++j){
            res = res*roots[j];
        }
        Q = res;

        R = getTranspuesta(res)*A;
    }
    template <typename T>
    anpi::Matrix<T> superMatriz(anpi::Matrix<T>& a,T n){
        T y = n*T(2);
        anpi::Matrix<T> b(n,y);

        for (int i = 0 ; i < n ; ++i) {
            for (int j = 0; j < n ; ++j) {
                b[i][j] = a[i][j];
            }
        }

        for (int i = 0 ; i < n ; ++i) {
            for (int j = n; j < n*2 ; ++j) {
                if(j==n+i){
                    b[i][j]=1;
                }
                else{ b[i][j] = 0;
                }}
        }


        return b;
    }


    /**
     * Funcion encargada de calcular Gauss Jordan necesario para la descomposicion QR
     * @param b Matriz
     * @param n
     * @return
     */
    template <typename T>
    anpi::Matrix<T> gaussJordan(anpi::Matrix<T>& b,int n){
        anpi::Matrix<T> A = superMatriz(b,n);
        int i, j, k;
        anpi::Matrix<T> c;
        c.allocate(n,n);
        T  factor(0);

        for (k = 0; k < n; k++) {

            factor = A[k][k];

            for (j = k; j < n*2; j++) {

                A[k][j] = A[k][j]/factor;

            }


            for (i = 0; i < n ; i++) {

                if (i != k) {

                    factor = A[i][k];

                    for (j = 0; j < 2*n; j++) {

                        A[i][j] = A[i][j] - factor * A[k][j];

                    }

                }

            }

        }
        for (int i = 0; i <n; ++i) {
            for (int j = 0; j <n ; ++j) {
                c[i][j] = A[i][j+n];
            }
        }
        return c;

    }

    /**
     * Funcion que convierte un vector en una matriz
     * @param x Vector a convertir
     * @param n Tamano del vector
     * @return El vector como matriz
     */

    template <typename T>
    anpi::Matrix<T> convertVector(std::vector<T> x, int n){
        anpi::Matrix<T> a;
        a.allocate(n,1);
        for (int i = 0; i < n; ++i) {
            a[i][0] = x[i];
        }
        return a;
    }

    /**
    * Funcion que convierte una matriz en un vector
    * @param x Matriz a convertir
    * @param n Tamano de la matriz
    * @return La matriz como vector
    */

    template <typename T>
    std::vector<T> convertMat(anpi::Matrix<T> a, int n){
        std::vector<T> x;
        for (int i = 0; i < n; ++i) {
            x.push_back(a[i][0]);
        }
        return x;
    }

    /**
     * Funcion que resuelve un sistema de ecuaciones lineales (Ax=b)
     * @param A Matriz
     * @param x Vector
     * @param b Vector
     * @return El resultado de Ax=b
     */
    template <typename T>
    bool solveQR(const anpi::Matrix<T>& A, std::vector<T> &x, std::vector<T> &b){

        if (A.rows() != A.cols())
        {
            throw std::invalid_argument("Para la descomposicon QR la matriz debe de ser cuadrada");
        }
        if (getDeterminante(A) == 0)
        {
            throw std::invalid_argument("Para ela descomposicon QR la matriz NO debe de ser singular");
        }
        anpi::Matrix<T> q = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        anpi::Matrix<T> r = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        anpi::Matrix<T> r_inv;
        anpi::Matrix<T> _x;
        _x = convertVector(x,x.size());
        anpi::Matrix<T> _b;
        _b = convertVector(b,b.size());
        anpi::Matrix<T> q_trans;
        qr(A,q,r);
        r_inv = gaussJordan(r,r.rows());
        q_trans = getTranspuesta(q);
        _x= (r_inv*q_trans);
        _x = _x*_b;
        x = convertMat(_x,_x.rows());

    }


    /**
     * Funcion encargada de imprimir una matriz
     * @param A Matriz a mostrar en consola
     */
    template <typename T>
    void imprimirMatrix(const anpi::Matrix<T>& A){
        for (int i = 0; i < A.rows(); ++i)
        {
            std::cout<<"[   ";
            for (int j = 0; j < A.cols(); ++j)
            {
                std::cout << A[i][j] << "   ";
            }
            std::cout<<"]"<<std::endl;
        }
    }


    /**
     * Funcion encargada de imprimir un vector
     * @param v Vector a imprimir en consola
     */
    template <typename T>
    void imprimirVector(std::vector<T>& vector){
        for (T v: vector)
        {
            std::cout<<"[ ";
            std::cout<<v<<" ]"<<std::endl;
        }

    }


    /**
	 *  Funcion usada para calcular la norma matricial Frobenius de una matriz
	 * @param anpi::Matrix<T> mat Matriz a calcular norma
	 * @return Retorna la norma matricial Frobenius
	 * */
    template <typename T>
    T normaFrob(anpi::Matrix<T> A)
    {

        int r = A.rows();
        int c = A.cols();

        T sum = 0;
        for (int i = 0; i < r; ++i)
        {
            for (int j = 0; j < c; ++j)
            {
                sum = sum + ((A[i][j]) * (A[i][j]));
            }
        }

        // std::cout << sum <<std::endl;
        sum = sqrt(sum);
        // std::cout << sum <<std::endl;
        return sum;
    }


}

#endif

