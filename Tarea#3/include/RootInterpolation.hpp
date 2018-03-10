/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>

#include "Exception.hpp"

#ifndef ANPI_ROOT_INTERPOLATION_HPP
#define ANPI_ROOT_INTERPOLATION_HPP

using namespace std;


namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], by means of the interpolation method.
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootInterpolation(const function<T(T)>& funct,T xl,T xu,const T eps = sqrt(numeric_limits<T>::epsilon())) {


        const T maxi = numeric_limits<T>::digits;

        T xr = xl; // Se define el xr la raiz que vamos a encontrar como x mas grande
        T fl = funct(xl); // Se evalua el x mas grande en la funcion
        T fu = funct(xu); // Se evalua el x mas pequeno en la funcion

        //T esp = sqrt(numeric_limits<T>::epsilon());
        //const int maxi = numeric_limits<T>::digits;

        T ea = T();

        int iu(0),il(0); // se le da valores a iu e il para ver si hay estancamiento

        if((funct(xl)*funct(xu)>T(0))){
            throw anpi::Exception("Unenclosed root"); // Si la raiz no se encuentra

        }

        if(xl >= xu){
            throw anpi::Exception("Inverted interval"); // Arroja el error si los intervalos estan invertidos
        }

        for(int i=maxi; i>0; --i){
            T xrold(xr); // Se define un x viejo para el calculo del error
            xr = xu-fu*(xl-xu)/(fl-fu);

            T fr = funct(xr);

            if(abs(xr)>numeric_limits<T>::epsilon()){ // Esto es para evitar una division entre 0
                ea = abs((xr-xrold)/xr)*T(100); // Se define el errpr aproximado
            }

            T cond = fl*fr; // Para saber cual subintervalo contine la raiz

            if(cond < T(0)){ // El lado menor tiene la raiz
                xu=xr;
                fu=fr;
                iu=0;
                il++;
                if(il>=2){
                    fl /= T(2);


                }

            }else if (cond > T(0)){ // El lado mayo contiene la raiz
                xl = xr;
                fl = fr;
                il=0;
                iu++;
                if(iu>=2){
                    fu/=2;
                }

            } else{
                ea = T(0); //NO HAY ERROR
                xr = (fl ==T(0)) ? xl : xu;
            }

            if(ea < eps){ // Si el error aproximado es menor que el epsilon entoneces llego a la aproximacion esperada.
                cout<< xr << endl;
                return xr;
            }

        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

}

#endif

