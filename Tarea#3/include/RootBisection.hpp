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

#include "Exception.hpp"

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

using namespace std;

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], using the bisection method.
     *
     * @param funct a std::function of the form "T funct(T x)"
     * @param xl xl interval limit
     * @param xu xu interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps = std::numeric_limits<T>::epsilon()) {

        const T maxi= numeric_limits<T>::digits*2;
        T xr=xl; //Se define un valor inicial para la raiz
        T xl_value=funct(xl); //Se calcula el valor inicial



        T error_aprox=T(); //Se iniciliza un error


        if(xl > xu){
            throw anpi::Exception("El intervalo esta invertido");
        }

        if((funct(xl)*funct(xu)>T(0))){
            throw anpi::Exception("Unenclosed root");
        }



        for (int i=maxi; i>0; --i){

            T xr_antes(xr); //se guarda el valor anterior de la posicion raiz estimada

            xr=(xl+xu)/2; //Se calcula un nuevo valor para llegar a la raiz

            T funcion_raiz=funct(xr); // se determina una nueva raiz

            if (abs(xr) > numeric_limits<T>::epsilon()) { //Se asegura de no dividir por cero
                error_aprox=abs((xr-xr_antes)/xr)*100; // se calcula un error relativo
            }

            T condicion = xl_value*funcion_raiz; // se define la condicion para tomar la siguiente accion

            if (condicion > 0) { //Se mueve el limite inferior
                xl=xr;
                xl_value=funcion_raiz;

            }
            else if (condicion < 0) { //Se mueve el limite superior
                xu=xr;
            }
            else {
                error_aprox = 0; // Se espera tener cero en uno de los limites
                xr = (abs(xl_value) < numeric_limits<T>::epsilon())
                     ? xl : xr;
            }

            if (error_aprox < eps){
                cout<<xr<<endl;
                return xr; //Declaracion de eps, cuando el error sea suficientemente bajo
            }
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN(); // no encontro la raiz!
    }

}


#endif

