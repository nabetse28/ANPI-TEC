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

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP


using namespace std;

namespace anpi {
  
  /**
   * Find a root of the function funct looking for it starting at xi
   * by means of the secant method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial position
   * @param xii second initial position 
   *
   * @return root found, or NaN if no root could be found
   */
  template<typename T>
  T rootSecant(const std::function<T(T)>& funct,T xi,T xii,const T es= sqrt(numeric_limits<T>::epsilon()) ) {

    const T maxi= numeric_limits<T>::digits;

    T xl, rts;
    T fl = funct(xi); // Se evalua xi en la funcion
    T f= funct(xii); // Se evalua xii en la funcion
    T temp; // Variable temporal para hacer un swap de variables

    if(abs(fl)<abs(f)){ // Si los datos estan invertidos
      rts = xi; // Esta es la que deberia de ser la raiz
      xl = xii; // Este es el x mas grande

      temp = fl; // Se hace el swap
      fl = f;
      f = temp;


    }else{
      xl = xi;
      rts = xii;

    }

    for(int i=0; i<maxi;i++){ // Se comienza a iterar hasta que el error sea menor que epsilon o que la funcion evaluada en el nueva raiz conveja
      T dx = (xl-rts)*f/(f-fl);
      xl = rts;
      fl=f;
      rts+=dx;
      f=funct(rts);

      if( abs(dx) < es || f == 0.0){ // Condicion de parada
        cout<<rts<< endl;
        return rts;

      }

    }
    
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

