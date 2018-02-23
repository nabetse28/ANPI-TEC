format long


#Funcion factorial 
function f = fact(num)
  f=1;
  if num<=1
    f=1;
  endif
  
  for i=1:num
    f=f*i;
  endfor;
  
endfunction


#Funcion principal del programa
function analisis(x,n)
  es=(x*(10^(2-n)));
  vr = e^(x);
  estimado=0;
  anterior=estimado;
  erpa = 2*es;
  erpv = 0;
  cont=0;
  while erpa>es
    estimado=estimado+(x^(cont))/(fact(cont));
    erpv = ((vr-estimado)/(vr))*100;
    erpa = ((estimado-anterior)/(estimado))*100;
    anterior=estimado;
    cont++;   
  endwhile
  disp("Termino");
  disp(cont+1);
  disp("Estimado:");
  disp(estimado);
  disp("Error Relativo Porcentual Verdadero:");
  disp(erpv);
  disp("Error Relativo Porcentual Aproximado:");
  disp(erpa);
  
  
endfunction


#Menu del Programa
function start() 
  n = input(">>> Ingrese la cantidad de cifras significativas: ");
  if(abs(mod(n,1))!=0)
    disp(">>> Por favor ingrese numero entero");
    start();
  elseif(n>15)
    disp("Ingrese un numero menor o igual a 15");
    start();
  else
    x = input(">>> Ingrese el valor de x para la funcion: ");
    disp("--------------------------------------");
    analisis(x,n);
  end
    
endfunction

start();
