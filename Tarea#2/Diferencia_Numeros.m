# -----------------------------------------
# Autor: Esteban Herrera Vargas
# Curso: Analisis Numerico para Ingenieria
# Tarea #2
# Profesor: Jose Pablo Alvarado
# Universidad: Instituto Tecnologico de Costa Rica
# Carrera: Ingenieria en Computadores
# -----------------------------------------



function menu()
  in=input("Ingrese 1 para correr la funcion con los parametros descritos o 2 para escoger los parametros: ");
   
  
  if(in==1)
    disp("Parametros a usar h = 1, x = 1 e iteraciones = 10000");
    h = 1;   
    x = 1;   
    n = 10000;
    start(h,x,n);
  elseif(in==2)
    h = input("Ingrese el valor de h: ");
    x = input("Ingrese el valor de x: ");
    n = input("Ingrese el valor de las iteraciones a realizar: ");
    start(h,x,n);
  else
    disp("Este dato no es valido");
  endif
  
endfunction

function start(h,x,n)

  funcion = @(x)(0.3*x*x*x*x - 0.15*x*x);


  funcionDev = @(x)(0.3*4*x*x*x - 2*0.15*x);


  difCentrada = @(x,h)((funcion(x+h)-funcion(x-h))/(2*h));


  difAdelante = @(x,h)((funcion(x+h)-funcion(x))/(h));


  difAtras = @(x,h)((funcion(x)-funcion(x-h))/(h));

  
  
  hpasos = zeros(size(n));    
  difCentradaPasos = zeros(size(n));    
  difAdelantePasos = zeros(size(n));   
  difAtrasPasos = zeros(size(n));   
  for i=1:n
    fdev = funcionDev(x);       
    difCentradaPasos(i) = (abs(difCentrada(x,h) - fdev))/fdev;    
    difAdelantePasos(i) = (abs(difAdelante(x,h) - fdev))/fdev;     
    difAtrasPasos(i) = (abs(difAtras(x,h) - fdev))/fdev;     
    hpasos(i) = h;
    h = (h)/(10**((1)/i));    
  endfor
  
 
  loglog(hpasos,difCentradaPasos,'-b');
  hold on
  loglog(hpasos,difAdelantePasos,'-r');
  hold on
  loglog(hpasos,difAtrasPasos,'-g');
  hold off
  legend("Centrada","Adelante","Atras");
  pause();
endfunction



menu();