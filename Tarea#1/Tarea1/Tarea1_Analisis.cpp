#include <iostream>
#include <math.h>
using namespace std;

template <class T>
void problema1(){
    int n;
    T x, error, vprox, vreal;
    vprox=0;
    cout<<">>> Ingrese un valor para n: "<< endl;
    cin>>n;
    cout<<">>> Ingrese un valor para x: "<< endl;
    cin>>x;

    //Este es el valor real
    vreal = n*x;

    //Para calcular el valor aproximado
    for(int i=0;i<n;i++){
        vprox+=x;
    }

    error = ((vreal-vprox)/(vreal))*100;
    cout<<"----------------------------------------------"<<endl;
    cout<<"Valor Real: "<<vreal<<endl;
    cout<<"Valor Aproximado: "<<vprox<<endl;
    cout<<"El Error Verdadero es "<<error<<"%"<<endl;
    cout<<"\n"<<endl;

}

template <class T>
void problema2(){
    T a, b, c, x1_T, x2_T, x1_A, x2_A, error_x1, error_x2,det;

    cout<<"Ingrese el valor de a: "<<endl;
    cin>>a;
    cout<<"Ingrese el valor de b"<<endl;
    cin>>b;
    cout<<"Ingrese el valor de c"<<endl;
    cin>>c;

    det= (b*b)-(4*a*c);
    if(det<0){
        cout<<"El determinante es menor a cero por lo que no se pueden realizar operaciones con numeros reales"<<endl;

    }else{
        //Para calcular las raices de forma tradicional

        x1_T = ((-b)+sqrt(b*b-(4*a*c)))/(2*a);
        x2_T = ((-b)-sqrt(b*b-(4*a*c)))/(2*a);

        //Para calcular las raices de forma alternativa

        x1_A = (-2*c)/(b+sqrt(b*b-(4-a*c)));
        x2_A = (-2*c)/(b-sqrt(b*b-(4-a*c)));

        //Para calcular los errores
        error_x1 = ((x1_T-x1_A)/(x1_T))*100;
        error_x2 = ((x2_T-x2_A)/(x2_T))*100;

        cout<<"Las raices tradicionales con a= "<<a<<" b= "<<b<<" c= "<<c<<" son:"<<endl;
        cout<<"X1= "<<x1_T<<endl;
        cout<<"X2= "<<x2_T<<endl;

        cout<<"Las raices alternativas con a= "<<a<<" b= "<<b<<" c= "<<c<<" son:"<<endl;
        cout<<"X1= "<<x1_A<<endl;
        cout<<"X2= "<<x2_A<<endl;

        cout<<"Error de X1 = "<<error_x1<<"%"<<endl;
        cout<<"Error de X2 = "<<error_x2<<"%"<<endl;

    }



}


int main(){
    int x,p;
    bool flag = true;

    while (flag){
        cout<<">>> Ingrese 2 o 3 si se quiere el problema 2 o 3, si quiere salir del programa ingrese 0: "<<endl;
        cin>>x;
        if (x==2){
            cout<<">>> Escoja 1 o 2 para precision simple o doble: "<<endl;
            cin>>p;
            if(p==1){
                cout<<"Se va a ejecutar el programa con Presicion Simple"<<endl;
                problema1<float>();
            }else if(p==2) {
                cout << "Se va a ejecutar el programa con Presicion Simple" << endl;
                problema1<float>();
            }else if(p==0){
                flag=true;
            }else{
                main();
            }
        }else if(x==3){
            cout<<">>> Escoja 1 o 2 para precision simple o doble: "<<endl;
            cin>>p;
            if(p==1){
                cout<<"Se va a ejecutar el programa con Presicion Simple"<<endl;
                problema2<float>();
            }else if(p==2) {
                cout << "Se va a ejecutar el programa con Presicion Simple" << endl;
                problema2<float>();
            }else if(p==0){
                flag=true;
            }else{
                main();
            }
        }else if(x==0){
            flag=true;
        }else{
            main();
        }
    }
    return 0;
}