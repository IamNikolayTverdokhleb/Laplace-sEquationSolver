//
// Created by Коля on 2019-05-29.
//

#include "Membrane.h"
using namespace std;
/*
 ****************               Вспомогательные функции               ****************
 */
vector<double> Membrane::progonka(vector<double> a, vector<double> b, vector<double> c, vector<double> d, int node){
    vector<double> alpha(node);
    vector<double> betta(node);
    vector<double> Solve(node);

    alpha[0] = 0.0;
    betta[0] = 0.0;

    alpha[1] = c[0] / b[0];
    betta[1] = d[0] / b[0];

    for (int i = 1; i < node-1; i++)
    {
        alpha[i + 1] = c[i] / (b[i] - a[i] * alpha[i]);
        betta[i + 1] = (d[i] + a[i] * betta[i]) / (b[i] - a[i] * alpha[i]);
    }

    Solve[node-1] = (d[node - 1] + (a[node - 1] * betta[node-1])) / (b[node - 1] - (a[node - 1] * alpha[node-1]));


    for (int i = node - 2; i >= 0; i--)
    {
        Solve[i] = alpha[i + 1] * Solve[i + 1] + betta[i + 1];
    }
    return Solve;
}
double Membrane::f(double t) {
    switch(test_index){
        case(1):{
            return 0.0;
        }
        case(2):{
            return 0.0;
        }
        case(3):{
            return -4.0;
        }
        case(4):{
            return 0.0;
        }
    }

}
double Membrane::g10(double t, double Mirror, double x1, double x2,int index) { /*Граничное условие при x1 = 0 */
    switch(test_index){
        case(1):{
            return 1.0;
        }
        case(2):{
            return 1.0 + x2;
        }
        case(3):{
            return x2*x2;
        }
        case(4):{
            return 0.0;
        }
    }
}
double Membrane::g1n(double t, double Mirror, double x1, double x2,int index) { /*Граничное условие при x1n = 0 */
    switch(test_index){
        case(1):{
            return 1.0;
        }
        case(2):{
            return 1.0 + x2;
        }
        case(3):{
            return 1.0 + x2*x2;
        }
        case(4):{
            return cos(x1)*sin(1)*exp(-2.0*t);
        }
    }
}
double Membrane::g20(double t, double Mirror, double x1, double x2,int index) { /*Граничное условие при x2 = 0 */
    switch(test_index){
        case(1):{
            return 1.0;
        }
        case(2):{
            return Mirror + index*tau;
        }
        case(3):{
            return Mirror;
        }
        case(4):{
            return sin(x2)*exp(-2.0*t);
        }
    }
}
double Membrane::g2n(double t, double Mirror, double x1, double x2,int index) { /*Граничное условие при x2n = 0 */
    switch(test_index){
        case(1):{
            return 1.0;
        }
        case(2):{
            return  Mirror + index*tau;
        }
        case(3):{
            return Mirror + index*tau;
        }
        case(4):{
            return cos(1.0)*sin(x2)*exp(-2.0*t);
        }
    }
}
void Membrane::description(){
    cout << "Длина мембраны по x1: " << L1 << endl;
    cout << "Длина мембраны по x2: " << L2 << endl;
    cout << "Шаг  по x1: " << h1 << endl;
    cout << "Шаг  по x2: " << h2 << endl;
    cout << "Шаг  по времени: " << tau << endl;
    cout << "Решается тест: " << test_index << endl;
}

void Membrane::fileOutput(tMatrix Solution){
    ofstream LayerOutput("../Results/Membrane.txt");
    for(int i = 0; i < node2; i++) {
        for (int j = 0; j < node1; j++) {
            LayerOutput << i * h1 << "\t\t" << j*h2 << "\t\t" << Solution[i][j] << endl;
        }
    }
}

void Membrane::fileOutputError(tMatrix Solution){
    ofstream LayerOutput("/Users/kola/Desktop/Uni/3Course/6Sem/6SemMethods/LAB4/Work/Results/VisualizeError.txt");
    cout << "Выполнено." << endl;
    for(int i = 0; i < node2; i++) {
        for (int j = 0; j < node1; j++) {
            LayerOutput << i * h1 << "\t\t" << j*h2 << "\t\t" << Solution[i][j] << endl;
        }
    }
}

void Membrane::fileOutputAnaliticalSolution(){
    auto s = to_string(test_index);
    tMatrix analitical_solution = new double* [node1];
    for (int i = 0; i < node2; i++) {
        analitical_solution[i] = new double[node2];
    }
    ofstream LayerOutput("../Results/MembraneAnaliticalSolution" + s + ".txt");
    switch(test_index) {
        case (1): {
            for (int i = 0; i < node2; i++) {
                for (int j = 0; j < node1; j++) {
                    analitical_solution[i][j] = 1.0;
                    LayerOutput << i * h1 << "\t\t" << j * h2 << "\t\t" << analitical_solution[i][j] << endl;
                }
            }
            break;
        }
        case (2): {
            for (int i = 0; i < node2; i++) {
                for (int j = 0; j < node1; j++) {
                    analitical_solution[i][j] = i*h1 + 1.0;
                    LayerOutput << i * h1 << "\t\t" << j * h2 << "\t\t" <<  analitical_solution[i][j] << endl;
                }
            }
            break;
        }
        case (3): {
            for (int i = 0; i < node2; i++) {
                for (int j = 0; j < node1; j++) {
                    analitical_solution[i][j] = i*i*h1*h1 + j*j*h2*h2;
                    LayerOutput << i * h1 << "\t\t" << j * h2 << "\t\t" <<  analitical_solution[i][j] << endl;
                }
            }
            break;
        }
        case (4): {
            for (int i = 0; i < node2; i++) {
                for (int j = 0; j < node1; j++) {
                    analitical_solution[i][j] = i*i*h1*h1 + j*j*h2*h2;
                    LayerOutput << i * h1 << "\t\t" << j * h2 << "\t\t" <<  analitical_solution[i][j] << endl;
                }
            }
            break;
        }
    }
}

void Membrane::errorCheck(tMatrix Solution, double t) {
    double error_norm = 0.0;
    double temp_error_norm;
    switch(test_index){
        case(1): {
            for (int i = 0; i < node2; i++) {
                temp_error_norm = 0.0;
                for (int j = 0; j < node1; j++) {
                    temp_error_norm += abs(Solution[i][j] - 1.0);
                    if (temp_error_norm > error_norm)
                    error_norm = temp_error_norm;
                }
            }
            break;
        }

        case(2):{

            for (int i = 0; i < node2; i++) {
                temp_error_norm = 0.0;
                for (int j = 0; j < node1; j++) {
                    temp_error_norm += abs(Solution[i][j] - j*h1 - 1.0);
                    if (temp_error_norm > error_norm)
                        error_norm = temp_error_norm;
                }
            }
            break;
        }
        case(3):{
            for (int i = 0; i < node2; i++) {
                temp_error_norm = 0.0;
                for (int j = 0; j < node1; j++) {
                    temp_error_norm += abs(Solution[i][j] - i*i*h1*h1 - j*j*h2*h2);
                    if (temp_error_norm > error_norm)
                        error_norm = temp_error_norm;
                }
            }
            break;
        }
        case(4):{
            tMatrix visualizeError = new double*[node1];
            for (int i = 0; i < node2; i++)
            {
                visualizeError[i]=new double[node2];
                for (int j=0; j < node1; j++)
                {
                    visualizeError[i][j]=0.0;
                }
            }
            
            for (int i = 1; i < node2-1; i++) {
                temp_error_norm = 0.0;
                for (int j = 1; j < node1-1; j++) {
                    visualizeError[i][j] = abs(Solution[i][j] - cos(j*h1)*sin(i*h2)*exp(-2.0*t));
                    temp_error_norm = abs(Solution[i][j] - cos(j*h1)*sin(i*h2)*exp(-2.0*t));
                    if (temp_error_norm > error_norm)
                        error_norm = temp_error_norm;
                }
             
            }
            fileOutputError(visualizeError);
            break;
        }
    }
    cout << "Норма разности точного и численного решений равна " << error_norm <<endl;
}
void Membrane::assemble_boundary_conditions(){
    switch(test_index) {
        case (1):
            this->bottom_flow_index = 0;
            this->top_flow_index = 0;
            this->left_flow_index = 0;
            this->right_flow_index = 0;
            break;
        case (2):
            this->bottom_flow_index = 0;
            this->top_flow_index = 0;
            this->left_flow_index = 1;
            this->right_flow_index = 1;
            break;
        case (3):
            this->bottom_flow_index = 0;
            this->top_flow_index = 0;
            this->left_flow_index = 1;
            this->right_flow_index = 1;
            break;

        case (4):
            this->bottom_flow_index = 0;
            this->top_flow_index = 0;
            this->left_flow_index = 0;
            this->right_flow_index = 0;
            break;
    }
}

