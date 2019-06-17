//
// Created by Коля on 2019-05-29.
//
#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#pragma once
#ifndef WORK_MEMBRANE_H
#define WORK_MEMBRANE_H

using namespace std;
typedef double *tVector ;
typedef double** tMatrix;
struct ParametrsStruct{
    double tau, h1, h2, L1,L2, t0;
    int node1, node2, test_index;
};

class Membrane {
public:
    Membrane(ParametrsStruct parametrsStruct) {
        this->tau = parametrsStruct.tau;
        this->h1 = parametrsStruct.h1;
        this->h2 = parametrsStruct.h2;
        this->L1 = parametrsStruct.L1;
        this->L2 = parametrsStruct.L2;
        this->t0 = parametrsStruct.t0;
        this->node1 = parametrsStruct.node1;
        this->node2 = parametrsStruct.node2;
        this->test_index = parametrsStruct.test_index;
    }
    ~Membrane(){};
    vector<double> progonka(vector<double> a, vector<double> b, vector<double> c, vector<double> d, int node);
    double f(double t);
    double g10(double t, double Mirror, double x1, double x2, int index); /*g - функции для определения ГУ первого рода */
    double g1n(double t, double Mirror, double x1, double x2, int index);
    double g20(double t, double Mirror, double x1, double x2, int index);
    double g2n(double t, double Mirror, double x1, double x2, int index);
    void description();
    void assemble_boundary_conditions();
    void fileOutput(tMatrix Solution);
    void fileOutputError(tMatrix Solution);
    void fileOutputAnaliticalSolution();
    void errorCheck(tMatrix Solution, double t);
/*
 ****************               Задаем начальные условия                ****************
 */
    tMatrix initialLayer(){
        /*Создание матрицы начального слоя*/
        tMatrix initial_Layer = new double*[node1];
        for (int i = 0; i < node2; i++)
        {
            initial_Layer[i]=new double[node2];
            for (int j=0; j < node1; j++)
            {
                initial_Layer[i][j]=0.0;
            }
        }
        /*Заполнение матрицы начального слоя в зависимости от теста*/
        switch(test_index) {
            case (1): {
                for (int i = 0; i < node2; i++) {
                    for (int j = 0; j < node1; j++) {
                        initial_Layer[i][j] = 0.0;
                    }
                }
                break;
            }
            case (2): {
                for (int i = 0; i < node2; i++) {
                    for (int j = 0; j < node1; j++) {
                        initial_Layer[i][j] = 1.0;
                    }
                }
                break;
            }
            case (3): {
                for (int i = 0; i < node2; i++) {
                    for (int j = 0; j < node1; j++) {
                        initial_Layer[i][j] = 1.0;
                    }
                }
                break;
            }
            case (4): {
                for (int i = 0; i < node2; i++) {
                    for (int j = 0; j < node1; j++) {
                        initial_Layer[i][j] = cos(j*h1)*sin(i*h2);
                    }
                }
                break;
            }
        }
        return initial_Layer;
    }

    /*
   ****************             Реализация расчета               ****************
   */
    void membraneCalculations()
    {
        int layer_num = 1;
        double t = t0 + tau,
                Lambda1 = 0.0,
                Lambda2 = 0.0,
                dif_matrix_norm = 0.0;

        vector<double> a1(node1);
        vector<double> b1(node1);
        vector<double> c1(node1);
        vector<double> d1(node1);
        vector<double> a2(node2);
        vector<double> b2(node2);
        vector<double> c2(node2);
        vector<double> d2(node2);
        vector<double> temp1(node2);
        vector<double> temp2(node2);
        tMatrix y_current = new double* [node1];
        tMatrix y_half = new double* [node1];
        tMatrix y_next = new double* [node1];
        for (int i = 0; i < node2; i++){
            y_current[i]=new double[node2];
            y_half[i]=new double[node2];
            y_next[i]=new double[node2];
            for (int j=0; j < node1; j++)
            {
                y_current[i][j]=0.0;
                y_half[i][j]=0.0;
                y_next[i][j]=0.0;
            }
        }
        y_current = initialLayer();

        do{
            t += 0.5*tau;
            /*! Первая прогонка по направлению j*/
            for(int i = 1; i < node2 - 1; i++){
                for(int j = 1; j < node1 - 1; j++){
                    a1[j] = 1.0/(h1*h1);
                    b1[j] = 2.0*((1.0/(h1*h1)) + 1.0/tau);
                    c1[j] = 1.0/(h1*h1);
                    Lambda2 = (y_current[i][j+1] - 2.0*y_current[i][j] + y_current[i][j-1])/(h2*h2);
                    d1[j] = (2.0/tau)*y_current[i][j] + f(t) + Lambda2;
                }
                a1[0] = 0.0;
                b1[0] = 1.0;
                c1[0] = left_flow_index;
                d1[0] = g20(t, 0.0, 0.0,i*h2, 1);
                a1[node1 - 1] = left_flow_index;
                b1[node1 - 1] = 1.0;
                c1[node1 - 1] = 0.0;
                d1[node1 - 1] = g2n(t, 0.0,L2, i*h2, 1);
                temp1 = progonka(a1, b1, c1, d1, node1);
                for(int j = 0; j < node1; j++){
                    y_half[i][j] = temp1[j];
                }
            }
            for(int j = 0; j < node1; j++){
                y_half[0][j] = g10(t, y_half[1][j], j*h2,0.0,0);
                y_half[node2 - 1][j] = g1n(t, y_half[node2 - 2][j], j*h2, L2,0);
            }

            t += 0.5*tau;

            /*! Вторая прогонка по направлению i*/

            for(int j = 1; j < node2 - 1; j++){
                for(int i = 1; i < node1 - 1; i++){
                    a2[i] = 1.0/(h2*h2);
                    b2[i] = 2.0*((1.0/(h2*h2)) + 1.0/tau);
                    c2[i] = 1.0/(h2*h2);
                    Lambda1 = (y_half[i+1][j] - 2.0*y_half[i][j] + y_half[i-1][j])/(h2*h2);
                    d2[i] = (2.0/tau)*y_half[i][j] + f(t) + Lambda1;
                }
                a2[0] = 0.0;
                b2[0] = 1.0;
                c2[0] = bottom_flow_index;
                d2[0] = g10(t,y_next[1][j],0.0,j*h1, 0);
                a2[node2 - 1] = top_flow_index;
                b2[node2 - 1] = 1.0;
                c2[node2 - 1] = 0.0;
                d2[node2 - 1] = g1n(t,y_next[node2 - 1][j],j*h2,0.0, 0);
                temp2 = progonka(a2, b2, c2, d2, node2);
                for(int i = 0; i < node1; i++){
                    y_next[i][j] = temp2[i];
                }
            }
            for(int i = 0; i < node1; i++){
                y_next[i][0] = g20(t, y_next[i][1], 0.0,i*h2, 1);
                y_next[i][node2 - 1] = g2n(t, y_next[i][node2 - 2], 0.0,i*h2,1);
            }

            layer_num++;

            /*Копирование решения для перехода на новый слой и вычисление матричных норм*/
            dif_matrix_norm = 0.0;
            double temp_difference_norm;
            for(int i = 0; i < node2; i++){
                temp_difference_norm = 0.0;
                for(int j = 0; j < node1; j++) {
                    temp_difference_norm += abs(y_next[i][j] - y_current[i][j]);
                    y_current[i][j] = y_next[i][j];
                    if (temp_difference_norm > dif_matrix_norm)
                        dif_matrix_norm = temp_difference_norm;
                }
            }

            cout << "На " << layer_num << " слое норма разности решений равна " << dif_matrix_norm << endl;

            //if (dif_matrix_norm < tau*eps) {
            if(abs(t - 1.0) < eps){
                fileOutput(y_next);
                //cout << "Решение вышло на стационар за " << t << " секунды." << endl;
                cout << "Погрешность в момент времени " << t << endl;
                errorCheck(y_next, t);
                break;
            }
        } while(true);

    }
    /*
   ****************               Запуск алгоритма               ****************
   */
    void run(){
        description();
        assemble_boundary_conditions();
        membraneCalculations();
        fileOutputAnaliticalSolution();
    }

private:
    double tau, h1, h2, L1, L2, t0;
    int node1, node2, test_index, left_flow_index, right_flow_index,top_flow_index,bottom_flow_index;
    const  double eps = pow(10,-3);
};


#endif //LAB4_MEMBRANE_H
