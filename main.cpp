#include <iostream>
#include "Membrane.h"
#include "Scripts.h"
#include <string>
#include <vector>

using namespace std;
int main()
{
    // Выбираем параметры расчетов
    double tau = 0.003125  ,
            h1 = 0.003125   ,
            h2 = 0.003125  ,
            L1 = 1.0,
            L2 = 1.0,
            t0 = 0.0;

    int node1 = int(L1/h1) + 1,
            node2 = int(L2/h2) + 1,
            test_index = 4;
    /*!
     * @param test_index
     * 1 - температура равна 1 на краях мембраны
     * 2 - второй тест из методички
     * 3 - третий тест из методички
     * 4 - уравение теплопроводности для анализа порядка аппроксимации
     */
    Scripts::clearAll();
    struct ParametrsStruct parametrsStruct = {tau, h1, h2, L1, L2, t0, node1, node2, test_index};
    Membrane membrane(parametrsStruct);
    membrane.run();
    Scripts::BuildMesh();
    return 0;
}
