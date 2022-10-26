// IS lab 3.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Algorithms.h"

double f(double x, double y) {
    return 1.0 * (std::pow(x - 2, 2) + std::pow(y + 1, 2));
}

double neg_f(double x, double y) {
    return -1.0 * (std::pow(x - 2, 2) + std::pow(y + 1, 2));
}

int main()
{
    std::cout << "Hello World!\n";
    srand(time(NULL));
    double (*func)(double, double);
    func = f;
    auto alg = Annealing<double>(1000000.0, 0.000005, 0.25, func);
    std::cout << "Annealing algorithm:\n";
    auto res = alg.find(-500.0, 500.0, ExtremumType::min);
    std::cout << res.first << " " << res.second << std::endl;

    /*func = neg_f;
    auto alg2 = Genetic<double>(100, 1000, func);
    std::cout << "Genetic algorithm:\n";
    res = alg2.find(-5.0, 5.0, ExtremumType::max);
    std::cout << res.first << " " << res.second << std::endl;*/

    func = f;
    auto alg3 = HillClimbing<double>(0.01, func);
    std::cout << "Hill climbing greedy algorithm:\n";
    res = alg3.find(-5.0, 5.0, ExtremumType::min);
    std::cout << res.first << " " << res.second << std::endl;
    std::cout << func(res.first, res.second) << std::endl;
    std::cout << "Hill climbing first algorithm:\n";
    res = alg3.find_first(-5.0, 5.0, ExtremumType::min);
    std::cout << res.first << " " << res.second << std::endl;
    std::cout << func(res.first, res.second) << std::endl;
    std::cout << "Hill climbing first with restart algorithm:\n";
    res = alg3.find_first_with_restart(-5.0, 5.0, ExtremumType::min);
    std::cout << res.first << " " << res.second << std::endl;
    std::cout << func(res.first, res.second) << std::endl;
    std::cout << "Hill climbing stochastic algorithm:\n";
    res = alg3.find_stochastic(-5.0, 5.0, ExtremumType::min);
    std::cout << res.first << " " << res.second << std::endl;
    std::cout << func(res.first, res.second) << std::endl;
    std::cout << "Hill climbing ray algorithm:\n";
    res = alg3.find_ray(-5.0, 5.0, ExtremumType::min);
    std::cout << res.first << " " << res.second << std::endl;
    std::cout << func(res.first, res.second) << std::endl;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
