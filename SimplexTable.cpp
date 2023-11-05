#include "SimplexTable.h" // Включаем заголовочный файл, который объявляет ваши классы
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Simplex {

// 2 СЛОЙ
    ConstraintSystem::ConstraintSystem(std::vector<std::vector<double>> matrix) : matrix(matrix) {
        // Инициализируем матрицу ограничений
    }
// 3 СЛОЙ
    void ConstraintSystem::calculate(int index_for_first, std::vector<double> k) {
        // Пересчет матрицы системы ограничений
        // для index_for_first домножаем на k[index_for_first]
        for (size_t i = 0; i < matrix[index_for_first].size(); i++) {
            matrix[index_for_first][i] /= k[index_for_first];
        }
        // для остальные вычитаем полученную первую строку помноженную на k[index_for_first]
        for (size_t i = 0; i < matrix.size(); i++) {
            if (i != index_for_first) {
                for (size_t j = 0; j < matrix[i].size(); j++) {
                    matrix[i][j] -= matrix[index_for_first][j]*k[i]; // для index_for_first-стобца это 1*k[i] по идее..
                }
            }
        }
    }
// 2 СЛОЙ
    TargetFunc::TargetFunc(std::vector<double> c, Optimization target) : c(c), target(target) {
        // Инициализируем целевую функцию
        int size_of_x = c.size(); 
        this->x = std::vector<double>(size_of_x, 0.0);
    }
// ? СЛОЙ
    double TargetFunc::sum() {
        // Реализация метода sum
        double result = 0.0;
        for (double value : x) {
            result += value;
        }
        return result;
    }

// 1 СЛОЙ
    SimplexTable::SimplexTable(ConstraintSystem constraints, TargetFunc Lx)
        : constraints(constraints), Lx(Lx) {
        // Инициализация полей класса SimplexTable
        this->x_plan = std::vector<double>(Lx.c.size(), 0.0);
        this->z = std::vector<double>(Lx.c.size());
        this->delta = std::vector<double>(Lx.c.size());

    }

// ? СЛОЙ
    double SimplexTable::get_current_Lx() {
        // Реализация метода get_current_Lx        
        return Lx.sum();
    }


// 2 СЛОЙ
    BasisCoefficient SimplexTable::get_currentA_row(int literal_numb_A) { // Находим ЧТО заменить (ЭТО заменяем)
        // Реализация метода get_currentA_row
        double dteta = std::numeric_limits<double>::max(); // MAXFLOAT;
        
        BasisCoefficient teta = BasisCoefficient();
        teta.C_basis = std::numeric_limits<double>::max(); // MAXFLOAT;
        teta.literal_numb_A = 0; // NAN;
        for (int i = 0; i<constraints.matrix.size(); ++i ) {
            dteta = constraints.matrix[i][0]/constraints.matrix[i][literal_numb_A];
            if (dteta<teta.C_basis) {
                teta.C_basis = dteta;
                teta.literal_numb_A = i;
            }
        }
        return teta;
    }

// 2 СЛОЙ
    BasisCoefficient SimplexTable::get_newA_colRow() { // Находим НА ЧТО заменить (ЭТИМ заменяем)
        // Реализация метода get_currentA_col

        double min_delta = std::numeric_limits<double>::max();
        int min_delta_index = -1;

        for (size_t i = 0; i < delta.size(); ++i) {
            if (delta[i] < min_delta) {
                min_delta = delta[i];
                min_delta_index = i;
            }
        }

        return BasisCoefficient{ min_delta_index, min_delta };
    }
// 2 СЛОЙ
    void SimplexTable::change_beginer_basis(BasisCoefficient currentA, BasisCoefficient newA) {
        // Замена элемента в базисе
        BeginerBasis[currentA.literal_numb_A] = newA;

        // Пересчет матрицы ограничений
        std::vector<double> k(constraints.matrix.size(), 0.0);

        for (size_t i = 0; i < constraints.matrix.size(); ++i) {
            k[i] = constraints.matrix[i][currentA.literal_numb_A] / constraints.matrix[currentA.literal_numb_A][currentA.literal_numb_A];
        }

        constraints.calculate(currentA.literal_numb_A, k); // помним что currentA теперь newA

        // Обновление опорного плана
        for (size_t i = 0; i < constraints.matrix.size(); ++i) {
            x_plan[i] = constraints.matrix[i][0]; // Обновление свободных членов в опорном плане
        }
        for (size_t i = 0; i < Lx.x.size(); ++i){
            for (BasisCoefficient b : BeginerBasis){
                // если литерал А имеется в начальном наборе показателей то учитываем его
                if (i == b.literal_numb_A){
                    Lx.x[i] = x_plan[i];
                }
                // иначе он не нужен в нашем итоговом плане
            }
        }
    }
// 2 СЛОЙ
    bool SimplexTable::check_delta() {
        // Реализация метода check_delta
        if (Lx.target == TargetFunc::MAXIMUM) {
            // проверка что все дельты неотрицательны
            for (double d: delta) {
                if (d < 0) return false; // найден отрицательный компонент - плане не оптимален
            }
            return true; // все элементы неотрицательны - найден оптимальный план
        }
        else if (Lx.target == TargetFunc::MINIMUM) {
            // проверка что все дельты неположительны
            for (double d: delta) {
                if (d>0) {
                    return false;
                }
            }
            return true;
        }
        else {
            // цель неопределена:
            throw std::logic_error("Цель неопределена");
        }

    }
    
// 2 СЛОЙ
    std::vector<double> SimplexTable::calculate_deltas() {
        std::vector<double> deltas(Lx.c.size(), 0.0);
        for (size_t i = 0; i < Lx.c.size(); ++i) {
            double delta = 0.0;
            for (size_t j = 0; j < BeginerBasis.size(); ++j) {
                std::cout << constraints.matrix[j][i + 1] << " ";
                delta += BeginerBasis[j].C_basis * constraints.matrix[j][i + 1];
            }
            deltas[i] = Lx.c[i] - delta;
        }
        return deltas;
    }

// 1 СЛОЙ - ТОЧКА ВХОДА
    void SimplexTable::iterating() {
        get_beginer(); // Сначала ищем начальный базис и план
        int iteration = 0;
        while (iteration < MAX_ITERATIONS) {
            delta = calculate_deltas();

            // Проверка на оптимальность
            if (check_delta()) {
                break; // Оптимальное решение найдено
            }

            auto leaving = get_newA_colRow(); // выводимый базис

            
            // Проверка на неограниченность решения
            bool unbounded = true;
            for (size_t i = 0; i < constraints.matrix.size(); ++i) {
                if (constraints.matrix[i][leaving.literal_numb_A] > 0) { 
                    unbounded = false; 
                    break; 
                    // если есть хотя бы один отрицательный - все нормально
                }
            }

            if (unbounded) {
                std::cout << "Решение неограничено." << std::endl;
                return;
            }

            auto entering = get_currentA_row(leaving.literal_numb_A);
            change_beginer_basis(leaving, entering);
            ++iteration;

            if (iteration == MAX_ITERATIONS) {
                std::cout << "Maximum iterations reached. No optimal solution found." << std::endl;
                return;
            }
        }

        
    }

// 2 СЛОЙ
    void SimplexTable::get_beginer(){
        // Найти начальный базис (пример: первые M переменных)
        BeginerBasis.clear();
        x_plan.clear();
        for (int i = 0; i < constraints.matrix[0].size(); i++) { //i-й столбец |||...|
            BasisCoefficient basis_coeff;
            basis_coeff.literal_numb_A = i; // Присваиваем номер переменной в базисе

            // Проверяем, является ли столбец единичным
            bool is_basic = true;  // Изначально предполагаем, что столбец не является единичным
            for (int j = 0; j < constraints.matrix.size(); j++) { // j-я строка __
                if (constraints.matrix[j][i] != 1.0 && constraints.matrix[j][i] != 0.0) {
                    is_basic = false;
                }
            }
            int one_count = 0; // Счетчик единиц в столбце
            for (int j = 0; j < constraints.matrix.size(); j++) {
                if (constraints.matrix[j][i] == 1.0) {
                    one_count++;
                    if (one_count > 1) {
                        is_basic = false;
                        break;
                    }
                } 
            }

            if (is_basic) {
                if (basis_coeff.literal_numb_A < Lx.c.size()){
                    basis_coeff.C_basis = Lx.c[basis_coeff.literal_numb_A];
                }
                basis_coeff.C_basis = 0.0; // Начальное значение коэффициента в базисе
            }

            BeginerBasis.push_back(basis_coeff);

            // Определить начальный опорный план на основе свободных членов
            x_plan.push_back(constraints.matrix[i][0]);
        }
    }

} // Конец пространства имен Simplex
