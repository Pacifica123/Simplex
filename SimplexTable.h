#include <stdexcept>
#include <string>
#include <vector>
namespace Simplex {
    inline double scalar(std::vector<double> first, std::vector<double> second)
    {
        //TODO: ЭТО ВСЕ ДЕЛО НАДО ПЕРЕМЕСТИТЬ В ФАЙЛ РЕАЛИЗАЦИИ
        if (first.size() != second.size()) {
            // Векторы должны быть одинаковой длины для скалярного произведения.
            throw std::invalid_argument("Векторы должны иметь одинаковую длину для скалярного произведения.");
        }   
        double result = 0.0;
        for (size_t i = 0; i < first.size(); i++) {
            result += first[i] * second[i];
        }
        return result;    
    };
    class BasisCoefficient
    {
        public:
            int literal_numb_A;
            double C_basis;
    };
    class ConstraintSystem
    {
        public:
            std::vector<std::vector<double>> matrix; //[i][0] = b[i]

            ConstraintSystem(std::vector<std::vector<double>> matrix);
            void calculate(int index_for_first, std::vector<double> k);

    };
    class TargetFunc
    {
        public:
            enum Optimization{
                MAXIMUM, MINIMUM
            };
            Optimization target;
            std::vector<double> c;
            std::vector<double> x;

            TargetFunc(std::vector<double> c, Optimization target);
            double sum();

    };
    class SimplexTable
    {
        public:
            std::vector<BasisCoefficient> BeginerBasis;
            std::vector<double> x_plan;
            std::vector<double> z;
            std::vector<double> delta;
            ConstraintSystem constraints;
            TargetFunc Lx;
            const int MAX_ITERATIONS = 100;

            double get_current_Lx();
            BasisCoefficient get_currentA_row(int literal_numb_A); // находим что заменяем
            BasisCoefficient get_newA_colRow(); // через минимум-тета с условими отрицатеьлности А
            void change_beginer_basis(BasisCoefficient currentA, BasisCoefficient newA); // меняем местами
            bool check_delta(); // проверка условия окончания цикла поиска оптимума
            void iterating(); // упорядоченный перебор опорных планов
            void get_beginer(); // нахождение начального базиса и опорного плана
            SimplexTable(ConstraintSystem constraints, TargetFunc Lx);
            std::vector<double> calculate_deltas();
            void print_table(int iteration);
    };
    
}