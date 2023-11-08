#include <iostream>
#include <ostream>
#include <vector>
#include "SimplexTable.h"

int main(){
    // std::cout << "Hello world!" << std::endl;
    // std::vector<double> c = {11.0, 2.5, 2.1} // - коэффициенты при целевой
    Simplex::TargetFunc Lx = Simplex::TargetFunc({11.0, 2.5, 2.1}, Simplex::TargetFunc::MAXIMUM);

    // Ограничения
    Simplex::ConstraintSystem cs = Simplex::ConstraintSystem({
        {0.25, 0.05, 0.025},
        {2, 0.5, 0.4}
    });

    Simplex::SimplexTable table = Simplex::SimplexTable(cs, Lx);

    table.iterating();

    // std::cout << table.Lx.sum() << std::endl;
    return 0;
}