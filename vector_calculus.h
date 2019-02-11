#ifndef vector_calculus_GUARD
#define vector_calculus_GUARD

#include <vector>

void add(const std::vector<double>& ,const std::vector<double>& , std::vector<double>& );

void subtract(const std::vector<double>& ,const std::vector<double>& , std::vector<double>& );

void mult(std::vector<double>& , double );

double norm(std::vector<double>& );

double norm2(std::vector<double>& );

void range(std::vector<double>& , double , double );

std::ofstream write(const std::string, const std::vector<std::vector<double> >&);

std::ifstream read(const std::string , std::vector<std::vector<double> >& );

void print(const std::vector<double>);

void print(const std::vector<std::vector<double> > );
#endif // VECTOR_CALCULUS_H_INCLUDED
