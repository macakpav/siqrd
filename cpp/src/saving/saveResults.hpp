#ifndef SAVERESULTS_HPP
#define SAVERESULTS_HPP
/*
    funciton to save a matrix to file with first column containing increments of dT
*/

#include <fstream>
namespace saving
{
    template <typename T, typename matrix_type, typename string>
    void saveResults(const T dT, const matrix_type &variables_matrix, const string &file_name)
    {

#ifndef NINFO
        std::cout << "Writing results to " << file_name << "...       ";
#endif
        T time = 0.0;
        std::ofstream outputFile(file_name);
        for (decltype(variables_matrix.size2()) j = 0; j < variables_matrix.size2(); j++)
        {
            outputFile << time << "  \t";
            for (decltype(variables_matrix.size1()) i = 0; i < variables_matrix.size1(); i++)
            {
                outputFile << variables_matrix(i, j) << "  \t";
            }
            outputFile << std::endl;
            time += dT;
        }
#ifndef NINFO
        std::cout << "done." << std::endl
                  << std::endl;
#endif
    }
} // namespace saving

#endif