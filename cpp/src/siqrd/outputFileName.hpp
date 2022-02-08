#ifndef OUTPUTFILENAME_HPP
#define OUTPUTFILENAME_HPP
/*
    Function to create file name acording to convetion from 1st assignment.
*/

#include <string>
#include <sstream>

namespace siqrd
{
    template <typename Method, typename sys>
    std::string outputFileName(const sys &system)
    {
        std::stringstream parameters;
        parameters << std::floor(system.beta_ * 100) << "_" << std::floor(system.mu_ * 100) << "_"
                   << std::floor(system.gamma_ * 100) << "_" << std::floor(system.alpha_ * 100) << "_"
                   << std::floor(system.delta_ * 100);
        std::string fileName;
        parameters >> fileName;
        fileName = std::string(Method::method_name) + "_" + fileName + ".out";
        return fileName;
    }
} // namespace siqrd

#endif