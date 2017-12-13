#include <iostream>

#include "WrapSphere.hpp"
#include "WrapCylinder.hpp"
#include "WrapDoubleCylinder.hpp"

int main()
{
    /*
    WrapSphere ws(Eigen::Vector3f(2.0f, 0.7f, 0.0f),
                  Eigen::Vector3f(-1.0f, -1.0f, -1.0f),
                  Eigen::Vector3f(0.0f, 0.0f, 0.0f),
                  1);
    ws.compute();
    std::cout << ws.getLength() << std::endl;
    std::cout << ws.getStatus() << std::endl;
    std::cout << ws.getPoints(10) << std::endl;
    */
    /*
    WrapCylinder wc(Eigen::Vector3f(2.0f, 0.7f, 0.0f),
                    Eigen::Vector3f(-1.0f, -1.0f, -1.0f),
                    Eigen::Vector3f(0.0f, 0.0f, 0.0f),
                    Eigen::Vector3f(1.0f, -1.0f, 1.0f),
                    1);
    wc.compute();
    std::cout << wc.getLength() << std::endl;
    std::cout << wc.getStatus() << std::endl;
    std::cout << wc.getPoints(10) << std::endl;
    */
    
    WrapDoubleCylinder wdc(Eigen::Vector3f(-4.0f, 2.0f, -2.0f),
                           Eigen::Vector3f(5.0f, 10.0f, 3.0f),
                           Eigen::Vector3f(0.0f, 0.0f, 0.0f),
                           Eigen::Vector3f(0.0f, 1.0f, 0.0f),
                           2,
                           Eigen::Vector3f(4.0f, 6.0f, 3.0f),
                           Eigen::Vector3f(1.0f, -2.0f, 1.0f),
                           2);
    wdc.compute();
    std::cout << wdc.getLength() << std::endl;
    std::cout << wdc.getStatus() << std::endl;
    std::cout << wdc.getPoints(10) << std::endl;
    
    return 0;
}
