#include <iostream>

#include "WrapSphere.hpp"
#include "WrapCylinder.hpp"
#include "WrapDoubleCylinder.hpp"

int main()
{
    /*
    WrapSphere ws(Eigen::Vector3f(0.0f, 3.0f, 0.0f),
                  Eigen::Vector3f(5.0f, 0.75f, 0.0f),
                  Eigen::Vector3f(4.0f, 1.0f, 0.0f),
                  0.5f);
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
    
    WrapDoubleCylinder wdc(Eigen::Vector3f(2.81523f, 0.579136f, 0.0f),
                           Eigen::Vector3f(8.44334f, -3.15549f, 0.0f),
                           Eigen::Vector3f(4.96025f, -1.57699, 0.0f),
                           Eigen::Vector3f(0.0f, 0.0f, -1.0f),
                           0.1,
                           Eigen::Vector3f(7.57674f, -2.44909f, 0.0f),
                           Eigen::Vector3f(0.0f, 0.0f, 1.0f),
                           0.1);
    wdc.compute();
    std::cout << wdc.getLength() << std::endl;
    //std::cout << wdc.getStatus() << std::endl;   
    std::cout << wdc.get_status_u() << std::endl;
    std::cout << wdc.get_status_v() << std::endl;

    std::cout << wdc.getPoints(10) << std::endl;
    
    return 0;
}
