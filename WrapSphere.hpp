#ifndef WrapSphere_H
#define WrapSphere_H

/*
 * WrapSphere.hpp
 * 
 * Obstacle Set Algorithm Simulation for Sphere Obstacles
 *
 */

#include "WrapObst.hpp"

class WrapSphere : public WrapObst
{   
public:
    // default constructor
    WrapSphere()
    {
        type = sphere;
    }

    void setSphereConfig(const Eigen::Vector3f &O)
    {
        this->point_O = O;
    }

    // constructor
    WrapSphere(const Eigen::Vector3f &P, 
                const Eigen::Vector3f &S, 
                const Eigen::Vector3f &O, 
               const float R)
        : WrapObst(P, S, O, R)
    {
        type = sphere;
    }

    using WrapObst::compute;
    void compute();
    
    using WrapObst::getPoints;
    Eigen::MatrixXf getPoints(int num_points);
};

#endif
