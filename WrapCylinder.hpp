#ifndef WrapCylinder_H
#define WrapCylinder_H

/*
 * WrapCylinder.hpp
 * 
 * Obstacle Set Algorithm Simulation for Cylinder Obstacles
 *
 */

#include "WrapObst.hpp"

class WrapCylinder : public WrapObst
{
private:
    Eigen::Vector3f vec_z;      // Cylinder Positive z axis

public:
    // default constructor
    WrapCylinder()
    {
        vec_z = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        type = cylinder;
    }

    void setCylinderConfig(const Eigen::Vector3f &O, 
                           const Eigen::Vector3f &Z)
    {
        this->point_O = O;
        this->vec_z = Z;
    }

    // constructor
    WrapCylinder(const Eigen::Vector3f &P, 
                 const Eigen::Vector3f &S, 
                 const Eigen::Vector3f &O, 
                 const Eigen::Vector3f &Z,
                 const double R)
        : WrapObst(P, S, O, R), vec_z(Z)
    {
        type = cylinder;
    }
    
    using WrapObst::compute;
    void compute();

    using WrapObst::getPoints;
    Eigen::MatrixXf getPoints(int num_points);
};

#endif
