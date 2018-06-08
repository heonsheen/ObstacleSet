#ifndef WrapDoubleCylinder_H
#define WrapDoubleCylinder_H

/*
 * WrapDoubleCylinder.hpp
 * 
 * Obstacle Set Algorithm Simulation for Double Cylinder Obstacles
 *
 */

#include "WrapObst.hpp"

class WrapDoubleCylinder : public WrapObst
{
private:
    Eigen::Vector3f
        vec_z_U,      // U Cylinder Positive z axis
        point_U,      // U Cylinder Origin
        vec_z_V,      // V Cylinder Positive z axis
        point_V,      // V Cylinder Origin
        point_g,
        point_h;

    Eigen::MatrixXf 
        M_U,          // Obstacle Coord Transformation Matrix for U
        M_V;          // Obstacle Coord Transformation Matrix for V
    
    double 
        radius_U,     // U Cylinder Radius
        radius_V;     // V Cylinder Radius
    
    Status 
        status_U,     // U Wrapping Status
        status_V;     // V Wrapping Status

public:
    // default constructor
    WrapDoubleCylinder()
    {
        vec_z_U = point_U = vec_z_V = point_V = point_g = point_h = 
            Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        type = double_cylinder;
    }

    void setCylinderConfig(const Eigen::Vector3f &U, 
                           const Eigen::Vector3f &Z_U,
                           const Eigen::Vector3f &V,
                           const Eigen::Vector3f &Z_V)
    {
        this->point_U = U;
        this->vec_z_U = Z_U;
        this->point_V = V;
        this->vec_z_V = Z_V;
    }

    // constructor
    WrapDoubleCylinder(const Eigen::Vector3f &P, 
                       const Eigen::Vector3f &S,
                       const Eigen::Vector3f &U, 
                       const Eigen::Vector3f &Z_U,
                       const float R_U,
                       const Eigen::Vector3f &V,
                       const Eigen::Vector3f &Z_V,
                       const float R_V)
        : WrapObst(P, S, Eigen::Vector3f(0.0f, 0.0f, 0.0f), 0.0f), 
          point_U(U), vec_z_U(Z_U), radius_U(R_U),
          point_V(V), vec_z_V(Z_V), radius_V(R_V)
    {
        type = double_cylinder;
    }
    
    using WrapObst::compute;
    void compute();

    using WrapObst::getPoints;
    Eigen::MatrixXf getPoints(int num_points);

    Status get_status_u()
    {
        return status_U;
    }

    Status get_status_v()
    {
        return status_V;
    }
};

#endif
