#ifndef WrapObst_H
#define WrapObst_H

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

enum Status { wrap, inside_radius, no_wrap, empty };
enum Type { none, sphere, cylinder, double_cylinder }; 
#define PI 3.141593

class WrapObst
{
protected:
    Eigen::Vector3f 
        point_P,    // Bounding-Fixed Via Point 1
        point_S,    // Bounding-Fixed Via Point 2
        point_O,    // Obstacle Center Point
        point_q,    // Obstacle Via Point 1 in Obstacle Frame
        point_t;    // Obstacle Via Point 2 in Obstacle Frame

    Eigen::MatrixXf M;  // Obstacle Coord Transformation Matrix
    Status status;      // Wrapping Status
    Type type;          // Obstacle Type
    float path_length,  // Wrapping Path Length
          radius;       // obstacle sphere radius

public:
    // set muscle origin point
    void setOrigin(const Eigen::Vector3f &P)
    {
        this->point_P = P;
    }

    // set muscle insertion point
    void setInsertion(const Eigen::Vector3f &S)
    {
        this->point_S = S;
    }

    // default constructor
    WrapObst()
    {
        point_P = point_S = point_O = point_q = point_t 
            = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        M = Eigen::MatrixXf(3,3);
        status = empty;
        path_length = 0.0f;
        radius = 0.0f;
        type = none;
    }

    // constructor
    WrapObst(const Eigen::Vector3f &P, 
             const Eigen::Vector3f &S, 
             const Eigen::Vector3f &O, 
             const float R):
        point_P(P), point_S(S), point_O(O), radius(R)
    {
        point_q = point_t = Eigen::Vector3f(0.0f, 0.0f, 0.f);
        M = Eigen::MatrixXf(3,3);
        status = empty;
        path_length = 0.0f;
        type = none;
    }

    // wrap calculation
    void compute() {}
    
    float getLength()
    {
        return this->path_length;
    }

    Status getStatus()
    {
        return this->status;
    }

    float getRadius()
    {
        return this->radius;
    }

    Eigen::MatrixXf getPoints() {}
    
};

#endif
