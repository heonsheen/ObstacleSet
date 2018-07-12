#include "WrapSphere.hpp"

void WrapSphere::compute()
{
    Eigen::Vector3f OS = this->point_S - this->point_O;
    OS = OS / OS.norm();
    Eigen::Vector3f OP = this->point_P - this->point_O;
    OP = OP / OP.norm();
    Eigen::Vector3f N = OP.cross(OS);
    N = N / N.norm();

    if (N.dot(Eigen::Vector3f(0.0f, 0.0f, 1.0f)) < 0)
        N = -N;
    
    this->M << OS.transpose(), N.cross(OS).transpose(), N.transpose();
    //std::cout << this->M << std::endl;

    Eigen::Vector3f p = this->M * (this->point_P - this->point_O);
    Eigen::Vector3f s = this->M * (this->point_S - this->point_O);

    double denom_q = p(0)*p(0) + p(1)*p(1);
    double denom_t = s(0)*s(0) + s(1)*s(1);
    double R = this->radius;

    this->status = wrap;

    if ((denom_q - R*R < 0.0f) || (denom_t - R*R < 0.0f))
    {
        this->status = inside_radius;
        
    }
        
    double root_q = sqrt(denom_q - R*R);
    double root_t = sqrt(denom_t - R*R);

    Eigen::Vector3f q(0.0f, 0.0f, 0.0f);
    Eigen::Vector3f t(0.0f, 0.0f, 0.0f);
    q(0) = (p(0) * R*R + R * p(1) * root_q) / denom_q;
    q(1) = (p(1) * R*R - R * p(0) * root_q) / denom_q;
    t(0) = (s(0) * R*R - R * s(1) * root_t) / denom_t;
    t(1) = (s(1) * R*R + R * s(0) * root_t) / denom_t;

    if (R * (q(0) * t(1) - q(1) * t(0)) > 0.0f)
    {
        this->status = no_wrap;
    }

    
    this->point_q = q;
    this->point_t = t;

    //std::cout << q << std::endl << t << std::endl;
    
    Eigen::Vector3f Q = this->M.transpose() * q + this->point_O;
    Eigen::Vector3f T = this->M.transpose() * t + this->point_O;

    //  std::cout << Q.transpose() << std::endl << T.transpose() << std::endl;

    this->path_length = R * acos(1.0f - 0.5f *
                                 ((Q(0) - T(0)) * (Q(0) - T(0))
                                  + (Q(1) - T(1)) * (Q(1) - T(1))) / (R*R));
}

Eigen::MatrixXf WrapSphere::getPoints(int num_points)
{
    double theta_q = atan(this->point_q(1) / this->point_q(0));
    if (this->point_q(0) < 0.0f)
        theta_q += PI;
    
    double theta_t = atan(this->point_t(1) / this->point_t(0));
    if (this->point_t(0) < 0.0f)
        theta_t += PI;
    
    Eigen::MatrixXf points(3, num_points + 1);
    
    double theta_s, theta_e;
    
    if (theta_q < theta_t)
    {
        theta_s = theta_q; theta_e = theta_t;
    }
    else
    {
        theta_s = theta_t; theta_e = theta_q;
    }

    if (theta_e - theta_s > theta_s + 2*PI - theta_e)
    {
        double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*PI;
    }

    int col = 0;    
    for (double i = theta_s; i <= theta_e + 0.001; 
         i += (theta_e - theta_s) / num_points)
    {
        if (col == num_points + 1) {
            break;
        }
        
        Eigen::Vector3f point = this->radius * this->M.transpose() *
            Eigen::Vector3f(cos(i), sin(i), 0.0f) + this->point_O;
        points.col(col++) = point;
    }

    return points;
}
