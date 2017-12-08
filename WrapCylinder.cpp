#include "WrapCylinder.hpp"

void WrapCylinder::compute()
{
    Eigen::Vector3f OP = this->point_P - this->point_O;
    OP = OP / OP.norm();
    Eigen::Vector3f vec_Z = vec_z / vec_z.norm();
    Eigen::Vector3f vec_X = vec_Z.cross(OP);
    vec_X = vec_X / vec_X.norm();
    Eigen::Vector3f vec_Y = vec_Z.cross(vec_X);
    vec_Y = vec_Y / vec_Y.norm();
    
    this->M << vec_X.transpose(), vec_Y.transpose(), vec_Z.transpose();

    Eigen::Vector3f p = this->M * (this->point_P - this->point_O);
    Eigen::Vector3f s = this->M * (this->point_S - this->point_O);

    float denom_q = p(0)*p(0) + p(1)*p(1);
    float denom_t = s(0)*s(0) + s(1)*s(1);
    float R = this->radius;

    if ((denom_q - R*R < 0.0f) || (denom_t - R*R < 0.0f))
    {
        this->status = inside_radius;
    }
        
    float root_q = sqrt(denom_q - R*R);
    float root_t = sqrt(denom_t - R*R);

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

    this->status = wrap;

    float qt_xy = R * acos(1.0f - 0.5f *
                           ((q(0) - t(0)) * (q(0) - t(0))
                            + (q(1) - t(1)) * (q(1) - t(1))) / (R*R));
    this->path_length = qt_xy;

    float pq_xy = sqrt((p(0)-q(0)) * (p(0)-q(0)) + 
                       (p(1)-q(1)) * (p(1)-q(1)));
    float ts_xy = sqrt((t(0)-s(0)) * (t(0)-s(0)) + 
                       (t(1)-s(1)) * (t(1)-s(1)));
    q(2) = p(2) + (s(2)-p(2)) * pq_xy / (pq_xy + qt_xy + ts_xy);
    t(2) = s(2) - (s(2)-p(2)) * ts_xy / (pq_xy + qt_xy + ts_xy);

    this->point_q = q;
    this->point_t = t;
    
    Eigen::Vector3f Q = this->M.transpose() * q + this->point_O;
    Eigen::Vector3f T = this->M.transpose() * t + this->point_O;

    // std::cout << Q.transpose() << std::endl << T.transpose() << std::endl;
}

Eigen::MatrixXf WrapCylinder::getPoints(int num_points)
{
    float theta_q = atan(this->point_q(1) / this->point_q(0));
    if (this->point_q(0) < 0.0f)
        theta_q += PI;
    
    float theta_t = atan(this->point_t(1) / this->point_t(0));
    if (this->point_t(0) < 0.0f)
        theta_t += PI;
    
    Eigen::MatrixXf points(3, num_points + 1);
    
    float theta_s, theta_e, z_s, z_e;
    
    if (theta_q < theta_t)
    {
        theta_s = theta_q; theta_e = theta_t;
        z_s = this->point_q(2); z_e = this->point_t(2);
    }
    else
    {
        theta_s = theta_t; theta_e = theta_q;
        z_s = this->point_t(2); z_e = this->point_q(2);
    }

    if (theta_e - theta_s > theta_s + 2*PI - theta_e)
    {
        float tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*PI;
        tmp = z_s; z_s = z_e; z_e = tmp;
    }

    int col = 0;
    float z_i = z_s, dz = (z_e - z_s) / num_points;
    for (float i = theta_s; i <= theta_e + 0.001; 
         i += (theta_e - theta_s) / num_points)
    {
        Eigen::Vector3f point = this->M.transpose() *
            Eigen::Vector3f(this->radius * cos(i), this->radius * sin(i), z_i) +
            this->point_O;
        z_i += dz;
        points.col(col++) = point;
    }

    return points;
}
