#include "WrapCylinder.hpp"

void WrapCylinder::compute()
{
    Eigen::Vector3f OP = this->point_P - this->point_U;
    OP = OP / OP.norm();
    Eigen::Vector3f vec_Z_U = vec_z_U / vec_z_U.norm();
    Eigen::Vector3f vec_X_U = vec_Z_U.cross(OP);
    vec_X_U = vec_X_U / vec_X_U.norm();
    Eigen::Vector3f vec_Y_U = vec_Z_U.cross(vec_X_U);
    vec_Y_U = vec_Y_U / vec_Y_U.norm();

    Eigen::Vector3f OS = this->point_S - this->point_V;
    OS = OS / OS.norm();
    Eigen::Vector3f vec_Z_V = vec_z_V / vec_z_V.norm();
    Eigen::Vector3f vec_X_V = vec_Z_V.cross(OP);
    vec_X_V = vec_X_V / vec_X_V.norm();
    Eigen::Vector3f vec_Y_V = vec_Z_V.cross(vec_X_V);
    vec_Y_V = vec_Y_V / vec_Y_V.norm();
    
    this->M_U << vec_X_U.transpose(), vec_Y_U.transpose(), vec_Z_U.transpose();
    this->M_V << vec_X_V.transpose(), vec_Y_V.transpose(), vec_Z_V.transpose();

    Eigen::Vector3f pv = this->M_V * (this->point_P - this->point_V);
    Eigen::Vector3f sv = this->M_V * (this->point_S - this->point_V);

    float denom_h = pv(0)*pv(0) + pv(1)*pv(1);
    float denom_t = sv(0)*sv(0) + sv(1)*sv(1);
    float Rv = this->radius_V;
        
    float root_h = sqrt(denom_h - Rv*Rv);
    float root_t = sqrt(denom_t - Rv*Rv);

    Eigen::Vector3f h(0.0f, 0.0f, 0.0f);
    Eigen::Vector3f t(0.0f, 0.0f, 0.0f);
    h(0) = (pv(0) * Rv*Rv + Rv * pv(1) * root_h) / denom_h;
    h(1) = (pv(1) * Rv*Rv - Rv * pv(0) * root_h) / denom_h;
    t(0) = (sv(0) * Rv*Rv - Rv * sv(1) * root_t) / denom_t;
    t(1) = (sv(1) * Rv*Rv + Rv * sv(0) * root_t) / denom_t;

    this->status = wrap;

    float ht_xy = Rv * acos(1.0f - 0.5f *
                            ((h(0) - t(0)) * (h(0) - t(0))
                             + (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv));
    float ph_xy = Rv * acos(1.0f - 0.5f *
                            ((pv(0) - h(0)) * (pv(0) - h(0))
                             + (pv(1) - h(1)) * (pv(1) - h(1))) / (Rv*Rv));
    float ts_xy = Rv * acos(1.0f - 0.5f *
                            ((t(0) - sv(0)) * (t(0) - sv(0))
                             + (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv));
    h(2) = pv(2) + (sv(2)-pv(2)) * ph_xy / (ph_xy + ht_xy + ts_xy);
    t(2) = sv(2) - (sv(2)-pv(2)) * ts_xy / (ph_xy + ht_xy + ts_xy);

    Eigen::Vector3f H = this->M_V.transpose() * h + this->point_V;
    Eigen::Vector3f T = this->M_V.transpose() * t + this->point_V;
    Eigen::Vector3f h0 = h;

    for (int i = 0; i < 100; i++)
    {
        Eigen::Vector3f pu = this->M_U * (this->point_P - this->point_V);
        Eigen::Vector3f su = this->M_U * (this->point_S - this->point_V);

        float denom_h = pv(0)*pv(0) + pv(1)*pv(1);
        float denom_t = sv(0)*sv(0) + sv(1)*sv(1);
        float Rv = this->radius_V;
        
        float root_h = sqrt(denom_h - Rv*Rv);
        float root_t = sqrt(denom_t - Rv*Rv);

        Eigen::Vector3f h(0.0f, 0.0f, 0.0f);
        Eigen::Vector3f t(0.0f, 0.0f, 0.0f);
        h(0) = (pv(0) * Rv*Rv + Rv * pv(1) * root_h) / denom_h;
        h(1) = (pv(1) * Rv*Rv - Rv * pv(0) * root_h) / denom_h;
        t(0) = (sv(0) * Rv*Rv - Rv * sv(1) * root_t) / denom_t;
        t(1) = (sv(1) * Rv*Rv + Rv * sv(0) * root_t) / denom_t;

        this->status = wrap;

        float ht_xy = Rv * acos(1.0f - 0.5f *
                                ((h(0) - t(0)) * (h(0) - t(0))
                                 + (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv));
        float ph_xy = Rv * acos(1.0f - 0.5f *
                                ((pv(0) - h(0)) * (pv(0) - h(0))
                                 + (pv(1) - h(1)) * (pv(1) - h(1))) / (Rv*Rv));
        float ts_xy = Rv * acos(1.0f - 0.5f *
                                ((t(0) - sv(0)) * (t(0) - sv(0))
                                 + (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv));
        h(2) = pv(2) + (sv(2)-pv(2)) * ph_xy / (ph_xy + ht_xy + ts_xy);
        t(2) = sv(2) - (sv(2)-pv(2)) * ts_xy / (ph_xy + ht_xy + ts_xy);

        Eigen::Vector3f H = this->M_V.transpose() * h + this->point_V;
        Eigen::Vector3f T = this->M_V.transpose() * t + this->point_V;
        Eigen::Vector3f h0 = h;
    }

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
