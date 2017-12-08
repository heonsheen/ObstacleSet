#include "WrapDoubleCylinder.hpp"

void WrapDoubleCylinder::compute()
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
    Eigen::Vector3f H0 = H;

    for (int i = 0; i < 100; i++)
    {
        // step 2
        Eigen::Vector3f pu = this->M_U * (this->point_P - this->point_U);
        Eigen::Vector3f hu = this->M_U * (H - this->point_U);

        float denom_q = pu(0)*pu(0) + pu(1)*pu(1);
        float denom_g = hu(0)*hu(0) + hu(1)*hu(1);
        float Ru = this->radius_U;
        
        float root_q = sqrt(denom_q - Ru*Ru);
        float root_g = sqrt(denom_g - Ru*Ru);

        Eigen::Vector3f q(0.0f, 0.0f, 0.0f);
        Eigen::Vector3f g(0.0f, 0.0f, 0.0f);
        q(0) = (pu(0) * Ru*Ru + Ru * pu(1) * root_q) / denom_q;
        q(1) = (pu(1) * Ru*Ru - Ru * pu(0) * root_q) / denom_q;
        g(0) = (hu(0) * Ru*Ru - Ru * hu(1) * root_g) / denom_g;
        g(1) = (hu(1) * Ru*Ru + Ru * hu(0) * root_g) / denom_g;

        float qg_xy = Ru * acos(1.0f - 0.5f *
                                ((q(0) - g(0)) * (q(0) - g(0))
                                 + (q(1) - g(1)) * (q(1) - g(1))) / (Ru*Ru));
        float pq_xy = Ru * acos(1.0f - 0.5f *
                                ((pu(0) - q(0)) * (pu(0) - q(0))
                                 + (pu(1) - q(1)) * (pu(1) - q(1))) / (Ru*Ru));
        float gh_xy = Ru * acos(1.0f - 0.5f *
                                ((g(0) - hu(0)) * (g(0) - hu(0))
                                 + (g(1) - hu(1)) * (g(1) - hu(1))) / (Ru*Ru));
        q(2) = pu(2) + (hu(2)-pu(2)) * pq_xy / (pq_xy + qg_xy + gh_xy);
        g(2) = hu(2) - (hu(2)-pu(2)) * gh_xy / (pq_xy + qg_xy + gh_xy);

        Eigen::Vector3f Q = this->M_U.transpose() * q + this->point_U;
        Eigen::Vector3f G = this->M_U.transpose() * g + this->point_U;
        
        // step 3
        Eigen::Vector3f gv = this->M_V * (this->point_G - this->point_V);

        float denom_h = gv(0)*gv(0) + gv(1)*gv(1);
        float root_h = sqrt(denom_h - Rv*Rv);

        Eigen::Vector3f h(0.0f, 0.0f, 0.0f);
        h(0) = (gv(0) * Rv*Rv + Rv * gv(1) * root_h) / denom_h;
        h(1) = (gv(1) * Rv*Rv - Rv * gv(0) * root_h) / denom_h;
        
        float ht_xy = Rv * acos(1.0f - 0.5f *
                                ((h(0) - t(0)) * (h(0) - t(0))
                                 + (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv));
        float gh_xy = Rv * acos(1.0f - 0.5f *
                                ((gv(0) - h(0)) * (gv(0) - h(0))
                                 + (gv(1) - h(1)) * (gv(1) - h(1))) / (Rv*Rv));
        float ts_xy = Rv * acos(1.0f - 0.5f *
                                ((t(0) - sv(0)) * (t(0) - sv(0))
                                 + (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv));
        h(2) = gv(2) + (sv(2)-gv(2)) * gh_xy / (gh_xy + ht_xy + ts_xy);
        t(2) = sv(2) - (sv(2)-gv(2)) * ts_xy / (gh_xy + ht_xy + ts_xy);

        Eigen::Vector3f H = this->M_V.transpose() * h + this->point_V;
        Eigen::Vector3f T = this->M_V.transpose() * t + this->point_V;

        double dist = (H - H0).norm();
        if (dist == 0) break;
        
        Eigen::Vector3f H0 = H;
    }

    this->q = q;
    this->g = g;
    this->h = h;
    this->t = t;
    std::cout << Q.transpose() << std::endl << G.transpose() << std::endl
              << H.transpose() << std::endl << T.transpose() << std::endl;
}

Eigen::MatrixXf WrapDoubleCylinder::getPoints(int num_points)
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
