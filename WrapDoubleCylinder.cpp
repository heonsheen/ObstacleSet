#include "WrapDoubleCylinder.hpp"

void WrapDoubleCylinder::compute()
{
    // compute Matrix U and V
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
    Eigen::Vector3f vec_X_V = vec_Z_V.cross(OS);
    vec_X_V = vec_X_V / vec_X_V.norm();
    Eigen::Vector3f vec_Y_V = vec_Z_V.cross(vec_X_V);
    vec_Y_V = vec_Y_V / vec_Y_V.norm();

    this->M_U.resize(3,3);
    this->M_U.row(0) = vec_X_U.transpose();
    this->M_U.row(1) = vec_Y_U.transpose();
    this->M_U.row(2) = vec_Z_U.transpose();

    this->M_V.resize(3,3);
    this->M_V.row(0) = vec_X_V.transpose();
    this->M_V.row(1) = vec_Y_V.transpose();
    this->M_V.row(2) = vec_Z_V.transpose();

    // step 1: compute H and T
    Eigen::Vector3f pv = this->M_V * (this->point_P - this->point_V);
    Eigen::Vector3f sv = this->M_V * (this->point_S - this->point_V);

    double denom_h = pv(0)*pv(0) + pv(1)*pv(1);
    double denom_t = sv(0)*sv(0) + sv(1)*sv(1);
    double Rv = this->radius_V;
        
    double root_h = sqrt(denom_h - Rv*Rv);
    double root_t = sqrt(denom_t - Rv*Rv);

    Eigen::Vector3f h(0.0f, 0.0f, 0.0f);
    Eigen::Vector3f t(0.0f, 0.0f, 0.0f);
    h(0) = (pv(0) * Rv*Rv + Rv * pv(1) * root_h) / denom_h;
    h(1) = (pv(1) * Rv*Rv - Rv * pv(0) * root_h) / denom_h;
    t(0) = (sv(0) * Rv*Rv - Rv * sv(1) * root_t) / denom_t;
    t(1) = (sv(1) * Rv*Rv + Rv * sv(0) * root_t) / denom_t;

    if (Rv * (h(0) * t(1) - h(1) * t(0)) > 0.0f)
    {
        status_V = no_wrap;
        h(0) = sv(0);
        h(1) = sv(1);
    }

    this->status = wrap;

    std::complex<double> ht_i = 1.0f - 0.5f *
        ((h(0) - t(0)) * (h(0) - t(0))
         + (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
    std::complex<double> ph_i = 1.0f - 0.5f *
        ((pv(0) - h(0)) * (pv(0) - h(0))
         + (pv(1) - h(1)) * (pv(1) - h(1))) / (Rv*Rv);
    std::complex<double> ts_i = 1.0f - 0.5f *
        ((t(0) - sv(0)) * (t(0) - sv(0))
         + (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

    double ht_xy = abs(Rv * acos(ht_i));
    double ph_xy = abs(Rv * acos(ph_i));
    double ts_xy = abs(Rv * acos(ts_i));

    h(2) = pv(2) + (sv(2)-pv(2)) * ph_xy / (ph_xy + ht_xy + ts_xy);
    t(2) = sv(2) - (sv(2)-pv(2)) * ts_xy / (ph_xy + ht_xy + ts_xy);

    Eigen::Vector3f H = this->M_V.transpose() * h + this->point_V;
    Eigen::Vector3f T = this->M_V.transpose() * t + this->point_V;
    Eigen::Vector3f H0 = H;

    Eigen::Vector3f q(0.0f, 0.0f, 0.0f);
    Eigen::Vector3f g(0.0f, 0.0f, 0.0f);
    Eigen::Vector3f Q, G;

    double len = 0.0f;

    Eigen::Vector3f pu = this->M_U * (this->point_P - this->point_U);
    for (int i = 0; i < 30; i++)
    {
        len = 0.0f;

        // step 2: compute Q and G
        Eigen::Vector3f hu = this->M_U * (H - this->point_U);

        double denom_q = pu(0)*pu(0) + pu(1)*pu(1);
        double denom_g = hu(0)*hu(0) + hu(1)*hu(1);
        double Ru = - this->radius_U;
        
        double root_q = sqrt(denom_q - Ru*Ru);
        double root_g = sqrt(denom_g - Ru*Ru);

        q(0) = (pu(0) * Ru*Ru + Ru * pu(1) * root_q) / denom_q;
        q(1) = (pu(1) * Ru*Ru - Ru * pu(0) * root_q) / denom_q;
        g(0) = (hu(0) * Ru*Ru - Ru * hu(1) * root_g) / denom_g;
        g(1) = (hu(1) * Ru*Ru + Ru * hu(0) * root_g) / denom_g;

        if (Ru * (q(0) * g(1) - q(1) * g(0)) > 0.0f)
        {
            status_U = no_wrap;
            g(0) = pu(0);
            g(1) = pu(1);
        }else {
            status_U = wrap;
        }

        std::complex<double> qg_i = 1.0f - 0.5f *
            ((q(0) - g(0)) * (q(0) - g(0))
             + (q(1) - g(1)) * (q(1) - g(1))) / (Ru*Ru);
        std::complex<double> pq_i = 1.0f - 0.5f *
            ((pu(0) - q(0)) * (pu(0) - q(0))
             + (pu(1) - q(1)) * (pu(1) - q(1))) / (Ru*Ru);
        std::complex<double> gh_i = 1.0f - 0.5f *
            ((g(0) - hu(0)) * (g(0) - hu(0))
             + (g(1) - hu(1)) * (g(1) - hu(1))) / (Ru*Ru);

        double qg_xy = abs(Rv * acos(qg_i));
        double pq_xy = abs(Rv * acos(pq_i));
        double gh_xy = abs(Rv * acos(gh_i));
        len += qg_xy;

        q(2) = pu(2) + (hu(2)-pu(2)) * pq_xy / (pq_xy + qg_xy + gh_xy);
        g(2) = hu(2) - (hu(2)-pu(2)) * gh_xy / (pq_xy + qg_xy + gh_xy);

        Q = this->M_U.transpose() * q + this->point_U;
        G = this->M_U.transpose() * g + this->point_U;

        // step 3: compute H based on G and T
        Eigen::Vector3f gv = this->M_V * (G - this->point_V);

        double denom_h = gv(0)*gv(0) + gv(1)*gv(1);
        double root_h = sqrt(denom_h - Rv*Rv);

        h = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        h(0) = (gv(0) * Rv*Rv + Rv * gv(1) * root_h) / denom_h;
        h(1) = (gv(1) * Rv*Rv - Rv * gv(0) * root_h) / denom_h;

        std::complex<double> ht_i = 1.0f - 0.5f *
            ((h(0) - t(0)) * (h(0) - t(0))
             + (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
        gh_i = 1.0f - 0.5f *
            ((gv(0) - h(0)) * (gv(0) - h(0))
             + (gv(1) - h(1)) * (gv(1) - h(1))) / (Rv*Rv);
        std::complex<double> ts_i = 1.0f - 0.5f *
            ((t(0) - sv(0)) * (t(0) - sv(0))
             + (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

        double ht_xy = abs(Rv * acos(ht_i));
        gh_xy = abs(Rv * acos(gh_i));
        double ts_xy = abs(Rv * acos(ts_i));
        len += ht_xy;

        h(2) = gv(2) + (sv(2)-gv(2)) * gh_xy / (gh_xy + ht_xy + ts_xy);
        t(2) = sv(2) - (sv(2)-gv(2)) * ts_xy / (gh_xy + ht_xy + ts_xy);

        if (Rv * (h(0) * t(1) - h(1) * t(0)) > 0.0f)
        {
            status_V = no_wrap;
            h(0) = sv(0);
            h(1) = sv(1);
        }else {
            status_V = wrap;
        }

        H = this->M_V.transpose() * h + this->point_V;
        T = this->M_V.transpose() * t + this->point_V;
        
        len += (G - H).norm();

        double dist = (H - H0).norm();
        if (dist == 0) break;
        
        H0 = H;
    }

    if (status_V == no_wrap)
    {
        t = sv;
        T = this->M_V.transpose() * t + this->point_V;
    }
    else
    {
        status_V = wrap;
    }
    
    if (status_U == no_wrap)
    {
        q = pu;
        Q = this->M_U.transpose() * q + this->point_U;
    }
    else
    {
        status_U = wrap;
    }

    this->path_length = len;
    this->point_q = q;
    this->point_g = g;
    this->point_h = h;
    this->point_t = t;
    
    /*
    std::cout << Q.transpose() << std::endl << G.transpose() << std::endl
              << H.transpose() << std::endl << T.transpose() << std::endl;
    */
}


Eigen::MatrixXf WrapDoubleCylinder::getPoints(int num_points)
{
    int col = 0;
    
    double theta_s, theta_e, theta_q, theta_g, theta_h, theta_t;
    double z_i, dz, z_s, z_e;

    Eigen::MatrixXf points;

    if (status_U == wrap && status_V == wrap)
        points = Eigen::MatrixXf(3, 3 * num_points + 1);
    else if (status_U == wrap || status_V == wrap)
        points = Eigen::MatrixXf(3, 2 * num_points + 1);
    else
        points = Eigen::MatrixXf(3, 1 * num_points + 1);

    theta_q = atan(this->point_q(1) / this->point_q(0));
    if (this->point_q(0) < 0.0f)
        theta_q += PI;
        
    theta_g = atan(this->point_g(1) / this->point_g(0));
    if (this->point_g(0) < 0.0f)
        theta_g += PI;
    
    theta_h = atan(this->point_h(1) / this->point_h(0));
    if (this->point_h(0) < 0.0f)
        theta_h += PI;
    
    theta_t = atan(this->point_t(1) / this->point_t(0));
    if (this->point_t(0) < 0.0f)
        theta_t += PI;

    // q to g
    if (status_U == wrap)
    {
        if (theta_q < theta_g)
        {
            theta_s = theta_q; theta_e = theta_g;
            z_s = this->point_q(2); z_e = this->point_g(2);
        }
        else
        {
            theta_s = theta_g; theta_e = theta_q;
            z_s = this->point_g(2); z_e = this->point_q(2);
        }

        if (theta_e - theta_s > theta_s + 2*PI - theta_e)
        {
            double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*PI;
            tmp = z_s; z_s = z_e; z_e = tmp;
        }

        z_i = z_s;
        dz = (z_e - z_s) / num_points;
        for (double i = theta_s; i <= theta_e + 0.001; 
             i += (theta_e - theta_s) / num_points)
        {
            if (col == num_points + 1) {
                break;  // During simulation, col sometimes is out of range and error occurs.
            }

            Eigen::Vector3f point = this->M_U.transpose() *
                Eigen::Vector3f(this->radius_U * cos(i), 
                                this->radius_U * sin(i), z_i) +
                this->point_U;
            z_i += dz;
            points.col(col++) = point;
        }
    }
    else
    {
        points.col(col++) = point_P;
    }
    std::cout << col << std::endl;
    
    // g to h
    Eigen::Vector3f G = this->M_U.transpose() * this->point_g + this->point_U;
    Eigen::Vector3f H = this->M_V.transpose() * this->point_h + this->point_V;
    Eigen::Vector3f diff = H - G;

    for (int i = 1; i < num_points; i++)
    {
        Eigen::Vector3f point = G + diff / num_points * i;
        points.col(col++) = point;
    }

    if (status_V == wrap)
    {
        // h to t
        if (theta_h < theta_t)
        {
            theta_s = theta_h; theta_e = theta_t;
            z_s = this->point_h(2); z_e = this->point_t(2);
        }
        else
        {
            theta_s = theta_t; theta_e = theta_h;
            z_s = this->point_t(2); z_e = this->point_h(2);
        }

        if (theta_e - theta_s > theta_s + 2*PI - theta_e)
        {
            double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*PI;
            tmp = z_s; z_s = z_e; z_e = tmp;
        }

        z_i = z_s;
        dz = (z_e - z_s) / num_points;

        int flag = 0;
        // I reverse the order of points, so it is h to t
        for (double i = theta_s; i <= theta_e + 0.001; 
             i += (theta_e - theta_s) / num_points)
        {
            if (flag == num_points + 1) {
                break; // in case it is out of range
            }

            Eigen::Vector3f point = this->M_V.transpose() *
                Eigen::Vector3f(this->radius_V * cos(i), 
                                this->radius_V * sin(i), z_i) +
                this->point_V;
            z_i += dz;
            points.col(col--) = point;
            flag++;
        }
    }
    else
    {
        points.col(col++) = point_V;
    }
    
    return points;
}

