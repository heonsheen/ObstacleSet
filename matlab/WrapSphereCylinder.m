classdef WrapSphereCylinder
    
    properties
        
        point_P         % bounding-fixed via point 1
        point_S         % bounding-fixed via point 2
        point_O         % Obstacle center point
        vec_z           % vector negative toward cylinder
        radius          % radius of the cylinder
        wrap_path_len   % wrap path length
        point_Q         % obstacle via point 1
        point_T         % obstacle via point 2
        point_W         % obstacle via midpoint
        mat             % transformation matrix
        point_q         % obstacle via point 1 obstacle frame
        point_t         % obstacle via point 2 obstacle frame
        point_ws        % obstacle via midpoint sphere frame
        point_wc        % obstacle via midpoint cylinder frame
        status          % {0: sphere-cylinder, 1: inside radius, 2: no wrap
                        %  3: single sphere, 4: single cylinder }
    end
    
    methods
        
        function obj = WrapSphereCylinder(point_P, point_S, point_O, vec_z, radius)
            
            obj.point_P = point_P.';
            obj.point_S = point_S.';
            obj.point_O = point_O.';
            obj.point_Q = [0.0, 0.0, 0.0];
            obj.point_T = [0.0, 0.0, 0.0];
            obj.point_W = [0.0, 0.0, 0.0];
            obj.vec_z = vec_z.';
            obj.radius = radius;
            obj.wrap_path_len = 0;
            obj.status = 0;
            
            obj = wrap_line(obj);
            display_points(obj);
        end
        
        function res = wrap_line(obj)
             
            vec_OS = obj.point_S - obj.point_O;
            vec_OS = vec_OS / norm(vec_OS);
            vec_OP = obj.point_P - obj.point_O;
            vec_OP = vec_OP / norm(vec_OP);
            
            vec_Z = obj.vec_z / norm(obj.vec_z);
            sum_vec = vec_OS + vec_OP;
            sum_vec = sum_vec / norm(sum_vec);
            if (sum_vec(1) * obj.radius < 0.0)
                sum_vec = -sum_vec;
            end
            
            proj_vec = sum_vec - dot(sum_vec, vec_Z) / norm(vec_Z)^2 * vec_Z;
            
            vec_X = proj_vec / norm(proj_vec);
            vec_Y = cross(vec_Z, vec_X);
            vec_Y = vec_Y / norm(vec_Y);
            
            mat_M = [vec_X.'; vec_Y.'; vec_Z.'];  
            obj.mat = mat_M;
                 
            p = mat_M * (obj.point_P - obj.point_O);
            s = mat_M * (obj.point_S - obj.point_O);
            R = obj.radius;
            
            if (p(3) >= 0.0) && (s(3) >= 0.0)
                % single sphere
                disp('single sphere');
                obj2 = wrap_sphere (obj);
                q = mat_M * (obj2.point_Q - obj2.point_O);
                t = mat_M * (obj2.point_S - obj2.point_O);
                if (q(3) < 0.0) || (t(3) < 0.0)
                    obj2.radius = -obj2.radius;
                    obj2 = wrap_sphere (obj2);
                end
                if obj2.status == 0
                    obj2.status = 3;
                end
                res = obj2;
            elseif (p(3) <= 0.0) && (s(3) <= 0.0)
                % single cylinder
                disp('single cylinder');
                obj2 = wrap_cylinder (obj, p, s, mat_M);
                if obj2.status == 0
                    obj2.status = 4;
                end
                res = obj2;
            else
                if (p(3) < 0.0)
                    t = p; p = s; s = t;
                    t = obj.point_P; obj.point_P = obj.point_S; obj.point_S = t;
%                     obj.radius = -obj.radius;
%                     R = -R;
                end
                
                if (s(1)^2 + s(2)^2 >= R^2)
                    obj2 = wrap_cylinder (obj, p, s, mat_M);
                    q = mat_M * (obj2.point_Q - obj2.point_O);
                    t = mat_M * (obj2.point_S - obj2.point_O);
                    if (obj2.status == 2)...
                        || ( (q(3) <= 0.0) && (t(3) <= 0.0) )
                        disp('case 1');
                        if obj2.status == 0
                            obj2.status = 4;
                        end
                        res = obj2;
                        return;
                    end
                end
                
                obj2 = wrap_sphere (obj);
                q = mat_M * (obj2.point_Q - obj2.point_O);
                t = mat_M * (obj2.point_S - obj2.point_O);
                if (obj2.status == 2)...
                    || ( (q(3) <= 0.0) && (t(3) <= 0.0) )
                    disp('case 2');
                    if obj2.status == 0
                        obj2.status = 3;
                    end
                    res = obj2;
                    return;
                end
                
                p_x = p(1); p_y = p(2); p_z = p(3);
                s_x = s(1); s_y = s(2); s_z = s(3);
                denom_q = p_x^2 + p_y^2;
                if denom_q < R^2
                    obj.status = 1;
                    res = obj;
                    return;
                end
                root_Q = sqrt(denom_q - R^2);
                q_x = (p_x*R^2 + R*p_y*root_Q) / denom_q;
                q_y = (p_y*R^2 - R*p_x*root_Q) / denom_q;
                theta_s = atan(s_y/s_x);
                theta_q = atan(q_y/q_x);
                theta_p = atan(p_y/p_x);
                theta_w = theta_p - (theta_s-theta_p)*p_z/(s_z-p_z);
                disp([theta_s theta_q theta_p theta_w]);
                disp(-p_z*sqrt(s_x^2+s_y^2)*sin(theta_w-theta_s) ...
                    + (s_z*abs(R))*(theta_w-theta_q) ...
                    + (s_z*sqrt((p_x-q_x)^2+(p_y-q_y)^2))*R/abs(R));
                theta_w = fzero(@(theta_w)(...
                    (-p_z*sqrt(s_x^2+s_y^2)*sin(theta_w-theta_s) ...
                    + (s_z*abs(R))*(theta_w-theta_q) ...
                    + (s_z*sqrt((p_x-q_x)^2+(p_y-q_y)^2))*R/abs(R))), ...
                    theta_p - (theta_s-theta_p)*p_z/(s_z-p_z));
                
                if ((R > 0) && (theta_w <= theta_q)) ...
                    || ((R < 0) && (theta_w >= theta_q))
                    disp ('case 3');
                    obj2 = wrap_sphere(obj);
                    q = mat_M * (obj2.point_Q - obj2.point_O);
                    t = mat_M * (obj2.point_S - obj2.point_O);
                    if (q(3) < 0.0) || (t(3) < 0.0)
                        obj2.radius = -obj2.radius;
                        obj2 = wrap_sphere (obj2);
                    end
                    if obj2.status == 0
                        obj2.status = 3;
                    end
                    res = obj2;
                else
                    disp ('case 4');
                    w_x = abs(R) * cos(theta_w);
                    w_y = abs(R) * sin(theta_w);
                    w_z = 0;
                    abs_val = abs(R * theta_w - theta_q);
                    pq_xy = sqrt((p_x-q_x)^2 + (p_y-q_y)^2);
                    q_z = abs_val * p_z / (abs_val + pq_xy);
                    
                    vec_w = [w_x; w_y; w_z];
                    vec_w = vec_w / norm(vec_w);
                    vec_s = s / norm(vec_OP);
                    vec_N = cross(vec_s, vec_w);
                    vec_N = vec_N / norm(vec_N);
                    mat_M2 = [vec_w.'; cross(vec_N, vec_w).'; vec_N.'];
                    vs = mat_M2 * s;
                    s_x = vs(1); s_y = vs(2);
                    denom_t = s_x^2 + s_y^2;
                    root_T = sqrt(denom_t - R^2);
                    t_x = (s_x*R^2 - R*s_y*root_T) / denom_t;
                    t_y = (s_y*R^2 + R*s_x*root_T) / denom_t;
                    vt = [t_x; t_y; 0.0];
                    t = mat_M2.' * vt;
                    
                    obj2.point_Q = mat_M.' * [q_x; q_y; q_z] + obj.point_O;
                    obj2.point_T = mat_M.' * t + obj.point_O;

                    res = obj2;
                end
            end
            
%             disp(obj.point_Q);
%             disp(obj.point_T);
%             disp(obj.wrap_path_len);
        end
        
        function res = wrap_sphere(obj)
            vec_OS = obj.point_S - obj.point_O;
            vec_OS = vec_OS / norm(vec_OS);
            vec_OP = obj.point_P - obj.point_O;
            vec_OP = vec_OP / norm(vec_OP);
            vec_N = cross(vec_OP, vec_OS);
            vec_N = vec_N / norm(vec_N);
            mat_M = [transpose(vec_OS); 
                     transpose(cross(vec_N, vec_OS)); 
                     transpose(vec_N)];
            obj.mat = mat_M;
                 
            p = mat_M * (obj.point_P - obj.point_O);
            s = mat_M * (obj.point_S - obj.point_O);
            
            p_x = p(1); p_y = p(2);
            s_x = s(1); s_y = s(2);
            R = obj.radius;
            
            denom_q = p_x^2 + p_y^2;
            denom_t = s_x^2 + s_y^2;
            if (denom_q - R^2 < 0.0) || (denom_t - R^2 < 0.0)
                disp(denom_q - R^2)
                disp(denom_t - R^2)
                obj.status = 1;
                res = obj;
                return;
            end
            
            root_Q = sqrt(denom_q - R^2);
            root_T = sqrt(denom_t - R^2);
            
            q_x = (p_x*R^2 + R*p_y*root_Q) / denom_q;
            q_y = (p_y*R^2 - R*p_x*root_Q) / denom_q;
            t_x = (s_x*R^2 - R*s_y*root_T) / denom_t;
            t_y = (s_y*R^2 + R*s_x*root_T) / denom_t;
            
            q = [q_x; q_y; 0.0];
            t = [t_x; t_y; 0.0];
            Q = transpose(mat_M) * q + obj.point_O;
            T = transpose(mat_M) * t + obj.point_O;
            Q_x = Q(1); Q_y = Q(2);
            T_x = T(1); T_y = T(2);
            
            obj.point_q = q;
            obj.point_t = t;
            obj.point_Q = Q;
            obj.point_T = T;
            
            if obj.radius * (Q_x * T_y - Q_y * T_x) > 0.0
                disp (Q);
                disp (T);
                disp(obj.radius * (Q_x * T_y - Q_y * T_x));
                obj.status = 2;
                res = obj;
                return;
            end
            
            QT = obj.radius * acos(1.0 - 0.5 * ((Q_x-T_x)^2 + (Q_y-T_y)^2) / obj.radius^2);
            if QT < 0.0
                QT = -QT;
            end
              
            obj.wrap_path_len = QT;
            
            res = obj;
        end
            
        function res = wrap_cylinder(obj, p, s, mat_M)
            p_x = p(1); p_y = p(2); p_z = p(3);
            s_x = s(1); s_y = s(2); s_z = s(3);
            R = obj.radius;
            
            denom_q = p_x^2 + p_y^2;
            denom_t = s_x^2 + s_y^2;
            if (denom_q - R^2 < 0.0) || (denom_t - R^2 < 0.0)
                disp(denom_q - R^2)
                disp(denom_t - R^2)
                disp('insideRadius');
                obj.status = 1;
                res = obj;
                return;
            end
            
            root_Q = sqrt(denom_q - R^2);
            root_T = sqrt(denom_t - R^2);
            
            q_x = (p_x*R^2 + R*p_y*root_Q) / denom_q;
            q_y = (p_y*R^2 - R*p_x*root_Q) / denom_q;
            t_x = (s_x*R^2 - R*s_y*root_T) / denom_t;
            t_y = (s_y*R^2 + R*s_x*root_T) / denom_t;
            
            disp([R q_x q_y t_x t_y]);
            if R * (q_x * t_y - q_y * t_x) > 0.0
                disp(R * (q_x * t_y - q_y * t_x));
                obj.status = 2;
                res = obj;
                return;
            end
            
            QT_xy = R * acos(1.0 - ...
                    0.5 * ((q_x-t_x)^2 + (q_y-t_y)^2) / R^2);
            if QT_xy < 0.0
                QT_xy = -QT_xy;
            end
              
            obj.wrap_path_len = QT_xy;

            PQ_xy = sqrt((p_x-q_x)^2 + (p_y-q_y)^2);
            TS_xy = sqrt((t_x-s_x)^2 + (t_y-s_y)^2);

            q_z = p_z + (s_z-p_z) * PQ_xy / (PQ_xy + QT_xy + TS_xy);
            t_z = s_z - (s_z-p_z) * TS_xy / (PQ_xy + QT_xy + TS_xy);

            pQ = transpose(mat_M) * [q_x; q_y; q_z] + obj.point_O;
            pT = transpose(mat_M) * [t_x; t_y; t_z] + obj.point_O;

            obj.point_q = [q_x; q_y; q_z];
            obj.point_t = [t_x; t_y; t_z];
            obj.point_Q = [pQ(1), pQ(2), pQ(3)];
            obj.point_T = [pT(1), pT(2), pT(3)];
            
            res = obj;
        end
        
        function display_points(obj)            
            
            figure
            x = [obj.point_P(1), obj.point_S(1), obj.point_Q(1), obj.point_T(1), obj.point_W(1), obj.point_O(1)];
            y = [obj.point_P(2), obj.point_S(2), obj.point_Q(2), obj.point_T(2), obj.point_W(2), obj.point_O(2)];
            z = [obj.point_P(3), obj.point_S(3), obj.point_Q(3), obj.point_T(3), obj.point_W(3), obj.point_O(3)];
            
            hold on
            
            r = obj.radius;
            rot = vrrotvec([0,0,1], obj.vec_z);

            [sx, sy, sz] = sphere();
            [cx, cy, cz] = cylinder();
            x0 = obj.point_O(1);
            y0 = obj.point_O(2);
            z0 = obj.point_O(3);
         
            h = mesh(cx*r,cy*r,2*cz*abs(r)-abs(r));
            rotate(h, [rot(1) rot(2) rot(3)], rot(4) / pi * 180);
            vz = obj.vec_z / norm(obj.vec_z);
            h.XData = h.XData + abs(r)*vz(1) + x0;
            h.YData = h.YData + abs(r)*vz(2) + y0;
            h.ZData = h.ZData + abs(r)*vz(3) + z0;
            
            h2 = mesh(sx*r, sy*r, sz*r);
            %rotate(h2, [rot(1) rot(2) rot(3)], rot(4) / pi * 180);
            h2.XData = h2.XData + x0;
            h2.YData = h2.YData + y0;
            h2.ZData = h2.ZData + z0;
            
            grid on
            axis equal
            
            text(obj.point_P(1), obj.point_P(2), obj.point_P(3), 'P');
            text(obj.point_S(1), obj.point_S(2), obj.point_S(3), 'S');
            text(obj.point_O(1), obj.point_O(2), obj.point_O(3), 'O');
            plot3(x(:,[1 2 6]), y(:,[1 2 6]), z(:,[1 2 6]), '.', 'MarkerSize', 10);
            
            if obj.status == 2
                plot3(x(:,[1 2]), y(:,[1 2]), z(:,[1 2]), 'r');
            elseif obj.status == 0
                plot3(x(:,[3 4]), y(:,[3 4]), z(:,[3 4]), '.', 'MarkerSize', 10);

                theta_q = atan(obj.point_q(2)/obj.point_q(1));
                if obj.point_q(1) < 0
                    theta_q = theta_q + pi;
                end
                theta_t = atan(obj.point_t(2)/obj.point_t(1));
                if obj.point_t(1) < 0
                    theta_t = theta_t + pi;
                end

                Ps = [];

                if theta_q > theta_t
                    theta_s = theta_t; theta_e = theta_q;
                    z_s = obj.point_t(3); z_e = obj.point_q(3);
                else
                    theta_e = theta_t; theta_s = theta_q;
                    z_e = obj.point_t(3); z_s = obj.point_q(3);
                end
                
                if (theta_e - theta_s > theta_s + 2*pi - theta_e)
                    tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*pi;
                    tmp = z_s; z_s = z_e; z_e = tmp;
                end
                dz = (z_e - z_s)/100;
                zi = z_s;
                for i = theta_s:(theta_e-theta_s)/100:theta_e
                    P = obj.mat.' * [obj.radius * cos(i), obj.radius * sin(i), zi]'...
                                + obj.point_O;
                    zi = zi + dz;
                    Ps(:,end+1) = P;
                end

                text(obj.point_Q(1), obj.point_Q(2), obj.point_Q(3), 'Q');
                text(obj.point_T(1), obj.point_T(2), obj.point_T(3), 'T');
                plot3(Ps(1,:), Ps(2,:), Ps(3,:), 'r');
                plot3(x(:,[1 3]), y(:,[1 3]), z(:,[1 3]), 'r');
                plot3(x(:,[2 4]), y(:,[2 4]), z(:,[2 4]), 'r');
            end
            
        end
        
    end
end


