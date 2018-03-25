classdef WrapCylinder
    
    properties
        
        point_P         % bounding-fixed via point 1
        point_S         % bounding-fixed via point 2
        point_O         % Obstacle center point
        vec_z           % z-axis of the cylinder
        radius          % radius of the cylinder
        wrap_path_len   % wrap path length
        point_Q         % obstacle via point 1
        point_T         % obstacle via point 2
        mat             % transformation matrix
        point_q         % Q in obstacle frame
        point_t         % T in obstacle frame
        status          % {0: working path, 1: inside radius, 2: no wrap}
    end
    
    methods
        
        function obj = WrapCylinder(point_P, point_S, point_O, vec_z, radius)
            
            obj.point_P = point_P.';
            obj.point_S = point_S.';
            obj.point_O = point_O.';
            obj.point_Q = [0.0, 0.0, 0.0];
            obj.point_T = [0.0, 0.0, 0.0];
            obj.point_q = [0.0, 0.0, 0.0];
            obj.point_t = [0.0, 0.0, 0.0];
            obj.radius = radius;
            obj.vec_z = vec_z.';
            obj.wrap_path_len = 0;
            obj.status = 0;
            obj = wrap_line(obj);
%             display_points(obj);
        end
        
        function res = wrap_line(obj)
             
            vec_OP = obj.point_P - obj.point_O;
            vec_OP = vec_OP / norm(vec_OP);
            vec_Z = obj.vec_z / norm(obj.vec_z);
            vec_X = cross(vec_Z, vec_OP);
            vec_X = vec_X / norm(vec_X);
            vec_Y = cross(vec_Z, vec_X);
            vec_Y = vec_Y / norm(vec_Y);
            
            mat_M = [vec_X.'; vec_Y.'; vec_Z.'];
            obj.mat = mat_M;
                 
            p = mat_M * (obj.point_P - obj.point_O);
            s = mat_M * (obj.point_S - obj.point_O);
            
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
                          
            if R * (q_x * t_y - q_y * t_x) > 0.0
                %disp(R * (q_x * t_y - q_y * t_x));
                obj.status = 2;
                %res = obj;
                %return;
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

            obj.point_q = [q_x, q_y, q_z]';
            obj.point_t = [t_x, t_y, t_z]';

            pQ = mat_M.' * [q_x; q_y; q_z] + obj.point_O;
            pT = mat_M.' * [t_x; t_y; t_z] + obj.point_O;

            obj.point_Q = [pQ(1), pQ(2), pQ(3)];
            obj.point_T = [pT(1), pT(2), pT(3)];
            
            res = obj;
            
%             disp(obj.point_Q);
%             disp(obj.point_T);
%             disp(obj.wrap_path_len);
        end
        
        function display_points(obj)
            figure
            x = [obj.point_P(1), obj.point_S(1), obj.point_Q(1), obj.point_T(1), obj.point_O(1)];
            y = [obj.point_P(2), obj.point_S(2), obj.point_Q(2), obj.point_T(2), obj.point_O(2)];
            z = [obj.point_P(3), obj.point_S(3), obj.point_Q(3), obj.point_T(3), obj.point_O(3)];
            
            hold on
            
            r = obj.radius;
            rot = vrrotvec([0,0,1], obj.vec_z);

            [cx, cy, cz] = cylinder();
            x0 = obj.point_O(1);
            y0 = obj.point_O(2);
            z0 = obj.point_O(3);

            h = mesh(cx*r,cy*r,4*cz*abs(r)-2*abs(r));
            rotate(h, [rot(1) rot(2) rot(3)], rot(4) / pi * 180);
            %vz = obj.vec_z / norm(obj.vec_z);
            vz = [0, 0, 0];
            h.XData = h.XData + abs(r)*vz(1) + x0;
            h.YData = h.YData + abs(r)*vz(2) + y0;
            h.ZData = h.ZData + abs(r)*vz(3) + z0;
            
            grid on
            axis equal
            
            text(obj.point_P(1), obj.point_P(2), obj.point_P(3), 'P');
            text(obj.point_S(1), obj.point_S(2), obj.point_S(3), 'S');
            text(obj.point_O(1), obj.point_O(2), obj.point_O(3), 'O');
            plot3(x(:,[1 2 5]), y(:,[1 2 5]), z(:,[1 2 5]), '.', 'MarkerSize', 10);
            
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
