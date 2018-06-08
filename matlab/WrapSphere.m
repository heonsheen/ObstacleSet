classdef WrapSphere
    
    properties
        
        point_P         % bounding-fixed via point 1
        point_S         % bounding-fixed via point 2
        point_O         % Obstacle center point
        radius          % radius of the sphere
        wrap_path_len   % wrap path length
        point_Q         % obstacle via point 1
        point_T         % obstacle via point 2
        mat             % transformation matrix
        point_q         % Q in obstacle frame
        point_t         % T in obstacle frame
        status          % {0: working path, 1: inside radius, 2: no wrap}
    end
    
    methods
        
        function obj = WrapSphere(point_P, point_S, point_O, radius)
            
            obj.point_P = transpose(point_P);
            obj.point_S = transpose(point_S);
            obj.point_O = transpose(point_O);
            obj.point_Q = [0.0, 0.0, 0.0];
            obj.point_T = [0.0, 0.0, 0.0];
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
            vec_N = cross(vec_OP, vec_OS);
            vec_N = vec_N / norm(vec_N);
            
            if (dot(vec_N, [0;0;1]) < 0)
                vec_N = -vec_N;
            end
            
            mat_M = [transpose(vec_OS);
                     transpose(cross(vec_N, vec_OS)); 
                     transpose(vec_N)];
            obj.mat = mat_M;
            disp(mat_M);
                 
            p = mat_M * (obj.point_P - obj.point_O);
            s = mat_M * (obj.point_S - obj.point_O);
            
            p_x = p(1); p_y = p(2);
            s_x = s(1); s_y = s(2);
            R = obj.radius;
            
            denom_q = p_x^2 + p_y^2;
            denom_t = s_x^2 + s_y^2;
            
            if (denom_q - R^2 < 0.0) || (denom_t - R^2 < 0.0)
                %disp(denom_q - R^2);
                %disp(denom_t - R^2);
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
            
            q = [q_x; q_y; 0.0];
            t = [t_x; t_y; 0.0];
            Q = mat_M.' * q + obj.point_O;
            T = mat_M.' * t + obj.point_O;
            Q_x = Q(1); Q_y = Q(2);
            T_x = T(1); T_y = T(2);
            
            obj.point_q = q;
            obj.point_t = t;
            
            obj.point_Q = Q;
            obj.point_T = T;
            
            disp(q);
            disp(t);
            if R * (q_x * t_y - q_y * t_x) > 0.0
                disp(obj.radius * (q_x * t_y - q_y * t_x));
                obj.status = 2;
                disp('nowrap');
                res = obj;
                return;
            end
            
            QT = obj.radius * acos(1.0 - 0.5 * ((Q_x-T_x)^2 + (Q_y-T_y)^2) / obj.radius^2);
            if QT < 0.0
                QT = -QT;
            end
              
            obj.wrap_path_len = QT;
            
            res = obj;
            
            %disp(obj.point_Q);
            %disp(obj.point_T);
            disp(obj.wrap_path_len);
        end
        
        function display_points(obj)
            figure
            x = [obj.point_P(1), obj.point_S(1), obj.point_Q(1), obj.point_T(1), obj.point_O(1)];
            y = [obj.point_P(2), obj.point_S(2), obj.point_Q(2), obj.point_T(2), obj.point_O(2)];
            z = [obj.point_P(3), obj.point_S(3), obj.point_Q(3), obj.point_T(3), obj.point_O(3)];
            
            hold on
            
            text(obj.point_P(1), obj.point_P(2), obj.point_P(3), 'P');
            text(obj.point_S(1), obj.point_S(2), obj.point_S(3), 'S');
            text(obj.point_O(1), obj.point_O(2), obj.point_O(3), 'O');
            plot3(x(:,[1 2 5]), y(:,[1 2 5]), z(:,[1 2 5]), '.', 'MarkerSize', 10);
            
            r = obj.radius;
            [sx, sy, sz] = sphere(50);
            x0 = obj.point_O(1);
            y0 = obj.point_O(2);
            z0 = obj.point_O(3);
            sx = sx*r + x0;
            sy = sy*r + y0;
            sz = sz*r + z0;
            lightGrey = 0.8*[1 1 1];
            surface(sx,sy,sz,'FaceColor', 'none','EdgeColor',lightGrey)
            
            grid on
            axis equal
            
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

                theta_s = min(theta_q, theta_t);
                theta_e = max(theta_q, theta_t);
                if (theta_e - theta_s > theta_s + 2*pi - theta_e)
                    tmp = theta_s;
                    theta_s = theta_e;
                    theta_e = tmp + 2*pi;
                end
                for i = theta_s:(theta_e-theta_s)/100:theta_e
                    P = obj.radius * obj.mat.' * [cos(i), sin(i), 0]' + obj.point_O;
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