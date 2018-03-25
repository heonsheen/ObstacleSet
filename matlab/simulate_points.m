function simulate_points(type, point_P, point_S)

    n_ps = [1 3 3];      % number of simulated point S
    
    nmax = n_ps(1)*n_ps(2)*n_ps(3);
    Ss = zeros(3,nmax);
    Qs = zeros(3,nmax);
    Ts = zeros(3,nmax);
    qs = zeros(3,nmax);
    ts = zeros(3,nmax);
    Ms = zeros(3,3,nmax);
    status = zeros(nmax);
    
    figure
    hold on
    
    if type == 0
        % sphere
        point_O = [0 0 1];
        R = 1;
        for i = 1: n_ps(1)
            for j = 1: n_ps(2)
                for k = 1: n_ps(3)
                    n = (i-1)*n_ps(2)*n_ps(3) + (j-1)*n_ps(3) + k;
                    dP = 0.3*[ceil(n_ps(1))/2-i ceil(n_ps(2))/2-j ceil(n_ps(3))/2-k];
                    %disp(dP);
                    Ss(:,n) = (point_S + dP).';
                    %disp(Ss(:,n).');
                    obj = WrapSphere(point_P, Ss(:,n).', point_O, R);
                    Ms(:,:,n) = obj.mat;
                    if obj.status == 0
                        Qs(:,n) = obj.point_Q;
                        Ts(:,n) = obj.point_T;
                        qs(:,n) = obj.point_q;
                        ts(:,n) = obj.point_t;
                    end
                    status(n) = obj.status;
                end
            end
        end
        
        text(point_O(1), point_O(2), point_O(3), 'O');
        text(point_P(1), point_P(2), point_P(3), 'P');
        plot3([point_O(1) point_P(1)], [point_O(2) point_P(2)], ...
                [point_O(3) point_P(3)], '.', 'MarkerSize', 10);
        
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

        for i = 1: nmax
            text(Ss(1,i), Ss(2,i), Ss(3,i), strcat('S',int2str(i)));
            plot3([Ss(1,i) Qs(1,i) Ts(1,i)], ...
                  [Ss(2,i) Qs(2,i) Ts(2,i)], ...
                  [Ss(3,i) Qs(3,i) Ts(3,i)], '.', 'MarkerSize', 10);
              
            if status(i) == 2
                plot3([point_P(1) Ss(1,i)], [point_P(2) Ts(S,i)], ...
                      [point_P(3) Ss(3,i)], 'r');
            elseif status(i) == 0

                theta_q = atan(qs(2,i)/qs(1,i));
                if qs(1,i) < 0
                    theta_q = theta_q + pi;
                end
                theta_t = atan(ts(2,i)/ts(1,i));
                if ts(1,i) < 0
                    theta_t = theta_t + pi;
                end

                theta_s = min(theta_q, theta_t);
                theta_e = max(theta_q, theta_t);
                if (theta_e - theta_s > theta_s + 2*pi - theta_e)
                    tmp = theta_s;
                    theta_s = theta_e;
                    theta_e = tmp + 2*pi;
                end
                pos = [];
                for j = theta_s:(theta_e-theta_s)/99:theta_e
                    po = R * Ms(:,:,i).' * [cos(j), sin(j), 0]' + point_O.';
                    pos(:,end+1) = po;
                end

                text(Qs(1,i), Qs(2,i), Qs(3,i), strcat('Q',int2str(i)));
                text(Ts(1,i), Ts(2,i), Ts(3,i), strcat('T',int2str(i)));
                plot3(pos(1,:), pos(2,:), pos(3,:), 'r');
                plot3([point_P(1) Qs(1,i)], [point_P(2), Qs(2,i)], [point_P(3), Qs(3,i)], 'r');
                plot3([Ts(1,i) Ss(1,i)], [Ts(2,i), Ss(2,i)], [Ts(3,i), Ss(3,i)], 'r');
            end
        end
        
    elseif type == 1
        % cylinder
        point_O = [0 0 1];
        vec_z = [1 -1 1];
        R = 1;
        for i = 1: n_ps(1)
            for j = 1: n_ps(2)
                for k = 1: n_ps(3)
                    n = (i-1)*n_ps(2)*n_ps(3) + (j-1)*n_ps(3) + k;
                    dP = [ceil(n_ps(1))/2-i ceil(n_ps(2))/2-j ceil(n_ps(3))/2-k];
                    Ss(:,n) = (point_S + dP).';
                    obj = WrapCylinder(point_P, Ss(:,n).', point_O, vec_z, R);
                    Qs(:,n) = obj.point_Q;
                    Ts(:,n) = obj.point_T;
                    qs(:,n) = obj.point_q;
                    ts(:,n) = obj.point_t;
                    Ms(:,:,n) = obj.mat;
                    status(n) = obj.status;
                end
            end
        end
        
        text(point_O(1), point_O(2), point_O(3), 'O');
        text(point_P(1), point_P(2), point_P(3), 'P');
        plot3([point_O(1) point_P(1)], [point_O(2) point_P(2)], ...
                [point_O(3) point_P(3)], '.', 'MarkerSize', 10);
        
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

        for i = 1: nmax
            text(Ss(1,i), Ss(2,i), Ss(3,i), strcat('S',int2str(i)));
            plot3([Ss(1,i) Qs(1,i) Ts(1,i)], ...
                  [Ss(2,i) Qs(2,i) Ts(2,i)], ...
                  [Ss(3,i) Qs(3,i) Ts(3,i)], '.', 'MarkerSize', 10);
              
            if status(i) == 2
                plot3([point_P(1) Ss(1,i)], [point_P(2) Ss(2,i)], ...
                      [point_P(3) Ss(3,i)], 'r');
            elseif status(i) == 0
                theta_q = atan(qs(2,i)/qs(1,i));
                if qs(1,i) < 0
                    theta_q = theta_q + pi;
                end
                theta_t = atan(ts(2,i)/ts(1,i));
                if ts(1,i) < 0
                    theta_t = theta_t + pi;
                end

                pos = [];

                if theta_q > theta_t
                    theta_s = theta_t; theta_e = theta_q;
                    z_s = ts(3,i); z_e = qs(3,i);
                else
                    theta_e = theta_t; theta_s = theta_q;
                    z_e = ts(3,i); z_s = qs(3,i);
                end
                
                if (theta_e - theta_s > theta_s + 2*pi - theta_e)
                    tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*pi;
                    tmp = z_s; z_s = z_e; z_e = tmp;
                end
                dz = (z_e - z_s)/100;
                zj = z_s;
                disp([theta_s theta_e]);
                for j = theta_s:(theta_e-theta_s)/100:theta_e
                    po = Ms(:,:,i).' * [R * cos(j), R * sin(j), zj]' + point_O.';
                    zj = zj + dz;
                    pos(:,end+1) = po;
                end

                text(Qs(1,i), Qs(2,i), Qs(3,i), 'Q');
                text(Ts(1,i), Ts(2,i), Ts(3,i), 'T');
                plot3(pos(1,:), pos(2,:), pos(3,:), 'r');
                plot3([point_P(1) Qs(1,i)], [point_P(2), Qs(2,i)], [point_P(3), Qs(3,i)], 'r');
                plot3([Ts(1,i) Ss(1,i)], [Ts(2,i), Ss(2,i)], [Ts(3,i), Ss(3,i)], 'r');
            end
        end
    elseif type == 2
        % sphere-cylinder
        Ts = zeros(n_ps(1)*n_ps(2)*n_ps(3),3);
    elseif type == 3
        % double cylinder
        Gs = zeros(n_ps(1)*n_ps(2)*n_ps(3),3);
        Hs = zeros(n_ps(1)*n_ps(2)*n_ps(3),3);
        
        point_O = [0 0 1];
        vec_z = [1 -1 1];
        R = 1;
        for i = 1: n_ps(1)
            for j = 1: n_ps(2)
                for k = 1: n_ps(3)
                    n = (i-1)*n_ps(2)*n_ps(3) + (j-1)*n_ps(3) + k;
                    dP = [ceil(n_ps(1))/2-i ceil(n_ps(2))/2-j ceil(n_ps(3))/2-k];
                    Ss(:,n) = (point_S + dP).';
                    obj = WrapCylinder(point_P, Ss(:,n).', point_O, vec_z, R);
                    Qs(:,n) = obj.point_Q;
                    Ts(:,n) = obj.point_T;
                    qs(:,n) = obj.point_q;
                    ts(:,n) = obj.point_t;
                    Ms(:,:,n) = obj.mat;
                    status(n) = obj.status;
                end
            end
        end
        
        text(point_O(1), point_O(2), point_O(3), 'O');
        text(point_P(1), point_P(2), point_P(3), 'P');
        plot3([point_O(1) point_P(1)], [point_O(2) point_P(2)], ...
                [point_O(3) point_P(3)], '.', 'MarkerSize', 10);
        
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

        for i = 1: nmax
            text(Ss(1,i), Ss(2,i), Ss(3,i), strcat('S',int2str(i)));
            plot3([Ss(1,i) Qs(1,i) Ts(1,i)], ...
                  [Ss(2,i) Qs(2,i) Ts(2,i)], ...
                  [Ss(3,i) Qs(3,i) Ts(3,i)], '.', 'MarkerSize', 10);
              
            if status(i) == 2
                plot3([point_P(1) Ss(1,i)], [point_P(2) Ss(2,i)], ...
                      [point_P(3) Ss(3,i)], 'r');
            elseif status(i) == 0
                theta_q = atan(qs(2,i)/qs(1,i));
                if qs(1,i) < 0
                    theta_q = theta_q + pi;
                end
                theta_t = atan(ts(2,i)/ts(1,i));
                if ts(1,i) < 0
                    theta_t = theta_t + pi;
                end

                pos = [];

                if theta_q > theta_t
                    theta_s = theta_t; theta_e = theta_q;
                    z_s = ts(3,i); z_e = qs(3,i);
                else
                    theta_e = theta_t; theta_s = theta_q;
                    z_e = ts(3,i); z_s = qs(3,i);
                end
                
                if (theta_e - theta_s > theta_s + 2*pi - theta_e)
                    tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*pi;
                    tmp = z_s; z_s = z_e; z_e = tmp;
                end
                dz = (z_e - z_s)/100;
                zj = z_s;
                disp([theta_s theta_e]);
                for j = theta_s:(theta_e-theta_s)/100:theta_e
                    po = Ms(:,:,i).' * [R * cos(j), R * sin(j), zj]' + point_O.';
                    zj = zj + dz;
                    pos(:,end+1) = po;
                end

                text(Qs(1,i), Qs(2,i), Qs(3,i), 'Q');
                text(Ts(1,i), Ts(2,i), Ts(3,i), 'T');
                plot3(pos(1,:), pos(2,:), pos(3,:), 'r');
                plot3([point_P(1) Qs(1,i)], [point_P(2), Qs(2,i)], [point_P(3), Qs(3,i)], 'r');
                plot3([Ts(1,i) Ss(1,i)], [Ts(2,i), Ss(2,i)], [Ts(3,i), Ss(3,i)], 'r');
            end
        end
    elseif type == 4    % muscle simulation
        % cylinder
        point_O = [0 0 1];
        vec_z = [-1 -1 1];
        vec_dir = [1 -1 1];
        R = 1;
        t0 = 0.05;
        n = 1;
        prev_len = 0;
        
        for i = 0: t0: 4
%             n = (i-1)*n_ps(2)*n_ps(3) + (j-1)*n_ps(3) + k;
            cla;

            rot = vrrotvec([0,0,1], vec_z);

            [cx, cy, cz] = cylinder();
            x0 = point_O(1);
            y0 = point_O(2);
            z0 = point_O(3);

            h = mesh(cx*R,cy*R,4*cz*abs(R)-2*abs(R));
            rotate(h, [rot(1) rot(2) rot(3)], rot(4) / pi * 180);
            %vz = obj.vec_z / norm(obj.vec_z);
            vz = [0, 0, 0];
            h.XData = h.XData + abs(R)*vz(1) + x0;
            h.YData = h.YData + abs(R)*vz(2) + y0;
            h.ZData = h.ZData + abs(R)*vz(3) + z0;

            grid on
            axis equal

            text(point_O(1), point_O(2), point_O(3), 'O');
            text(point_P(1), point_P(2), point_P(3), 'P');
            plot3([point_O(1) point_P(1)], [point_O(2) point_P(2)], ...
                    [point_O(3) point_P(3)], '.', 'MarkerSize', 10);
            
            dP = i * vec_dir;
            Ss(:,n) = (point_S + dP).';
            obj = WrapCylinder(point_P, Ss(:,n).', point_O, vec_z, R);
            Qs(:,n) = obj.point_Q;
            Ts(:,n) = obj.point_T;
            qs(:,n) = obj.point_q;
            ts(:,n) = obj.point_t;
            Ms(:,:,n) = obj.mat;
            status(n) = obj.status;
            
            text(Ss(1,n), Ss(2,n), Ss(3,n), 'S');

            if i == 0
                prev_len = obj.wrap_path_len;
            end
            
            if (status(n) == 2) && (...
                    (abs(obj.wrap_path_len - prev_len) >= abs(prev_len)) || ...
                    (i == 0))
                plot3([point_P(1) Ss(1,n)], [point_P(2) Ss(2,n)], ...
                      [point_P(3) Ss(3,n)], 'r');
            else
                plot3([Ss(1,n) Qs(1,n) Ts(1,n)], ...
                  [Ss(2,n) Qs(2,n) Ts(2,n)], ...
                  [Ss(3,n) Qs(3,n) Ts(3,n)], '.', 'MarkerSize', 10);
                 
                theta_q = atan(qs(2,n)/qs(1,n));
                if qs(1,n) < 0
                    theta_q = theta_q + pi;
                end
                theta_t = atan(ts(2,n)/ts(1,n));
                if ts(1,n) < 0
                    theta_t = theta_t + pi;
                end

                pos = [];

                if theta_q > theta_t
                    theta_s = theta_t; theta_e = theta_q;
                    z_s = ts(3,n); z_e = qs(3,n);
                else
                    theta_e = theta_t; theta_s = theta_q;
                    z_e = ts(3,n); z_s = qs(3,n);
                end

                if (theta_e - theta_s > theta_s + 2*pi - theta_e)
                    tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2*pi;
                    tmp = z_s; z_s = z_e; z_e = tmp;
                end
                dz = (z_e - z_s)/100;
                zj = z_s;
                %disp([theta_s theta_e]);
                for j = theta_s:(theta_e-theta_s)/100:theta_e
                    po = Ms(:,:,n).' * [R * cos(j), R * sin(j), zj]' + point_O.';
                    zj = zj + dz;
                    pos(:,end+1) = po;
                end

                text(Qs(1,n), Qs(2,n), Qs(3,n), 'Q');
                text(Ts(1,n), Ts(2,n), Ts(3,n), 'T');
                plot3(pos(1,:), pos(2,:), pos(3,:), 'r');
                plot3([point_P(1) Qs(1,n)], [point_P(2), Qs(2,n)], [point_P(3), Qs(3,n)], 'r');
                plot3([Ts(1,n) Ss(1,n)], [Ts(2,n), Ss(2,n)], [Ts(3,n), Ss(3,n)], 'r');
                prev_len = obj.wrap_path_len;
            end
            drawnow;
            
        end
    end
end