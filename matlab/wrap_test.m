for i = 1:100
    p1 = 20 * rand(3, 1);
    p2 = 20 * rand(3, 1);
    p0 = 20 * rand(3, 1);
    R = 20 * rand(1, 1);

    %disp([p1(1), p1(2), p1(3), R(1)]);
    ws(i) = WrapSphere([p1(1), p1(2), p1(3)], ...
                       [p2(1), p2(2), p2(3)], ...
                       [p0(1), p0(2), p0(3)], ...
                       R(1));
end