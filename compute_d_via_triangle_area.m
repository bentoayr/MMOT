% compute d using counter example 

function area = compute_d_via_triangle_area(x,gamma)

    Ax = x(1,1);
    Ay = x(2,1);
    Bx = x(1,2);
    By = x(2,2);
    Cx = x(1,3);
    Cy = x(2,3);

    area = abs(Ax*(By - Cy) + Bx*(Cy - Ay) + Cx*(Ay - By))/2;
    d12 = norm(x(:,1) - x(:,2));
    d13 = norm(x(:,1) - x(:,3));
    d23 = norm(x(:,2) - x(:,3));
    
    if (     ( (d12 < 10*eps)      &&    (d13 > 10*eps) )     ||    ( (d13 < 10*eps) && (d12 > 10*eps) )         ||          (     (d23 < 10*eps)  && (d13 > 10*eps)      ) )
        area = gamma;
    end
    
end