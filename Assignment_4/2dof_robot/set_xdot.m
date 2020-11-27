function [ xdot ] = set_xdot( M,h,g,q_dot, tau )

    q_ddot = inv(M)* (tau - h - g);

    xdot(1)=q_dot(1);
    xdot(2)=q_dot(2);
    xdot(3)=q_ddot(1);
    xdot(4)=q_ddot(2);
    
    % just to make sure that x_dot is a column vector (not a row)
    xdot = [xdot(1);xdot(2);xdot(3);xdot(4)]; 

end