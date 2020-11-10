function [p, J]=forwardKinPlanar2DOF(robot_param, q)

    l1=robot_param.l1;
    l2=robot_param.l2; 
    q1 = q(1);
    q2 = q(2);
    
    p(1) = l1*cos(q1) + l2*cos(q1+q2);
    p(2) = l1*sin(q1) + l2*sin(q1+q2);
    p=[p(1); p(2)];
    
    J11 = -l1*sin(q1) - l2*sin(q1+q2);
    J12 = -l2*sin(q1+q2);
    J21 = l1*cos(q1)+l2*cos(q1+q2);
    J22 = l2*cos(q1+q2);    
    J=[J11 J12; J21 J22];
    
    

end