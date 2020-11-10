function [M,h,g]=getRobotDyn_planar2DOF(robot_param, q, qdot)
    m1=robot_param.m1; %kg
    m2=robot_param.m2; %kg

    l1=robot_param.l1;
    l2=robot_param.l2; 

  
    g_acc=robot_param.g_acc;
    
    q1=q(1);
    q2=q(2);
    
%     M11= m1*r1^2 + m1*l1*2/12 + m2*(l1^2 + r2^2 +2*l1*r2*cos(q(2))) + m2*l2*2/12;
%     M22= m2*r2^2 + m2*l2*2/12;
%     M12= m2*(r2^2 + l1*r2*cos(q(2))) + m2*l2*2/12;
%     M21=M12;

    M11 = (m1*l1^2 + m2*l1^2 + m2*l2^2) + 2*m2*l1*l2*cos(q2);
    M12 = l2^2*m2 + l1*l2*m2*cos(q2);
    M21 = M12;
    M22 = m2*l2^2;
    
    M=[M11 M12; M21 M22];
    
    h=zeros(2,1);
%     h122 = -m2*l1*r2*sin(q(2));
%     h112 = h122;
%     h211 = -h122;
%     h(1)= h122*qdot(2)^2+2*h112*qdot(1)*qdot(2);
%     h(2)= h211*qdot(1)^2;
    h(1) = -m2*l1*l2*(2*qdot(1)*qdot(2) + qdot(2)^2)*sin(q2);
    h(2) = m2*l1*l2*sin(q2)*qdot(1)^2;
 
    g=zeros(2,1);
%     g(1)= m1*g_acc*r1*cos(q(1)) + m2*g_acc*(l1*cos(q(1))+r2*cos(q(1)+q(2)));
%     g(2)= m2*g_acc*r2*cos(q(1)+q(2));
    g(1) = (m1+m2)*g_acc*l1*cos(q1) + m2*g_acc*l2*cos(q1+q2);
    g(2) = m2*g_acc*l2*cos(q1+q2);

end