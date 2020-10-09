k_p = 0.1;

start_R = [0 0 0];
end_R = [0 0 pi/2];

start_r = [0; 0; 0];
end_r = [0; 10; 5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate desired R
desired_R_euls = interpolator(start_R, end_R, 10);
% eul2rotm
desired_R_rotms = eul2rotm(eulpoints);
% interpolate desired r
desired_r_vectors = interpolator(start_r, end_r, 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 시작 r
r_b = [0; 0; 0];
% 목표 r
r_d = [0; 10; 5];
% 시작 R
rot0 = eul2rotm(start_R);
R_b = rot0;
% 목표 R
rotd = eul2rotm(end_R);
R_d = rotd;

dot_r_d = transpose([0 0 0]);
u = transpose(R_b)*(dot_r_d + k_p*(r_d - r_b))

% omega_d = vee(transpose(R_d)*dot(R_d))이지만, 목표위치의 변화가 없으므로 dot(R_d) = I.
omega_d = vee(transpose(R_d));
R_tilde = transpose(R_b) * R_d;
b_R_d = R_tilde;

% 행렬 logarithm은 log 함수가 아닌 expm 함수를 사용해야됨
V = -k_p*vee(expm(R_tilde));
U = b_R_d*(omega_d + V)

% hat 연산자 구현 함수
function ss = hat(vec)
switch length(vec)
    case 3
        ss = [...
            0, -vec(3), vec(2);...
            vec(3), 0, -vec(1);...
            -vec(2), vec(1), 0];       
    case 1
        ss = [...
            0, vec; ...
            -vec, 0];
end
end

% vee 연산자 구현 함수
function vec = vee(ss)
switch numel(ss)
    case 4
        if ~isequal(ss(1,2), -ss(2,1))
            warning('The provided matrix is not skew symmetric')
        end
        vec = ss(1,2);
    case 9
        if ~isequal(ss(3,2),-ss(2,3)) || ~isequal(ss(1,3),-ss(3,1)) || ~isequal(ss(2,1),-ss(1,2))
            warning('The provided matrix is not skew symmetric.')
        end
        vec = [ss(3,2); ss(1,3); ss(2,1)];
end
end

function points = interpolator(start_point, end_point, numbers)
q1 = linspace(start_point(1), end_point(1), numbers);
q2 = linspace(start_point(2), end_point(2), numbers);
q3 = linspace(start_point(3), end_point(3), numbers);
points = [transpose(q1), transpose(q2), transpose(q3)];
end
