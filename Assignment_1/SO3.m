k_p = 0.1;

% 목표위치 rpy 0 0 pi/2
eul = [0 0 pi/2];
% eul = [0 0 0];
rotd = eul2rotm(eul);
R_d = rotd

% 목표의 시작 위치
start_point = [0 0 0];
% 목표의 종료 위치
end_point = [0 0 pi/2];
% 목표 점 interpolation
desired_R_euls = interpolator(start_point, end_point, 10);
% 목표점 오일러각을 회전행렬로 변환
desired_R_matrices = eul2rotm(eulpoints);

% 초기위치 rpy 0 0 0
start_point = [0 0 0];
rot0 = eul2rotm(start_point);
R_b = rot0

% omega_d = vee(transpose(R_d)*dot(R_d))이지만, 목표위치의 변화가 없으므로 dot(R_d) = I.
omega_d = vee(transpose(R_d))

R_tilde = transpose(R_b) * R_d
b_R_d = R_tilde;

% 행렬 logarithm은 log 함수가 아닌 expm 함수를 사용해야됨
V = -k_p*vee(expm(R_tilde))
U = b_R_d*(omega_d + V)
% dot_R_b = R_b * hat(U);

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