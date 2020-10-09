clear all;

k_v = 10.0;
k_omega = 5.0;
k_d = 0.1;

dt = 0.01;
interpol_size = 500;
time_vec = 0.0:0.01:0.01*(interpol_size-1);

start_R = [0 0 0];
end_R = [pi/2 pi/2 pi/4];

start_r = [0; 0; 0];
end_r = [1; 1; 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate desired R
desired_R_euls = interpolator(start_R, end_R, interpol_size);
% eul2rotm
desired_R_rotms = eul2rotm(desired_R_euls);
% interpolate desired r
desired_r_vectors = transpose(interpolator(start_r, end_r, interpol_size));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 시작 r
r_b = start_r;

% 시작 R
R_b = eul2rotm(start_R);

r_b_list = zeros(3, interpol_size);
R_b_list = zeros(3, 3, interpol_size);
for i = 1:interpol_size
    if i == 1
        dot_r_d = [0; 0; 0];
        dot_R_d = [0 0 0; 0 0 0; 0 0 0];
        omega_b = 0;
        v_b = 0;
    else
        dot_r_d = (desired_r_vectors(:, i) - desired_r_vectors(:, i-1)) / dt;
        dot_R_d = (desired_R_rotms(:, :, i) - desired_R_rotms(:, :, i-1))./dt;
    end
    
    R_tilde = transpose(R_b) * desired_R_rotms(:, :, i);
    b_R_d = R_tilde;
    U = (k_omega*vee(expm(R_tilde)) + k_d*omega_b);
    omega_b = U;
    R_b = R_b*expm(dt*hat(omega_b));
    R_b_list(:, :, i) = R_b;
    
    u = transpose(R_b)*(dot_r_d + k_v*(desired_r_vectors(:, i) - r_b)) + k_d*v_b;
    v_b = u;
    r_b = r_b +dt*u;
    r_b_list(:, i) = r_b;
end

subplot(2,3,1);
plot(time_vec, desired_r_vectors(1,:), time_vec, r_b_list(1, :));
title('Axis X');
legend('Desired r', 'Body r');

subplot(2,3,2);
plot(time_vec, desired_r_vectors(2,:), time_vec, r_b_list(2, :));
title('Axis Y');
legend('Desired r', 'Body r');

subplot(2,3,3);
plot(time_vec, desired_r_vectors(3,:), time_vec, r_b_list(3, :));
title('Axis Z');
legend('Desired r', 'Body r');

eul_list = rotm2eul(R_b_list);
subplot(2,3,4);
plot(time_vec, desired_R_euls(:, 1), time_vec, eul_list(:, 1));
title('Roll');
legend('Desired r', 'Body r');

subplot(2,3,5);
plot(time_vec, desired_R_euls(:, 2), time_vec, eul_list(:, 2));
title('Pitch');
legend('Desired r', 'Body r');

subplot(2,3,6);
plot(time_vec, desired_R_euls(:, 3), time_vec, eul_list(:, 3));
title('Yaw');
legend('Desired r', 'Body r');

% hat 연산자 함수
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

% vee 연산자 함수
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
