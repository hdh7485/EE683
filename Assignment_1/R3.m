k_p = 0.1;

r_b = [0; 0; 0];
r_d = [0; 10; 5];

% initialize rpy 0 0 0
start_R = [0 0 0];
rot0 = eul2rotm(start_R);
R_b = rot0;

% start r point
start_r = [0; 0; 0];
% end r point
end_r = [0; 10; 5];
% interpolate desired point
desired_r_vectors = interpolator(start_r, end_r, 10)

dot_r_d = transpose([0 0 0]);

u = transpose(R_b)*(dot_r_d + k_p*(r_d - r_b))

function points = interpolator(start_point, end_point, numbers)
q1 = linspace(start_point(1), end_point(1), numbers);
q2 = linspace(start_point(2), end_point(2), numbers);
q3 = linspace(start_point(3), end_point(3), numbers);
points = [transpose(q1), transpose(q2), transpose(q3)];
end