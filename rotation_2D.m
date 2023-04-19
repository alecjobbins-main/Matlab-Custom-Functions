function [rotated_meshgrid_x, rotated_meshgrid_y] = rotation_2D(input_2D_coordinate_matrix, theta)
    %% Rotates a coordinate matrix (typically 2d matrix of standard x-y coordinates)
    % to some angle theta (radians)

    % Create rotation matrix
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    % Rotate your point(s)
    point = [3 5]'; % arbitrarily selected
    rotpoint = R*point;

end