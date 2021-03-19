function varargout = flows_direction_similarity(V_1, V_2, num_modes)
% Uses average dot product of normalized flow vectors as measure of
% pattern similarity
for i = 1:num_modes
    
    V_1_x = V_1.vx(:, i);
    V_1_y = V_1.vy(:, i);
    V_1_z = V_1.vz(:, i);
    
    V_2_x = V_2.vx(:, i);
    V_2_y = V_2.vy(:, i);
    V_2_z = V_2.vz(:, i);
    
    V_1_mags = sqrt(V_1_x.*V_1_x + V_1_y.*V_1_y + V_1_z.*V_1_z);
    V_2_mags = sqrt(V_2_x.*V_2_x + V_2_y.*V_2_y + V_2_z.*V_2_z);
    
    V_1_x = V_1_x./V_1_mags;
    V_1_y = V_1_y./V_1_mags;
    V_1_z = V_1_z./V_1_mags;
    
    V_2_x = V_2_x./V_2_mags;
    V_2_y = V_2_y./V_2_mags;
    V_2_z = V_2_z./V_2_mags;
    
    dot_products = V_1_x.*V_2_x + V_1_y.*V_2_y + V_1_z.*V_2_z;
    
    varargout{i} = abs(mean(dot_products));
    
end
end
