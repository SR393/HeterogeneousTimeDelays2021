function varargout = flows_magnitude_similarity(V_1, V_2, num_modes)
% Uses sum-of-squares-of-differences measure on flow vector magnitudes 
% as measure of pattern dissimilarity
for i = 1:num_modes
    
    V_1_x = V_1.vx(:, i);
    V_1_y = V_1.vy(:, i);
    V_1_z = V_1.vz(:, i);
    
    V_2_x = V_2.vx(:, i);
    V_2_y = V_2.vy(:, i);
    V_2_z = V_2.vz(:, i);
    
    V_1_mags = sqrt(V_1_x.*V_1_x + V_1_y.*V_1_y + V_1_z.*V_1_z);
    V_2_mags = sqrt(V_2_x.*V_2_x + V_2_y.*V_2_y + V_2_z.*V_2_z);
    
    varargout{i} = sum((V_2_mags - V_1_mags).^2);
    
end
end
