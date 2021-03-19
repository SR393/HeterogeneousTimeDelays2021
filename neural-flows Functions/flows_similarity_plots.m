function flows_similarity_plots(path, h)
%Plots of flows_direction_similarity and flows_magnitude_similarity
if nargin < 2
    h = 1;
end

step = 5;
start = 0;
stop = 1000;
time = start:step:stop;
p = length(time);
num_modes = 1;
dir_1_matr = zeros(num_modes, p);
dir_matr = zeros(num_modes, p - length(h));

mag_matr = zeros(num_modes, p);

for i = 1:p
    
    i
    flows_start = load(strcat('Output\\Data\\', path, '\\', 't', num2str(start), '_flows.mat'));
    flows_start = flows_start.flow_modes.V;
    flows_i = load(strcat('Output\\Data\\', path, '\\', 't', num2str(start+i*step), '_flows.mat'));
    flows_i = flows_i.flow_modes.V;
    flows_imin1 = load(strcat('Output\\Data\\', path, '\\', 't', num2str(start+step*(i-1)), '_flows.mat'));
    flows_imin1 = flows_imin1.flow_modes.V;
    
    [s1] = flows_direction_similarity(flows_start, flows_i, num_modes);
    dir_1_matr(:, i) = [s1];
    [s1] = flows_magnitude_similarity(flows_imin1, flows_i, num_modes);
    mag_matr(:, i) = [s1];
    
    if i > length(h)
        
        flows_t = load(strcat('Output\\Data\\', path, '\\', 't', num2str(start+step*(i-1)), '_flows.mat'));
        flows_t = flows_t.flow_modes.V;
        weighted_avg.vx = h(1)*flows_t.vx;
        weighted_avg.vy = h(1)*flows_t.vy;
        weighted_avg.vz = h(1)*flows_t.vz;
        
        for j = 2:length(h)
            
            flows_t = load(strcat('Output\\Data\\', path, '\\', 't', num2str(start+step*(i-j)), '_flows.mat'));
            flows_t = flows_t.flow_modes.V;
            weighted_avg.vx = weighted_avg.vx + h(j)*flows_t.vx;
            weighted_avg.vy = weighted_avg.vy + h(j)*flows_t.vy;
            weighted_avg.vz = weighted_avg.vz + h(j)*flows_t.vz;
            
        end
        [s1] = flows_direction_similarity(weighted_avg, flows_i, num_modes);
        dir_matr(:, i - length(h)) = [s1];
    end
end
figures_path = strcat('Output\\Figures\\', path);
if ~isfolder(figures_path)
    mkdir(figures_path)
end

plot(step*length(h):step:stop, dir_matr)
savefig(strcat(figures_path, '\\', 'Direction_Similarity.fig'))

plot(time, mag_matr)
savefig(strcat(figures_path, '\\', 'Magnitude_Difference.fig'))

plot(time, dir_1_matr)
savefig(strcat(figures_path, '\\', 'Direction1_Similarity.fig'))

end
