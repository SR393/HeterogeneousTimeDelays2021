function window_data(data1, locs, start, finish, winlen, data_dt, window_dt, save_loc)
%Saves windows of length winlen taken from data at intervals window_dt in
%interval [start, finish]

if ~isfolder(save_loc)
    mkdir(save_loc)
end

ht = data_dt;
window_step = window_dt/data_dt;
start = start/data_dt;
finish = finish/data_dt;
winlen = winlen/data_dt;

for i = start:window_step:finish
    
    data = data1(i+1:i+winlen+1, :);
    save(strcat(save_loc, '\t', num2str(i*data_dt), '.mat'), 'data', 'locs', 'ht');
    
end
end