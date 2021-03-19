function del_tmstamps(path)

path = strcat('Output\\Data\\', path);

files = dir(strcat(path, '\*.mat'));

for j = 1:length(files)
        filename = files(j).name;
        filename1 = filename(1:length(filename)-33);
        filename1 = strcat(filename1, '.mat');
        movefile(strcat(path, '\', filename), strcat(path, '\', filename1))
end

end
