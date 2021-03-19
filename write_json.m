function write_json(data_loc, filename)
%Replaces template.json with json file specific to a particular
%neural-flows analysis
nfs_dir = '';
data_dir = strcat(nfs_dir, 'Data\\');
jsons_dir = strcat(nfs_dir, 'jsons\\');

if ~isfolder(strcat(jsons_dir, data_loc))
    mkdir(strcat(jsons_dir, data_loc))
end

template_json = strcat(nfs_dir, 'template.json');
fid = fopen(template_json, 'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

%replace .json paths and filenames
A{6} = sprintf('%s', strcat('            "name": "', filename, '_in.json",'));
A{7} = sprintf('%s', strcat('            "dir":' ,' "', jsons_dir, data_loc, '"'));
A{12} = sprintf('%s', strcat('            "name": "', filename, '_out.json",')); 
A{13} = sprintf('%s', strcat('            "dir":' ,' "', jsons_dir, data_loc, '"}'));
A{18} = sprintf('%s', strcat('      "dir_out": "', nfs_dir, 'Output\\Data\\', data_loc, '",'));
A{19} = sprintf('%s', strcat('      "dir_tmp":','"/', filename, '",'));

%replace input and output data paths and filenames
A{33} = sprintf('%s', strcat('      "label":', ' "', filename, '_interp",'));
A{59} = sprintf('%s', strcat('      "label":', ' "', filename, '_flows"'));
A{152} = sprintf('%s', strcat('      "dir": "', data_dir, data_loc, '",')); 
A{153} = sprintf('%s', strcat('      "name":', ' "', filename, '.mat"'));

fid = fopen(strcat(jsons_dir, data_loc, '\\', filename, '.json'), 'w+');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
