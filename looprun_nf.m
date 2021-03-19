function looprun_nf(path_ending)

datafolder = strcat('Data\\', path_ending);
if ~isfolder(datafolder)
    mkdir(datafolder)
end
jsonfilefolder = strcat('jsons\\', path_ending);
if ~isfolder(jsonfilefolder)
    mkdir(jsonfilefolder)
end
savefigloc = strcat('Output\\Figures\\', path_ending);
if ~isfolder(savefigloc)
    mkdir(savefigloc)
end
outdatafolder = strcat('Output\\Data\\', path_ending);
if ~isfolder(outdatafolder)
    mkdir(outdatafolder)
end

files = dir(strcat(datafolder, '\*.mat'));

for k = 1:length(files)
    filename = files(k).name;
    filename = filename(1:length(filename)-4);
    write_json(path_ending, filename);
    run_nf(jsonfilefolder, filename, savefigloc);
end
end
    