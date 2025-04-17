function [streams] = load_project_data(ID,session,runs,type,N)

if ID == 321
    ID = 'Subject 321';
    k=321;
elseif ID == 322
    ID = 'Subject 322';
    k=322;
elseif ID == 323
    ID = 'Subject 323';
    k=323;
elseif ID == 621
    ID = 'Subject 621';
    k=621;
elseif ID == 622
    ID = 'Subject 622';
    k=622;
elseif ID == 623
    ID = 'Subject 623';
    k=623;
elseif ID == 921
    ID = 'Subject 921';
    k=921;
elseif ID == 922
    ID = 'Subject 922';
    k=922;
elseif ID == 923
    ID = 'Subject 923';
    k=923;
elseif ID == 1221
    ID = 'Subject 1221';
    k=1221;
elseif ID == 1222
    ID = 'Subject 1222';
    k=1222;
elseif ID == 1223
    ID = 'Subject 1223';
    k=1223;
elseif ID == 1224
    ID = 'Subject 1224';
    k=1224;
elseif ID == 1321
    ID = 'Subject 1321';
    k=1321;
elseif ID == 1322
    ID = 'Subject 1322';
    k=1322;
elseif ID == 1323
    ID = 'Subject 1323';
    k=1323;
end


base_path = fullfile('/Users/alexandra/Desktop/PhD/Project/Experiment/Data', ID,session);
streams = cell(1, runs);

if (strcmp(session, 'Relaxation'))
    load('ErrP_cap_chan_file.mat');
    if (strcmp(type,'EOG'))
      for i = 1:1:runs
            fprintf('Session is EOG \n'); % Print the session type
            filename = sprintf('EOG_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path,type, filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
      end  

    elseif (strcmp(type,'Eyes Closed pre'))
        for i = 1:1:runs
            fprintf('Session is Eyes Closed pre Relaxation\n'); % Print the session type
            filename = sprintf('ECpre_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path,type, filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  

    elseif (strcmp(type, 'Eyes Closed post'))
       fprintf('Session is Eyes Closed post Relaxation \n'); % Print the session type

        for i = 1:1:runs
            filename = sprintf('ECpost_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path,type, filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  
    elseif(strcmp(type,'Relax'))
        fprintf('Session is Relaxation \n'); % Print the session type

        for i = 1:1:runs
            filename = sprintf('relax_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path, type,filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  

    elseif any(strcmp(type, {'Nback', 'Nback + relax'}))
        if N == 0
            pathn = ('N0');
        elseif N == 1
            pathn = ('N1');
        elseif N == 2
            pathn = ('N2');
        elseif N == 3
            pathn = ('N3');
        end
        fprintf('Session is %s with N = %s/n',session,N); % Print the session type

        for i = 1:1:runs
            filen = fullfile(base_path,type,pathn);
            filename = sprintf('nb_run-%03d_eeg.xdf',i);
            file_path = fullfile(filen, filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end
    end


elseif (strcmp(session, 'tACS'))
    load("ch32Locations.mat");
 if (strcmp(type,'EOG'))
      for i = 1:1:runs
            fprintf('Session is EOG \n'); % Print the session type
            filename = sprintf('EOG_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path, type,filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
      end 

 elseif (strcmp(type,'Eyes Closed find tACS'))
        for i = 1:1:runs
            fprintf('Session is Eyes Closed find tACS \n'); % Print the session type
            filename = sprintf('ECfind_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path,type, filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end 

    elseif (strcmp(type,'Eyes Closed pre nothing'))
        for i = 1:1:runs
            fprintf('Session is Eyes Closed pre Nothing \n'); % Print the session type
            filename = sprintf('ECpreN_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path, type,filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  

    elseif (strcmp(type, 'Eyes Closed post nothing'))
        for i = 1:1:runs
            fprintf('Session is Eyes Closed post Nothing\n'); % Print the session type
            filename = sprintf('ECpostN_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path, type,filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  
     elseif (strcmp(type, 'Eyes Closed pre tACS'))
        for i = 1:1:runs
            fprintf('Session is Eyes Closed pre tACS\n'); % Print the session type
            filename = sprintf('ECpretACS_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path, type,filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  
    elseif (strcmp(type, 'Eyes Closed post tACS'))
        for i = 1:1:runs
            fprintf('Session is Eyes Closed post tACS\n'); % Print the session type
            filename = sprintf('ECpostACS_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path,type, filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  
    elseif(strcmp(type,'Nothing'))
        for i = 1:1:runs
            fprintf('Session is Nothing\n'); % Print the session type
            filename = sprintf('Nothing_run-%03d_eeg.xdf',i);
            file_path = fullfile(base_path, type,filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end  

    elseif any(strcmp(type, {'Nback', 'Nback + tACS'}))
        if N == 0
            pathn = ('N0');
        elseif N == 1
            pathn = ('N1');
        elseif N == 2
            pathn = ('N2');
        elseif N == 3
            pathn = ('N3');
        end
        fprintf('Session is %s with N = %s/n',session,N); % Print the session type

        for i = 1:1:runs
            
            filen = fullfile(base_path,type,pathn);
            filename = sprintf('nb_run-%03d_eeg.xdf',i);
            file_path = fullfile(filen, filename);
            streams{i} = loadxdf(file_path);  % This is a cell array of 3 elements per file
            fprintf('Loaded file: %s\n', filename); % Print the loaded file
        end
  end

end




end


   
