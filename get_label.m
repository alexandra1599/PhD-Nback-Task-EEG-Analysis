function label = get_label(streams)

label = cell(1,size(streams,2)); 
l1 = cell(1,63);

for j = 1:1:size(streams,2)

for i = 1:1:size(streams{1,j},2)
    data_type = string(streams{1,j}{1,i}.info.type);
    format long g

        if strcmp(data_type, 'EEG')  % Markers data stream
            l = streams{1,j}{1,i}.info.desc.channels.channel;
            for k=1:1:size(l,2)
                label{1,k} = l{1,k}.label;
            end
            
       end
end

r=7;
label = label(1:end-r);  % Removes the last 7 columns
end



end