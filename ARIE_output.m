function ARIE_output(x, y, z)
% Output the ARIE data as .txt format

%% Set file name
curr_time = clock;
time_str = [num2str(curr_time(1)) '-' num2str(curr_time(2)) '-' ...
            num2str(curr_time(3)) '-ARIE_OUTPUT'];
log_filename = [time_str '.txt'];        
fin = fopen(log_filename, 'a', 'n', 'UTF-8');
% numeric format
format_double = '%4.6f\t';

%% Write to the .txt file
for i = 1:length(x)
    for j = 1:length(y)
       
        if ~isnan(z(j, i))
            s = [num2str(x(i),format_double) '   '...
                     num2str(y(j),format_double) '   '...
                     num2str(z(j, i),format_double)];

            fprintf(fin,'%s\n', s);
        else
            continue;
        end
        
    end
end

%% Close the file
fclose(fin);