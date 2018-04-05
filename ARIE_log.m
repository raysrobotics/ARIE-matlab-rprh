function ARIE_log(param, options, flag)
%  Log the current params, options and flag in a .txt file

%% Set file name
curr_time = clock;
time_str = [num2str(curr_time(1)) '-' num2str(curr_time(2)) '-' ...
            num2str(curr_time(3)) '-ARIE_LOG'];
log_filename = [time_str '.txt'];        
fin = fopen(log_filename, 'a', 'n', 'UTF-8');

format_double = '%4.4f\t';
format_int = '%5d\t';

%% Write to the .txt file
s1 = ['x:' num2str(param(1),format_double)...
      'y:' num2str(param(2),format_double)...
      'z:' num2str(param(3),format_double)];
s2 = ['theta_x:' num2str(options.theta_x,format_double)...
      'theta_y:' num2str(options.theta_y,format_double)];
s3 = ['outrange:' num2str(flag.outrange,format_int)...
      'inter:'    num2str(flag.inter   ,format_int)...
      'contact:'  num2str(flag.contact ,format_int)];

fprintf(fin,'%s\t%s\n%s\n\n', s1,s2,s3);

%% Close the file
fclose(fin);