function GPHMM_plot(SNPdatasource, resultssource, plotssource,configfile)
%this function is used to plots results of multiple samples.
%SNPdatasource: the directory containing all SNParray data
%resultssource: the directory containing all result files (one to one
%mapping assumed between SNParray file and result file
%plotssource: the directory that plots will be saved
if isdeployed
    if nargin<3
        error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.txt.'] );
    end
else
    if nargin<3
        error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.m.'] );
    end
end
current_version = '1.2';
disp(['GPHMM_plot (version ' current_version ') is loading...'])
tic
% if nargin == 4
% plot_mc = 0;
% end
datafilelist = dir(SNPdatasource);

%record all the time used for batch annotation
% disp(['-----Plot GPHMM results now, ' num2str(length(datafilelist)-2) ' SNP-array data files are found-----'])

if length(datafilelist)<3
    error('No SNP-array data files in the directory!');
else %now do batch annotation
    if nargin<4
        if length(datafilelist)>3
            disp('-----NO config file is found and GPHMM will plot ALL data files in directory-----')
            disp(['-----Plot GPHMM results now, ' num2str(length(datafilelist)-2) ' SNP-array data files are found-----'])

        else
            disp('-----NO config file is found and GPHMM will plot all data files-----')
            disp('-----Plot GPHMM results now, one SNP-array data file is found-----')
        end
        filename=cell(1,(length(datafilelist)-2));
        for i=3:length(datafilelist)
            filename{i-2}=datafilelist(i).name;
        end
        clear i
    else
        disp('-----Config file is found and GPHMM will ONLY plot the data files listed in config file-----')
        fp_conf=fopen(configfile,'r');
        if fp_conf==-1
            error('Can not open the config file, please check it again!');
        end
        filename=[];
        while 1
            confline=fgetl(fp_conf);
            if confline==-1
                break;
            end
            [temp, tumor_w1, tumor_w2]=strread(confline,'%s\t%f\t%f');
            filename=[filename {char(temp{1})}];

        end
        %         result=textscan(fp_conf,'%s');
        %         filename=result{1};
        fclose(fp_conf);
        clear result;
        if length(filename)>1
%             disp('-----NO config file is found and GPHMM will plot ALL data files in directory-----')
            disp(['-----Plot GPHMM results now, ' num2str(length(filename)) ' SNP-array data files are found-----'])

        else
%             disp('-----NO config file is found and GPHMM will plot all data files-----')
            disp('-----Plot GPHMM results now, ONE SNP-array data file is found-----')
        end
    end

    %----------------------load gcdata ------------------------%
%     fid = fopen(Gcdatafile,'r');
%     if fid == -1
%         error('Can not open gc file!');
%     end
%     results = textscan(fid,'%s %*s %*f %f','headerLines',1);
%     fclose(fid);
%     gc_pids = results{1}; %will be used throughout the procedure
%     %     gcs = (results{2}-41.6610)/4.7086; %will be used throughout the procedure
%     gcs = (results{2}-mean(results{2}))/std(results{2}); %will be used throughout the procedure
%     clear results;

    for fileindx = 1:length(filename)
        SNPdatafile = [SNPdatasource '/' filename{fileindx}];
        results = regexp(filename{fileindx},'^(.+)\.+.+','tokens','once');
        if isempty(results)
            fn_nosuffix = filename{fileindx};
        else
            fn_nosuffix = results{1};
            if ~isempty(strfind(fn_nosuffix,'.'))
                fn_nosuffix(strfind(fn_nosuffix,'.')) = '_';
            end
        end
        resultsfile = [resultssource '/' fn_nosuffix '.results'];
        plotsdir = [plotssource '/' fn_nosuffix];
        s = mkdir(plotsdir);
        if ~s %failed to make a directory
            error(['Can not make directory: ' plotsdir]);
        else
            GPHMM_plot_normalized_results(SNPdatafile,resultsfile,plotsdir,fn_nosuffix);
            disp ([num2str(fileindx) '. ' filename{fileindx} ' is done'])
        end
    end %for fileindx = 3:length(datafilelist)
end %if length(datafilelist)<3

t_all = toc;
disp(['-----Batch plotting is finished, totally ' num2str(t_all/60) ' minites were used-----'])
clear all