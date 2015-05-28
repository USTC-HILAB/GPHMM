function GPHMM(SNPdatasource,Pfbdatafile,Gcdatafile,Outputsource,genotypepath,sameformat,FileListfile)
%--------------------------------------------------------------------%
%------------------>       version 0.0.1       <---------------------
%--------------------------------------------------------------------%
%05/26/2010 by Ao
%This is the first version of the GPHMM method, which is based on the codes
%of GPHMM

%Input and output
%SNPdatasource: can be a directory or a file containing the path of every
%SNP array data
%Pfbdatafile: path of pdf file
%Gcdatafile: path of gc file
%Outputsource: directory of output files
%sourceflag: 1: directory otherwise: list file
%sameformat: if true for all SNPdatafile,then just map pbf and gc file
%into different chromosomes once.
global current_version
global NoSolutionFlag
current_version = '1.4';
sourceflag = 1;
% sameformat = 0;

if isdeployed
    if nargin<6
        error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.txt.'] );
    end
    if length(genotypepath) == 2 && genotypepath(1) == 'N' && genotypepath(2) == 'A'
        genotypepath = '';
    end
    if length(sameformat)~=1
        error('Please input 1 if the SNP array data files are in the same format and 0 otherwise')
    elseif  sameformat~='1' && sameformat~='0'
        error('Please input 1 if the SNP array data files are in the same format and 0 otherwise')
    end
    sameformat = str2double(sameformat);
else
    if nargin<6
        error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.m.'] );
    end
end
%parameters used in GPHMM
%===============================================
%---for cn up to 7---
depend_table_sc = [...
    %w1
    1 1 0.01 0.5;...
    2 1 1 1.0;...
    3 1 2 0.5;...
    4 1 2 1.0;...
    5 1 3 0.67;...
    6 1 3 1.0;...
    7 1 4 0.75;...
    8 1 4 0.5;...
    9 1 4 1.0;...
    10 1 5 0.8;...
    11 1 5 0.6;...
    12 1 5 1.0;...
    13 1 6 5/6;...
    14 1 6 4/6;...
    15 1 6 0.5;...
    16 1 6 1.0;...
    17 1 7 6/7;...
    18 1 7 5/7;...
    19 1 7 4/7;...
    20 1 7 1.0;...
    %     21 1 0.01 0.5...
    ];

depend_table = [...
    %w1
    1 1 0.01 0.5;...
    2 1 1 1.0;...
    3 0 2 0.5;...
    4 1 2 1.0;...
    5 1 3 0.67;...
    6 1 3 1.0;...
    7 1 4 0.75;...
    8 0 4 0.5;...
    9 1 4 1.0;...
    10 1 5 0.8;...
    11 1 5 0.6;...
    12 1 5 1.0;...
    1 2 0.01 0.5;...
    2 2 1 1.0;...
    3 0 2 0.5;...
    4 2 2 1.0;...
    5 2 3 0.67;...
    6 2 3 1.0;...
    7 2 4 0.75;...
    8 0 4 0.5;...
    9 2 4 1.0;...
    10 2 5 0.8;...
    11 2 5 0.6;...
    12 2 5 1.0;...
    1 3 0.01 0.5;...
    2 3 1 1.0;...
    3 3 2 0.5;...
    4 3 2 1.0;...
    5 3 3 0.67;...
    6 3 3 1.0;...
    7 3 4 0.75;...
    8 3 4 0.5;...
    9 3 4 1.0;...
    10 3 5 0.8;...
    11 3 5 0.6;...
    12 3 5 1.0;...
    13 3 6 5/6;...
    14 3 6 4/6;...
    15 3 6 0.5;...
    16 3 6 1.0;...
    17 3 7 6/7;...
    18 3 7 5/7;...
    19 3 7 4/7;...
    20 3 7 1.0;...
    %     21 3 0.01 0.5;...
    ];

%ABBBBB,AABBBB,AAABBB,BBBBBB
%ABBBBBB,AABBBBB,AAABBBB,BBBBBBB
%===============================================
screening_method = 2;
thres_EM = 1e-4;
max_iter1 = 100;
verbose1 = 0;
init_GPHMM_paras = [{[]},{[]},{[]},{[]},{[]},{[]},{[]},...
    {[]},{[]}]; %initial parameters will be assigned in the main function

%initialization of global variable
global data_lrr_sep
global data_baf_sep
global data_pb_sep
global data_gc_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global sin_or_mult %single clone sin_or_mult=0, else XX=1.
global tumor_range
sin_or_mult=0;

if sourceflag  == 1 % reading a directory
    disp(['GPHMM (version ' current_version ') is loading...'])
    if sameformat
        disp('Warning: you are assuming all SNP-array data files have EXACTLY the same format!')
        disp('It will lead to WRONG results if it is not the truth. You can set the parameter ')
        disp('"sameformat" to 0 to force SNP-id mapping on every data file.')
    end
    datafilelist = dir(SNPdatasource);
    if length(datafilelist)<3
        error('No SNP-array data files in the directory!');
    else %now do batch annotation
        
        if nargin<7
            disp('-----NO list file is found and GPHMM will analyse ALL data files in directory-----')
            filename=cell(1,(length(datafilelist)-2));
            for i=3:length(datafilelist)
                filename{i-2}=datafilelist(i).name;
            end
            tumor_range_table=[0.1*ones((length(datafilelist)-2),1) ones((length(datafilelist)-2),1)];
            clear i
        else
            disp('-----List file is found and GPHMM will ONLY analyse the data files in list-----')
            fp_conf=fopen(FileListfile,'r');
            if fp_conf==-1
                error('Can not open the list file, please check it again!');
            end
            filename=[];
            tumor_range_table=[];
            while 1
                confline=fgetl(fp_conf);
                if confline==-1
                    break;
                end
                [temp, tumor_w1, tumor_w2]=strread(confline,'%s\t%f\t%f');
                filename=[filename {char(temp{1})}];
                if (tumor_w1<0) ||(tumor_w1>1) ||(tumor_w2<0)||(tumor_w2>1)||(tumor_w2<tumor_w1)
                    disp(['Warning: tumor proportion range input is invalid: (' num2str(tumor_w1) ' , ' num2str(tumor_w2) ')'...
                        ', and they are modified to (0.1 , 1.0).']);
                    tumor_w1=0.1;
                    tumor_w2=1;
                end
                tumor_range_table=[tumor_range_table; [tumor_w1 tumor_w2]];
                
            end
            fclose(fp_conf);
            clear result;
        end
        %--------------perform GPHMM--------------------
        data_pb_sep = [];
        data_gc_sep = [];
        %this format works for Turin data set
        id_indx = 1;
        Chr_indx = 2;
        pos_indx = 3;
        baf_indx = 4;
        lrr_indx = 5;
        
        disp('GPHMM: Reading pfb file...');
        %--------------load pfbdata--------------------
        fid = fopen(Pfbdatafile,'r');
        if fid == -1
            error('Can not open pfb file!');
        end
        results = textscan(fid,'%s %*s %*f %f','headerLines',1);
        fclose(fid);
        pfb_pids = results{1}; %will be used throughout the procedure
        pfbs = results{2}; %will be used throughout the procedure
        clear results;
        
        %--------------load gcdata--------------------
        disp('GPHMM: Reading GC file...');
        fid = fopen(Gcdatafile,'r');
        if fid == -1
            error('Can not open gc file!');
        end
        results = textscan(fid,'%s %*s %*f %f','headerLines',1);
        fclose(fid);
        gc_pids = results{1}; %will be used throughout the procedure
        %         gcs = (results{2}-41.6610)/4.7086; %will be used throughout the procedure
        if std(results{2})==0
            gcs = ones(size(results{2})); % when all gc values are equal!
        else 
            gcs = (results{2}-mean(results{2}))/std(results{2}); %will be used throughout the procedure
        end
        clear results;
        
        %record all the time used for batch annotation
        t_all = 0;
        if length(filename)==1
            disp(['-----GPHMM batch annotation starts now, ONE SNP-array data file is found-----']);
        elseif length(filename)>1
            disp(['-----GPHMM batch annotation starts now, ' num2str(length(filename)) ' SNP-array data files are found-----']);
        end
        for fileindx = 1:length(filename)
            %clear global variables
            data_lrr_sep = [];
            data_baf_sep = [];
            
            gamma_sep = [];
            condi_probs_sep = [];
            condi_probs_fluct_sep = [];
            
            %record time cost
            tic
            tumor_range=tumor_range_table(fileindx,:);
            results = regexp(filename{fileindx},'^(.+)\.+.+','tokens','once');
            if isempty(results)
                fn_nosuffix = filename{fileindx};
            else
                fn_nosuffix = results{1};
                if ~isempty(strfind(fn_nosuffix,'.'))
                    fn_nosuffix(strfind(fn_nosuffix,'.')) = '_';
                end
            end
            %--------open result files --------------
            fid2 = fopen([Outputsource '/' fn_nosuffix '.results'],'w');
            if fid2 == -1
                error(['Can not open result file for ' filename{fileindx}]);
            end
            
            %--------------load SNPdata--------------------
            disp('GPHMM: Reading Tumor SNP-array file...');
            fid = fopen([SNPdatasource '/' filename{fileindx}],'r');
            if fid == -1
                error(['Can not open SNP array file: ' filename{fileindx}]);
            end
            temp=fgetl(fid);
            results = textscan(fid,'%s %f %f %f %f %*f %*f', 'treatAsEmpty', {'NA', 'na'});% only autosomes are considered here
            fclose(fid);
            
            Chr = results{Chr_indx};
            Chromosomes = reshape(intersect(unique(Chr),1:30),1,[]); % only use autosome i=reshape(unique(Chr),1,[]
            data_snpid_all = results{id_indx};
            data_pos_all = results{pos_indx};
            data_baf_all = results{baf_indx};
            data_lrr_all = results{lrr_indx};
            clear results temp;
            
            fprintf(1,'GPHMM: Total %d SNP probes are loaded from data file "%s".\n',length(data_lrr_all),filename{fileindx});
            if length(data_lrr_all)<10000
                fprintf(1,'Warning: the number of SNP probes loaded is too small, Pls check the format of data file!\n');
            end
            
            %use at least 30000 SNP for screening
            stepsize_ds = max(floor(length(Chr)/40000),10);
            
            if length(data_snpid_all) > length(Chr)
                data_snpid_all(end) = [];%remove last id in next row
            end
            if size(data_pos_all,1)>size(data_pos_all,2) %make sure it's 1byN
                data_pos_all = data_pos_all';
            end
            if size(data_baf_all,1)>size(data_baf_all,2) %make sure it's 1byN
                data_baf_all = data_baf_all';
            end
            if size(data_lrr_all,1)>size(data_lrr_all,2) %make sure it's 1byN
                data_lrr_all = data_lrr_all';
            end
            %preprocess LRR and BAF signals
            data_lrr_all = GPHMM_preprocess_lrr(data_lrr_all,9,1);
            data_baf_all(isnan(data_baf_all)) = 1;
            
            if (sum(data_baf_all>1)>0)
                disp(['Warning: certain BAF values are found more than 1 in file: "'  filename{fileindx} ...
                    '", and they are modified to 0.5 by deflaut.']);
                data_baf_all(data_baf_all>1)=0.5;
            end
            if (sum(data_baf_all<0)>0)
                disp(['Warning: certain BAF values are found less than 0 in file: "'  filename{fileindx} ...
                    '", and they are modified to 0.5 by deflaut.']);
                data_baf_all(data_baf_all<0)=0.5;
            end
            if ~sameformat || (sameformat&&(isempty(data_pb_sep)))
                data_pb_sep = [];
                data_gc_sep = [];
                %only load pfb and gc data when either data files are in different
                %formats or they are but need to load them for the first time
                %--------------map pfbdata--------------------
                disp('GPHMM: Preprocessing the data...');
                [tv,indx] = ismember(data_snpid_all,pfb_pids);
                data_pb_all = reshape(0.5*ones(size(Chr)),1,[]); %make sure it's 1byN
                data_pb_all(tv) = pfbs(indx(tv));
                if (sum(data_pb_all>1)>0)
                    disp(['Warning: certain PFB values are found more than 1 in file: "' Pfbdatafile ...
                        '", and they are modified to 0.5 by deflaut.']);
                    data_pb_all(data_pb_all>1)=0.5;
                elseif (sum(data_pb_all<0)>0)
                    disp(['Warning: certain PFB values are found less than 0 in file: "' Pfbdatafile ...
                        '", and they are modified to 0.5 by deflaut.']);
                    data_pb_all(data_pb_all<0)=0.5;
                end
                clear tv indx;
                
                %--------------map gcdata--------------------
                [tv,indx] = ismember(data_snpid_all,gc_pids);
                data_gc_all = reshape(0.0*ones(size(Chr)),1,[]); %make sure it's 1byN
                data_gc_all(tv) = gcs(indx(tv));
                clear tv indx;
                
                %-----divide into different chromosomes-----
                for i = Chromosomes
                    tv = ismember(Chr,i);
                    data_pb_sep = [data_pb_sep {data_pb_all(tv)}];
                    data_gc_sep = [data_gc_sep {data_gc_all(tv)}];
                end
                clear tv data_pb_all data_gc_all;
                
            end % if ~sameformat
            
            %-------divide into different chromosomes----------
            data_pos_sep = [];
            data_snpid_sep = [];
            data_bafA_sep = [];
            
            for i = Chromosomes
                tv = ismember(Chr,i);
                %LRR
                data_lrr_sep = [data_lrr_sep {data_lrr_all(tv)}];
                %BAF->mBAF
                temp = data_baf_all(tv);
                tvA=temp<0.5;
                data_bafA_sep = [data_bafA_sep {tvA}];
                clear tvA
                temp(temp<0.5) = 1-temp(temp<0.5);
                data_baf_sep = [data_baf_sep {temp}];
                %pos
                data_pos_sep = [data_pos_sep {data_pos_all(tv)}];
                %snpid
                data_snpid_sep = [data_snpid_sep {data_snpid_all(tv)}];
            end
            clear tv temp data_lrr_all data_baf_all data_pos_all...
                data_snpid_all Chr;
            
            disp('GPHMM: Computating the genotype (this may take a few minutes)...');
            %------------------ call GPHMM --------------------
            [LL_all,GPHMM_paras,best_indx,ds_info] = ...
                GPHMM_main(screening_method,init_GPHMM_paras,depend_table,depend_table_sc,stepsize_ds,thres_EM,verbose1,fn_nosuffix);
            
            %-------------- save ds_info ----------------------
            
            % eval(['save ' Outputsource '/' fn_nosuffix ' ds_info'])
            
            %-------------- process results ----------------------
            w_all = GPHMM_paras{5}{best_indx};
            if w_all(2)>0 %multiple clone
                [p_states,num_SNP,aCN,p_s,aCN_C1,aCN_C2,states_pre,state_seq,het_seq] = GPHMM_process_results_new(w_all,depend_table);
            else
                [p_states,num_SNP,aCN,p_s,aCN_C1,aCN_C2,states_pre,state_seq,het_seq] = GPHMM_process_results_new(w_all(1),depend_table_sc);
                Indicator_Score=GPHMM_postprocess(state_seq,het_seq,1-sum(w_all),GPHMM_paras{3}{best_indx},GPHMM_paras{4}{best_indx},...
                    sqrt(GPHMM_paras{7}{best_indx}(1)),sqrt(GPHMM_paras{7}{best_indx}(2)),sqrt(GPHMM_paras{6}{best_indx}));
            end
            
            listfile='LOG.txt';
            fp_list=fopen(listfile,'a+');
            if fp_list==-1
                warning('Can not open the file: "%s"!\n',listfile);
            else
                if fileindx==1
                    fprintf(fp_list,'Version:\tDate\tTime\tSample\tTumor content\tBaseline shift\tGC coef.\tACN\tBAF_Het Sigma\tLRR sigma\n');
                end
                fprintf(fp_list,'GPHMM%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n',...
                    current_version,datestr(clock,'mmm-dd-yyyy'),datestr(clock,'HH:MM:SS'),fn_nosuffix,1-sum(w_all),GPHMM_paras{3}{best_indx},...
                    GPHMM_paras{4}{best_indx},aCN_C1,sqrt(GPHMM_paras{7}{best_indx}(1)),sqrt(GPHMM_paras{6}{best_indx}));
            end
            fclose(fp_list);
            disp('GPHMM: Plotting the result...');
            GPHMM_plot_tumor(current_version,Outputsource,fn_nosuffix,Chromosomes,data_pos_sep,1-sum(w_all),GPHMM_paras{3}{best_indx},...
                GPHMM_paras{4}{best_indx},state_seq,data_bafA_sep,Indicator_Score);
           
            clear LL_all gamma_all ds_info
            
            %-------------- output summary of the results--------------
            disp('GPHMM: Writing the result file...');
            fprintf(fid2,'---------------------------------------------------------------\n');
            fprintf(fid2,['             Summary of GPHMM results (version ' current_version ')          \n']);
            if NoSolutionFlag
                fprintf(fid2,'Warning: Prediction results may be inaccurate due to the failure\n');
                fprintf(fid2,'in finding optimal initial global parameters!\n');
            end
            
            fprintf(fid2,'General information of this cancer sample:                      \n');
            fprintf(fid2,'   LRR correction factor: %6.4f\n',GPHMM_paras{3}{best_indx});
            fprintf(fid2,'   GC coefficient: %6.4f\n',GPHMM_paras{4}{best_indx});
            fprintf(fid2,'   Proportion of abnormal cells in the sample: %6.4f\n',1-sum(w_all));
            fprintf(fid2,'   Standard deviation of LRR signal: %6.4f\n',sqrt(GPHMM_paras{6}{best_indx}));
            fprintf(fid2,'   Standard deviation of BAF signal: %6.4f\n',sqrt(GPHMM_paras{7}{best_indx}(1)));
            fprintf(fid2,'   Proportion of all abnormal chromosomal regions: %6.4f\n',1-p_s(4));
            fprintf(fid2,'   Estimated average cancer DNA index: %6.4f\n',aCN(3)/2);
           
            fprintf(fid2,'---------------------------------------------------------------\n');
            fprintf(fid2,'\n');
            clear w_all p_states num_SNP aCN p_s p_AACR aCN_C1 p_ACR_C1 aCN_C2 p_ACR_C2
            
            %--------------output state assignment in segments--------------
            fprintf(fid2,'Chr\tStartPos\tEndPos\tState\tScore\tCN\tAI\tLength\tStartSNPid\tEndSNPid\tStartIndx\tEndIndx\n');
            cn = [0 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 ];
            Bcn= [0 1 1 2 2 3 3 2 4 4 3 5 5 4 3 6 6 5 4 7 ];
            AI_mapping = [1 1 3 2 3 2 3 3 2 3 3 2 3 3 3 2 3 3 3 2 1];
            Agtype=[{''} {'A'} {'AA'} {'AAA'} {'AAAA'} {'AAAAA'} {'AAAAAA'} {'AAAAAAA'}];
            Bgtype=[{''} {'B'} {'BB'} {'BBB'} {'BBBB'} {'BBBBB'} {'BBBBBB'} {'BBBBBBB'}];
            states_pre = states_pre{best_indx};
            %             MC_status = ['1' '2' 'B'];
            for i=1:size(states_pre,1)
                Chr_i = states_pre(i,1);
                St_i = states_pre(i,2);
                Ed_i = states_pre(i,3);
                S_C1 = states_pre(i,4);
                %                 S_C2 = states_pre(i,5);
                %                 MC_indx = states_pre(i,6);
                fprintf(fid2,'%d\t%d\t%d\t%d\t%1.4f\t%d\t%d\t%d\t%s\t%s\t%d\t%d\n',Chromosomes(Chr_i),data_pos_sep{Chr_i}(St_i),...
                    data_pos_sep{Chr_i}(Ed_i),S_C1,mean(gamma_sep{Chr_i}(S_C1,St_i:Ed_i)),cn(S_C1),AI_mapping(S_C1),Ed_i-St_i+1,data_snpid_sep{Chr_i}{St_i},data_snpid_sep{Chr_i}{Ed_i},St_i,Ed_i);
            end
            fclose(fid2);
            clear Chr_i St_i Ed_i S_C1  states_pre
            
            disp('GPHMM: Writing the genotype file...');
            if ~isnan(genotypepath)
                fpG=fopen([genotypepath '/' fn_nosuffix '.Gtype'],'w');
                if fpG==-1
                    error('Can not open the Genotype result file! Please check again. ');
                end
                
                fprintf(fpG,'Name\tChr\tposition\tCN\tB allele CN\tState\tGenotype\n');
                for i=1:length(state_seq)
                    
                    cn_seq=cn(state_seq{i});
                    Bcn_seq=Bcn(state_seq{i});
                    tv=(data_bafA_sep{i}==1 & het_seq{i});
                    Bcn_seq(tv)=cn_seq(tv)-Bcn_seq(tv);
                    tv=(data_bafA_sep{i}==0 & ~het_seq{i});
                    Bcn_seq(tv)=cn_seq(tv);
                    tv=(data_bafA_sep{i}==1 & ~het_seq{i});
                    Bcn_seq(tv)=0;
                    for j=1:length(state_seq{i})
                        fprintf(fpG,'%s\t%d\t%d\t%d\t%d\t%d',data_snpid_sep{i}{j},Chromosomes(i),data_pos_sep{i}(j),cn_seq(j),Bcn_seq(j),state_seq{i}(j));
                        if state_seq{i}(j)==1
                            fprintf(fpG,'\tNA\n');
                        else
                            fprintf(fpG,'\t%s',Agtype{cn_seq(j)-Bcn_seq(j)+1});
                            fprintf(fpG,'%s\n',Bgtype{Bcn_seq(j)+1});
                        end
                    end
                    clear cn_seq Bcn_seq
                end
                fclose(fpG);
                
            end
            clear GPHMM_paras best_indx data_pos_sep data_snpid_sep Chromosomes
            
            disp('GPHMM: Tumor sample analysis done!');
            
            %final display a brief report
            t = toc;
            t_all = t_all+t;
            disp ([num2str(fileindx) '. ' filename{fileindx} ' is done, time used: ' num2str(t)])
            
        end %for fileindx = 3:length(datafilelist)
    end % if length(datafilelist)<3
    
    disp(['-----GPHMM is finished, totally ' num2str(t_all/60) ' minites were used-----'])
    clear all
end %if sourceflag  == 1
