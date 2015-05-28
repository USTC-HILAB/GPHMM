function GPHMM_plot_normalized_results(SNPdatafile,resultsfile,plotsdir,barcode)
%this function is used to plot 22 figures given the SNP array data as well
%as the results from NormCNV, all these figures are stored in a specified
%directory
% if nargin==6
% plot_mc = 0;
% end
%%%
%----------------------read results  ------------------------%
fid = fopen(resultsfile,'r');
if fid == -1
    error(['Can not open result file: ' resultsfile]);
end

%get estimated global parameters from the first row of the result file
lrr_shift = [];
% gc_coef = [];
w = [];

while 1
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'StartPos')),break,end 
    %LRR baseline shift
    result1 = regexp(tline,'LRR correction factor:\s*(\S+)','tokens','once');
    result2 = regexp(tline,'LRR baseline shift:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        lrr_shift = [lrr_shift str2num(result1{1})];
    elseif ~isempty(result2)
        lrr_shift = [lrr_shift str2num(result2{1})];
    end
    %GC coefficient
    result1 = regexp(tline,'GC coefficient:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        gc_coef = str2num(result1{1});
    end
    %w1 and w2
    result1 = regexp(tline,'Proportion of abnormal cells in the sample:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        w = [w str2num(result1{1})];
    end
end
%report errors if these values are not parsed successfully
if isempty(lrr_shift)
    error(['Can not read estimated LRR shift from ',resultsfile]);
end
if isempty(gc_coef)
    error(['Can not read estimated GC coefficient from ',resultsfile]);
end
if length(w)<1
    error(['Can not read estimated tumor proportion from ',resultsfile]);
end

%then read the results
results = textscan(fid,'%d %*d %*d %d %*f %*d %*d %*d %*s %*s %d %d','headerLines',0, 'treatAsEmpty', {'NA', 'na'});
fclose(fid);
chr_seg = results{1};
state_C1_seg = results{2};
% state_C2_seg = results{3};
% mc_seg = ones(size(results{4}));
start_seg = results{3};
end_seg = results{4};
clear results;
%
%----------------------read array data ------------------------%
fid = fopen(SNPdatafile,'r');
if fid == -1
    error(['Can not open SNP array file: ' SNPdatafile]);
end
temp=fgetl(fid);
results = textscan(fid,'%s %f %f %f %f %*f %*f', 'treatAsEmpty', {'NA', 'na'});% only autosomes are considered here
fclose(fid);

%this format works for Turin data set
id_indx = 1;
Chr_indx = 2;
pos_indx = 3;
baf_indx = 4;
lrr_indx = 5;

Chr = results{Chr_indx};
data_snpid_all = results{id_indx};
if length(data_snpid_all) > length(Chr)
    data_snpid_all(end) = [];%remove last id in next row
end
data_pos_all = results{pos_indx};
data_baf_all = results{baf_indx};
data_lrr_all = results{lrr_indx};
clear results temp;

%--------------map gcdata--------------------
% [tv,indx] = ismember(data_snpid_all,gc_pids);
% data_gc_all = reshape(zeros(size(Chr)),1,[]); %make sure it's 1byN
% data_gc_all(tv) = gcs(indx(tv));

%--set median of LRR to 0 and replace NaNs in data--
if size(data_baf_all,1)>size(data_baf_all,2) %make sure it's 1byN
    data_baf_all = data_baf_all';
end
if size(data_lrr_all,1)>size(data_lrr_all,2) %make sure it's 1byN
    data_lrr_all = data_lrr_all';
end
%set median of LRR to 0 and replace NaNs in data
tv = isnan(data_lrr_all);
temp = median(data_lrr_all(~tv));
data_lrr_all(~tv) = data_lrr_all(~tv)-temp;
data_lrr_all(tv) = 0;
clear tv temp;
data_baf_all(isnan(data_baf_all)) = 0;


%
%---for cn up to 7---
depend_table = [...
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
    21 1 0.01 0.5];
CN_mapping = depend_table(:,3)';
Muc_mapping = depend_table(:,4)';
AI_mapping = [1 1 3 2 3 2 3 3 2 3 3 2 3 3 3 2 3 3 3 2 1];

% amplicon_info = [...
% 7   55054219	55242525;...
% 8   19840862	19869050;...
% 8   43200001    48100000;...
% 8	128409679	129230679;...
% 11  51400001    56400000;...
% 11	69021739 69321739;...
% 14  104512893	106012893;...
% 17	35011111 35237111;...
% 17  22100001 23200000;...
% 17  35743196 35903196;...
% 20	51465030	51785030;...
% 21	34562717	35862717;...
% ];
%--------------- plot figures ---------------------%

for i=reshape(intersect(unique(Chr),1:30),1,[]) %onlyl include autosome reshape(unique(Chr),1,[])
    %process SNParray data
    tv = ismember(Chr,i);
    data_lrr = data_lrr_all(tv);
    %imputation with midean filter and remove GC content bias
    data_baf = data_baf_all(tv);
    data_pos = data_pos_all(tv);
    
    indx1 = find(chr_seg==i);
    
    %plot
    clf reset
    marker_size = 3;
    %---plot LOH---
    subplot(4,1,1)
    hold on
    set(gca,'YGrid','on')
    axis ([-Inf Inf 0.9 3.1])
    set(gca,'YTick',[1:3])
    set(gca,'YTickLabel',{'Del','LOH','Het'});
    for j=reshape(indx1,1,[])
%         if mc_seg(j) == 1;
            AI = AI_mapping(state_C1_seg(j));
            line_style = 'r-';
%         elseif mc_seg{j} == '2';
%             AI = AI_mapping(state_C2_seg(j));
%             line_style = 'b-';
%         else %both
%             AI = AI_mapping(state_C2_seg(j));
%             line_style = 'm-';
%         end
        plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[AI AI], line_style, 'LineWidth',marker_size);
%         if plot_mc
%             if mc_seg(j) == 1
%                 if AI == 3
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[3 3], 'm-', 'LineWidth',marker_size);
%                 else
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[3 3], 'b-', 'LineWidth',marker_size);
%                 end
%             elseif mc_seg{j} == '2'
%                 if AI == 3
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[3 3], 'm-', 'LineWidth',marker_size);
%                 else
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[3 3], 'r-', 'LineWidth',marker_size);
%                 end
%             end
%         end
    end
    %replace '_' with '-' in barcode to display it correctly
    barcode_m = barcode;
    tmp = strfind(barcode_m,'_');
    barcode_m(tmp) = '-';
    title (['Chromosome ' num2str(i) ', ' barcode_m])

    %---plot CN---
    subplot(4,1,2)
    hold on
    set(gca,'YGrid','on')
    %     axis ([-Inf Inf -0.05 7.05])
    axis ([-Inf Inf -0.1 7.1])
    set(gca,'YTick',[0:1:7])
    set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
    ylabel('CN');
    for j=reshape(indx1,1,[])
%         if mc_seg(j) == 1;
            CN = CN_mapping(state_C1_seg(j));
%             uCN = CN_mapping(state_C1_seg(j))*Muc_mapping(state_C1_seg(j));
            line_style = 'r-';
%         elseif mc_seg{j} == '2';
%             CN = CN_mapping(state_C2_seg(j));
% %             uCN = CN_mapping(state_C2_seg(j))*Muc_mapping(state_C2_seg(j));
%             line_style = 'b-';
%         else %both
%             CN = CN_mapping(state_C2_seg(j));
% %             uCN = CN_mapping(state_C2_seg(j))*Muc_mapping(state_C2_seg(j));
%             line_style = 'm-';
%         end
        plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[CN CN], line_style, 'LineWidth',marker_size);
        %plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[uCN uCN], line_style, 'LineWidth',marker_size*0.3);
%         if plot_mc
%            if mc_seg(j) == 1
%                 if CN == 2
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[2 2], 'm-', 'LineWidth',marker_size);
%                 else
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[2 2], 'b-', 'LineWidth',marker_size);
%                 end
%             elseif mc_seg{j} == '2'
%                 if CN == 2
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[2 2], 'm-', 'LineWidth',marker_size);
%                 else
%                     plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[2 2], 'r-', 'LineWidth',marker_size);
%                 end
%             end
%         end
    end

    marker_size = 3;
    %---plot BAF---
    subplot(4,1,3)
    plot(data_pos,data_baf,'b.', 'MarkerSize',marker_size)
    hold on
    for j = 0:0.25:1
        plot([data_pos(1) data_pos(end)],[j j],'k-','LineWidth',0.5)
    end
    %plot expected BAF mean values
    for j=reshape(indx1,1,[])
        if state_C1_seg(j) == 21 %total deletion
            baf_mean = 0.5;
        else
%             if w(2)>0 %multiple clones
%                 temp1 = w(1)*CN_mapping(state_C1_seg(j))*Muc_mapping(state_C1_seg(j))+w(2)*CN_mapping(state_C2_seg(j))*...
%                     Muc_mapping(state_C2_seg(j))+(1-w(1)-w(2))*1;
%                 temp2 = w(1)*CN_mapping(state_C1_seg(j))+w(2)*CN_mapping(state_C2_seg(j))+(1-w(1)-w(2))*2;
%             else %single clone
                temp1 = w(1)*CN_mapping(state_C1_seg(j))*Muc_mapping(state_C1_seg(j))+(1-w(1))*1;
                temp2 = w(1)*CN_mapping(state_C1_seg(j))+(1-w(1))*2;
%             end
            baf_mean = temp1/temp2;
        end
        plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[baf_mean baf_mean],'r-','LineWidth',1);
        plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[1-baf_mean 1-baf_mean],'r-','LineWidth',1);
    end
    ylabel('BAF');
    axis ([-Inf Inf 0 1])
    set(gca,'YTick',[0 0.5 1])

    %---plot LRR---
    subplot(4,1,4)
    plot(data_pos,data_lrr,'b.', 'MarkerSize',marker_size)
    hold on
    for j = -1:0.5:1
        plot([data_pos(1) data_pos(end)],[j j],'k-','LineWidth',0.5)
    end
    %plot expected LRR mean values
    for j=reshape(indx1,1,[])
        if state_C1_seg(j) == 21 %total deletion
            lrr_mean = -3.5;
        else
%             if w(2)>0 %multiple clones
%                 lrr_mean = 2*log10((w(1)*CN_mapping(state_C1_seg(j))+w(2)*CN_mapping(state_C2_seg(j))+(1-w(1)-w(2))*2)/2)+lrr_shift;
%             else  %single clone
                lrr_mean = 2*log10((w(1)*CN_mapping(state_C1_seg(j))+(1-w(1))*2)/2)+lrr_shift;
%             end
        end
        plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[lrr_mean lrr_mean],'r-','LineWidth',1);
    end
    ylabel('LRR');
    axis ([-Inf Inf -2.5 1.5])
    set(gca,'YTick',[-2 -1 0 1 2])


    %plot amplicon regions
%     if ~isempty(amplicon_info)
%         indx2 = find(amplicon_info(:,1)==i);
%         if ~isempty(indx2)
%             for k = 1:length(indx2)
%                 subplot(4,1,1)
%                 plot([amplicon_info(indx2(k),2) amplicon_info(indx2(k),2)],[0.95 3.05],'k:','LineWidth',0.5) % line for start postion
%                 plot([amplicon_info(indx2(k),3) amplicon_info(indx2(k),3)],[0.95 3.05],'k--.','LineWidth',0.5)% line fore end postion
%                 subplot(4,1,2)
%                 % %                 plot([amplicon_info(indx2(k),2) amplicon_info(indx2(k),2)],[-0.05 5.05],'k:','LineWidth',0.5) % line for start postion
%                 % %                 plot([amplicon_info(indx2(k),3) amplicon_info(indx2(k),3)],[-0.05 5.05],'k--.','LineWidth',0.5)% line fore end postion
%                 plot([amplicon_info(indx2(k),2) amplicon_info(indx2(k),2)],[-0.05 7.05],'k:','LineWidth',0.5) % line for start postion
%                 plot([amplicon_info(indx2(k),3) amplicon_info(indx2(k),3)],[-0.05 7.05],'k--.','LineWidth',0.5)% line fore end postion
%                 subplot(4,1,3)
%                 plot([amplicon_info(indx2(k),2) amplicon_info(indx2(k),2)],[0 1],'k:','LineWidth',0.5) % line for start postion
%                 plot([amplicon_info(indx2(k),3) amplicon_info(indx2(k),3)],[0 1],'k--.','LineWidth',0.5)% line fore end postion
%                 subplot(4,1,4)
%                 plot([amplicon_info(indx2(k),2) amplicon_info(indx2(k),2)],[-2.5 1.5],'k:','LineWidth',0.5) % line for start postion
%                 plot([amplicon_info(indx2(k),3) amplicon_info(indx2(k),3)],[-2.5 1.5],'k--.','LineWidth',0.5)% line fore end postion
%             end
%         end
%     end

    %save figure
%     figpath = [plotsdir '\Chr_' num2str(i) '_' barcode];
    figpath = [plotsdir '/Chr_' num2str(i) '_' barcode '.png'];
    %     eval(['print -djpeg -r600 ' figpath ])
    eval(['print -dpng -r200 ' figpath ])

end
