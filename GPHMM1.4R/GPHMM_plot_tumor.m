function GPHMM_plot_tumor(Version,Outputsource,barcode,chr,data_pos_sep,w,lrr_shift,gc_coef,data_genotype_sep2,data_bafA_sep,score)
%this function is used to plot 22 figures given the SNP array data as well
%as the results from NormCNV, all these figures are stored in a specified
%directory
global data_lrr_sep
global data_baf_sep
global data_gc_sep

% load matlab2
%---for cn up to 7---
cn = [0.01 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 ];
Bcn= [0 1 1 2 2 3 3 2 4 4 3 5 5 4 3 6 6 5 4 7 ];
lrr= [0 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 ];
loh= [1 1 0 1 0 1 0 0 1 0 0 1 0 0 0 1 0 0 0 1 ];
cn_unique=unique(cn);
loh_unique=unique(loh);
lrr_unique=unique(lrr);
Colors=0.65.*jet(length(cn_unique));
Bcn= [0.005 1 1 2 2 3 3 2 4 4 3 5 5 4 3 6 6 5 4 7 ];
AI_mapping = [1 1 3 2 3 2 3 3 2 3 3 2 3 3 3 2 3 3 3 2 1];
Colors_lrr=[%153 192 13;...
    %20 123 57;...
    15 111 50;...
    33 172 71;...
    164 164 164;...
    %84 115 183;...
    226 29 4]./255;
Colors_baf=[164 164 164;...
    %20 123 57;...
    84 115 183;...
    ]./255;

%--------------- plot figures ---------------------%
plotsdir = [Outputsource '/' barcode ];
mkdir(plotsdir);
for j=1:length(data_genotype_sep2)
    clf reset
    set(gcf,'visible','off');
    
    
    marker_size = 1.5;
    data_baf=data_baf_sep{j};
    data_baf(data_bafA_sep{j})=1-data_baf_sep{j}(data_bafA_sep{j});
    %---plot BAF---
    subplot(3,1,1)
    hold on
    %plot(data_pos_sep{j},data_baf_sep{j},'b.', 'MarkerSize',marker_size)
    for i=1:length(loh_unique)
        tv=ismember(data_genotype_sep2{j},find(loh==loh_unique(i)));
        plot(data_pos_sep{j}(tv), data_baf(tv),'.', 'MarkerSize',marker_size*2,'Color',Colors_baf(i,:));
    end
    %plot(data_pos_sep{j},data_baf_sep{j},'b.', 'MarkerSize',marker_size)
    %plot expected BAF mean values
    
    temp=((1-w)+w.*Bcn(data_genotype_sep2{j}))./((1-w).*2+w.*cn(data_genotype_sep2{j}));
    line_style='r.';
    plot(data_pos_sep{j}, temp,line_style, 'MarkerSize',marker_size);
    plot(data_pos_sep{j}, 1-temp,line_style, 'MarkerSize',marker_size);
    
    ylabel('BAF');
    axis([-inf inf -0.08 1.08]);
    YLabelPos=[ 0 0.25 0.5 0.75 1];
    YLabel=[{'0.0'},{'0.25'},{'0.5'},{'0.75'},{'1.0'}];
    set(gca,'YTick',YLabelPos);
    set(gca,'YTickLabel',YLabel,'YGrid','on','FontSize',6,'Box','on');
    set(gca,'XTick',[]);
    
    temp_barcode=barcode;
    tmp = strfind(temp_barcode,'_');
    temp_barcode(tmp) = '-';
    title (['GPHMM' Version ': ' temp_barcode ', Tumor Content: ' num2str(int32(w*100)) '%'],'FontSize',8);
    clear tmp temp_barcode
    %axes('Position',[0.1 0.55 0.8 0.8]);
    
    %---plot LRR---
    subplot(3,1,2);
    hold on
    %plot(data_pos_sep{j},data_baf_sep{j},'b.', 'MarkerSize',marker_size)
    for i=1:length(lrr_unique)
        tv=ismember(data_genotype_sep2{j},find(lrr==lrr_unique(i)));
        plot(data_pos_sep{j}(tv), data_lrr_sep{j}(tv),'.', 'MarkerSize',marker_size*2,'Color',Colors_lrr(i,:));
    end
    %plot(data_pos_sep{j},data_lrr_sep{j},'b.', 'MarkerSize',marker_size);
    
    %plot expected LRR mean values
    temp = 2*log10(((1-w)*2+w.*cn(data_genotype_sep2{j}))./2)+lrr_shift+gc_coef.*data_gc_sep{j};
    line_style='k.';
    plot(data_pos_sep{j}, temp,line_style, 'MarkerSize',marker_size*2);
    clear temp
    
    ylabel('LRR');
    axis ([-Inf Inf -2.2 1.5]);
    YLabelPos=[-1:0.5:1];
    set(gca,'YTick',YLabelPos);
    set(gca,'YGrid','on','FontSize',6,'Box','on');
    set(gca,'XTick',[]);
    marker_size = 5;
    LineWid=3;
    %---plot CN---
    subplot(3,1,3)
    hold on
    %axes('Position',[0.1 0.1 0.8 0.3]);
    color_cn=jet(6);
    plot([data_pos_sep{j}(1) data_pos_sep{j}(end)], [-1 -1],'k-');
    tv=score{j}<0.05;
    plot(data_pos_sep{j}(tv), 7*(1-score{j}(tv)),'b.', 'MarkerSize',1);
    %%%% ---- Tumor CN plot ------
    startindx=1;
    while (startindx<length(data_genotype_sep2{j}))
        endindx=startindx+1;
        while(data_genotype_sep2{j}(startindx)==data_genotype_sep2{j}(endindx) && endindx<length(data_genotype_sep2{j}))
            endindx=endindx+1;
        end
        temp=cn(data_genotype_sep2{j}(startindx));
        line_style='-';
        plot([data_pos_sep{j}(startindx) data_pos_sep{j}(endindx-1)], [temp temp],line_style, 'LineWidth',LineWid/2,'Color',[1 0 0]);
        clear temp
        startindx=endindx;
        if startindx>=length(data_genotype_sep2{j})
            break;
        end
    end
    
    set(gca,'YGrid','on','Box','on')
    %     axis ([-Inf Inf -0.05 7.05])
    axis ([-Inf Inf -0.5 7.5])
    set(gca,'YTick',[0:2:6])
    ylabel('CN','FontSize',8);
    xlabel(['Chromosome ' num2str(chr(j))],'FontSize',8);
    set(gca,'YTickLabel',{'0','2','4','6'},'FontSize',5);
    
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 450 260]);
 
    figpath = [plotsdir '/Chr_' num2str(chr(j)) '_' barcode '.png'];
    eval(['print -dpng -r600 ' figpath ])
end
end
