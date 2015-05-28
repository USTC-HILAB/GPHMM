function GPHMM_result_summary(path_results,path_summary)

fid = fopen(path_summary,'w');
if fid == -1
    error('Can not open summary file!');
end
fprintf(fid,'ID\tLRRshift\tGCcoef\tP_Normal\tSD_L\tSD_B\tP_Alt\tDNA_I\tEntropy\tP_C1\tP_C2\tP_Alt_C1\tP_Alt_C2\tP_Alt_C1&2\n');

datafilelist = dir(path_results);

for i=3:length(datafilelist)
    if isempty(strfind(datafilelist(i).name,'.results')) 
        continue;
    end
    fid2 = fopen([path_results '\' datafilelist(i).name],'r');
    if fid2 == -1
        error(['Can not read: ' datafilelist(i).name]);
    end
    
    %get estimated global parameters from the first row of the result file
    ID = datafilelist(i).name(1:end-8);
    w = [];
    while 1
        tline = fgetl(fid2);
        if ~isempty(strfind(tline,'StartPos')),break,end
        %LRR baseline shift
        result1 = regexp(tline,'LRR correction factor:\s*(\S+)','tokens','once');
        result2 = regexp(tline,'LRR baseline shift:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            lrr_shift = str2double(result1{1});
        elseif ~isempty(result2)
            lrr_shift = str2double(result2{1});
        end
        %GC coefficient
        result1 = regexp(tline,'GC coefficient:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            gc_coef = str2double(result1{1});
        end
        %normal cell proportion
        result1 = regexp(tline,'Proportion of normal cells in the sample:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            normal_p = str2double(result1{1});
        end
        %LRR SD
        result1 = regexp(tline,'Standard deviation of LRR signal:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            LRR_sd = str2double(result1{1});
        end
        %BAF SD
        result1 = regexp(tline,'Standard deviation of BAF signal:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            BAF_sd = str2double(result1{1});
        end
        %porportion of abnormal regions 
        result1 = regexp(tline,'Proportion of all abnormal chromosomal regions:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            abnormal_p = str2double(result1{1});
        end
        %overall DNA index
        result1 = regexp(tline,'Estimated average cancer DNA index:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            DNA_i = str2double(result1{1});
        end
        %Shannon's entropy of the prediction results
        result1 = regexp(tline,'Shannon''s entropy of the prediction results:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            entropy = str2double(result1{1});
        end
        %w1 and w2
        result1 = regexp(tline,'Proportion in the sample:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            w = [w str2double(result1{1})];
        end
        %porportion of C1-specific abnormal regions 
        result1 = regexp(tline,'Cancer clone 1 only:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            p1 = str2double(result1{1});
        end
        %porportion of C1-specific abnormal regions
        result1 = regexp(tline,'Cancer clone 2 only:\s*(\S+)','tokens','once');
        if ~isempty(result1)
            p2 = str2double(result1{1});
        end
    end
    fclose(fid2);
    
    %report errors if these values are not parsed successfully
    if isempty(lrr_shift)
        error(['Can not read estimated LRR shift from ',resultsfile]);
    end
    if isempty(gc_coef)
        error(['Can not read estimated GC coefficient from ',resultsfile]);
    end
    if isempty(normal_p)
        error(['Can not read estimated porportion of normal cells from ',resultsfile]);
    end
    if isempty(LRR_sd)
        error(['Can not read estimated standard deviation of LRR signals from ',resultsfile]);
    end
    if isempty(BAF_sd)
        error(['Can not read estimated standard deviation of BAF signals from ',resultsfile]);
    end
    if isempty(abnormal_p)
        error(['Can not read estimated porportion of abnormal regions from ',resultsfile]);
    end
    if isempty(DNA_i)
        error(['Can not read estimated DNA index from ',resultsfile]);
    end
    if isempty(entropy)
        error(['Can not read estimated Shannon''s entropy from ',resultsfile]);
    end
    if isempty(p1)
        error(['Can not read estimated proportion of C1-specific abnormal regions from ',resultsfile]);
    end
    if isempty(p2)
        error(['Can not read estimated proportion of C1-specific abnormal regions from ',resultsfile]);
    end
    if length(w)<2
        error(['Can not read estimated tumor proportion from ',resultsfile]);
    end

    if w(2) == 0 %single clone
        w(1) = 0;
        p1 = 0;
    end
    %output results to summary file
    fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',ID,lrr_shift,gc_coef,normal_p,LRR_sd,BAF_sd,abnormal_p,DNA_i,entropy,w(1),w(2),p1,p2,p1+p2);
end
fclose(fid);


