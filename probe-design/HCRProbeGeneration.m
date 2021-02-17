% clear 
% close all 

numProbes = 6;
%probes = {};
probes = strings;
Files = dir(strcat(directory, '/HomologyRegionsHCR/*.fa'));

check_probe_number = false;
if probes_per_gene_file
    probes_per_gene = readtable(probes_per_gene_file, 'ReadVariableNames', false); % table with column names Var1 Var2 of gene name and probe number
    check_probe_number = true;
end

% open file for writing bad probes
fileBadGenes = fopen('genes_with_no_probes.txt','wt');

% Files
for j=1:length(Files)
    Files(j).name

    % get name
    name_arr = split(Files(j).name, '-');
    gene_name = join(name_arr(1:length(name_arr)-1),'-');

    numProbes_for_gene = numProbes;
    if check_probe_number
        if ~isempty(probes_per_gene{strcmp(probes_per_gene.Var1, gene_name), 'Var2'})
            numProbes_for_gene = probes_per_gene{strcmp(probes_per_gene.Var1, gene_name), 'Var2'}(1)
        end
    end

   % [Header, Sequence] = fastaread(Files(j).name);
    [Header, Sequence] = fastaread(strcat(directory, '/HomologyRegionsHCR/', Files(j).name));
    
    % check if empty
    if isa(Header, 'char') & string(Header) == 'WARNING: NO VALID PROBES'
        fprintf(fileBadGenes, '%s\n', Files(j).name);
        continue
    end

    if isa(Sequence, 'char')
        Sequence = string(Sequence); 
    end

    k = length(Sequence);
    k
    % if ~isempty(strfind(Files(j).name, 'intron'))
    %     k = min(k, 30) % on the off chance we got a bazillion intron probes, don't look for ones too far away from the TSS
    % end
    idx = floor(linspace(1,k,numProbes_for_gene));
    idx = unique(idx);
    Sequence(idx)

    %probes{j,1} = string(Files(j).name);
    probes(j,1,1:length(idx)) = strcat(regexprep(Files(j).name,'.fa','','ignorecase'),'-',string(idx));
    probes(j,2,1:length(idx)) = string(Sequence(idx));
    for m=1:length(idx)
        p = char(probes(j,2,m));
        p;
        probes(j,3,m) = strcat(string(p(1:25)),'TA','GAAGAGTCTTCCTTTACG'); % 3 prime
        probes(j,4,m) = strcat('/5Phos/GGAGGGCAGCAAACGG','AA',string(p(28:end)),'NNNNNNNNNN','CTCGACCGTTAGCAAAGCTC'); % 5 prime
    end
end
probes1 = cat(1,probes(:,:));
probes1(ismissing(probes1)) = string(); %"";
size(probes1)
probes1
filePh = fopen('probes.csv','w');
fprintf(filePh,strcat(repmat('%s,',1,length(Files)-1), '%s\n'),probes1{:});
fclose(filePh);
fclose(fileBadGenes);
% strcat(repmat('%s,',1,4*numProbes-1), '%s\n')
% writetable(probes1, 'probes.csv')

% xlswrite('probes.xlsx',probes1);
% string()