% Generates probes for multiplex experiment, does one gene at a time
% Comes from fasta file, if exon will be one entry, if intron will be many
% Doesn't matter tho, since we will iterate through all things in the fasta file!

% clear all;


probelength = 52;
spacerlength = 50;

GC_cutoff = [40 65];

% stuff that comes in from command line
% fasta_name (path to fasta)
% ref_seq (id for writing)
% genename (for checking hits against)
% fasta comes in from the command line 

[Header, Sequence] = fastaread(fasta_name);

length(Sequence);

if isa(Sequence, 'char')
    Sequence = string(Sequence); 
end

length(Sequence);

% class(Sequence)

List = {};
probenum = 1;

sequence_features = 0;
blast_failure = 0;
correct = 0;
total=0;

for y=1:length(Sequence)
    y
    S = char(Sequence(y));
    S

    tail = probelength + spacerlength + 1; %index to where the tail of the probe design is

    while tail < length(S)                
        %blast against excluded sequences
        TempSeq = S((tail - probelength - spacerlength): tail-spacerlength -1);
        total = total + 1;

        % disp(TempSeq);
        fastawrite(strcat('fastas/', ref_seq, '-', intron, '-blast-temp.fa'), TempSeq); 
  
        hits = false;
        blasthits=false;
        %check for repeats of same base
        if seqwordcount(TempSeq,'CCCCC') > 0
            hits = true;
        end
        if seqwordcount(TempSeq,'GGGGG') > 0
            hits = true;
        end
        if seqwordcount(TempSeq,'AAAAG') > 0
            hits = true;
        end
        if seqwordcount(TempSeq,'TTTTT') > 0
            hits = true;
        end
        %BSRDI
        %restriction site exclusion
%         if seqwordcount(TempSeq,'GCAATG') > 0
%           hits = true;
%         end
        
        %BSAI
        %restriction site exclusion
        if seqwordcount(TempSeq,'GGTCTC') > 0
          hits = true;
        end
        if seqwordcount(TempSeq,'GAGACC') > 0
          hits = true;
        end
        
        %check GC
        seqProperties = oligoprop(TempSeq, 'HPBase', 7, 'HPLoop', 6, 'DimerLength', 10);
        if seqProperties.GC < GC_cutoff(1)
            hits = true;
        end
        if seqProperties.GC > GC_cutoff(2)
            hits = true;
        end
        %check for hairpins or dimers
        if ~isempty(seqProperties.Hairpins) || ~isempty(seqProperties.Dimers)
            hits = true;
        end

        % check lowercase for repeats! can include 0-5 lowercase bases
        if numel(TempSeq) - nnz(TempSeq == upper(TempSeq)) > 5
            hits = true;
            TempSeq;
        end

        if length(TempSeq) < probelength
            TempSeq
            hits = true;
        end
        
        %only blast if sequence properties are good
        if hits == false 
            blastout = blastlocal('InputQuery',strcat('fastas/', ref_seq, '-', intron, '-blast-temp.fa'),'Program','blastn','DataBase', blastdb, 'blastargs', '-S 1 -F F');
            % for i = 1:length(blastout.Hits)
            % for i = 1:numel(blastout.Hits)
            for i = 1:length(blastout.Hits)
                if isempty(strfind(blastout.Hits(i).Name, genename)) %ignore same gene
                    if isempty(strfind(blastout.Hits(i).Name, 'PREDICTED')) && isempty(strfind(blastout.Hits(i).Name, 'predicted')) %ignore predicted
                        if isempty(strfind(blastout.Hits(i).Name, 'RIKEN'))
                            for j = 1:max(size(blastout.Hits(i).HSPs))
                                if ~isempty(strfind(blastout.Hits(i).HSPs(j).Strand, 'Plus/Plus')) %Look for only Plus/Plus
                                    if blastout.Hits(i).HSPs(j).Identities.Match >25 %25 bases of identity or greater
                                        blastout.Hits(i).HSPs(j).Identities.Match;
                                        if blastout.Hits(i).HSPs(j).QueryIndices(1) < 12 &&  blastout.Hits(i).HSPs(j).QueryIndices(2)>20 %well centered
                                            alignment = blastout.Hits(i).HSPs(j).Alignment;
                                            blastout.Hits(i).Name
                                            blasthits = true;
                                            blastout.Hits(i).HSPs(j)
                                            TempSeq;
                                        end
                                    end
                                end
                            end
                        end
                    else
                        blastout.Hits(i);
                    end
                end
            end
            blastres{probenum} = blastout;
        end
            
        if blasthits==true
            tail=tail+8;    
            blast_failure = blast_failure + 1;    
        elseif hits == true
            tail = tail + 5; %shift window by 2 if probes not found
            sequence_features = sequence_features + 1;
        else
            probeindx{probenum} = tail - probelength - spacerlength;
            tail = tail + probelength + spacerlength;
            List{probenum} = TempSeq;
            probenum = probenum + 1;
            % disp(TempSeq);
            TempSeq;
            correct = correct + 1;

        end
        delete(strcat('fastas/', ref_seq, '-', intron, '-blast-temp.fa'));
    end
end

size(List); 

sequence_features
blast_failure
correct
total

% no valid probes
if correct == 0
    fasta.Sequence = 'N';
    fasta.Header = 'WARNING: NO VALID PROBES';
    fastawrite(outfasta, fasta)
    % fclose(fopen(outfasta, 'w'));
    % also write to some universal empty file thing
end


for ii = 1:max(size(List))
    RevList{ii} = seqrcomplement(List{ii});
    % oldfolder = cd('HomologyRegionsHCR');
    fasta.Sequence = RevList{ii};
    fasta.Header = ['probe' num2str(ii) ' index ' num2str(probeindx{ii})];
    % fastawrite(strcat('HomologyRegionsHCR/', ref_seq, '-', intron.fa'), fasta);
    fastawrite(outfasta, fasta)
    % cd(oldfolder);
end

clear probeindx;
clear List;
