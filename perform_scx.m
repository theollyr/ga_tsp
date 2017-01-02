function [ NewChrom ] = perform_scx(OldChrom, Prob)
%PERFORM_SCX Performs sequential constructive crossover
%   OldChrom pool of individuals on which to perform xover
%   Prob probability of xover

    % using the global matrix with precalculated distances
    global Dist

    NewChrom = zeros(size(OldChrom));

    for t=2:size(OldChrom,1)
        if rand < Prob
            NewChrom(t-1,:) = seq_const_x([OldChrom(t-1,:); OldChrom(t,:)], Dist);
        else
            NewChrom(t-1,:) = OldChrom(t-1,:);
        end
    end

    NewChrom(end,:) = OldChrom(end,:);
end

