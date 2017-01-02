function offspring = seq_const_x(Chroms, Dist)
%SEQ_CONST_X Sequential constructive crossover operator
%   Chroms array with two parents
%   Dist matrix with precalculated distances
%   Outputs a new offspring

    OffLen = size(Chroms, 2);

    curr_node = 1;
    offspring = zeros(1, OffLen);
    offspring(1) = 1;

    Par1 = Chroms(1,:);
    Par2 = Chroms(2,:);

    Pool = 2:OffLen;

    for t=Pool
        en1 = find_legitimate_node(Par1);
        en2 = find_legitimate_node(Par2);

        if isnan(en1)
            LegitNodes = Pool(~ismember(Pool, offspring));
            en1 = LegitNodes(1);
        end

        if isnan(en2)
            LegitNodes = Pool(~ismember(Pool, offspring));
            en2 = LegitNodes(1);
        end

        if Dist(curr_node,en1) < Dist(curr_node,en2)
            offspring(t) = en1;
            curr_node = en1;
        else
            offspring(t) = en2;
            curr_node = en2;
        end
    end

    function LegitNode = find_legitimate_node(Chrom)
        Rest = Chrom((find(Chrom == curr_node) + 1):end);
        Rest = Rest(~ismember(Rest, offspring));

        if isempty(Rest)
            LegitNode = NaN;
        else
            LegitNode = Rest(1);
        end
    end
end

