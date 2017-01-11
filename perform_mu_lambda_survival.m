function [NewChrom, NewObjV] = perform_mu_lambda_survival(Chrom, SelCh, ObjV, ObjVSel)
    %perform_mu_lambda_survival: Mu+Lambda survivor selection on top of
    %elitism!
    %   Chrom is array with parents
    %   SelCh is array with offsprings
    %   ObjV is the fitness of parents
    %   ObjVSel is the fitness of offsprings

    mu = size(Chrom,1);
    lambda = size(SelCh,1);

    %the pool of candidates, first add the offsprings
    poolChrom = SelCh;
    poolObjV = ObjVSel;

    %if mu > lambda, the (mu-lambda) number of elites must be preserved
    eliteChrom = [];
    eliteObjV = [];
    if(mu > lambda)
        [~,sortedIndex] = sort(ObjV);

        %first mu-lambda are kept
        eliteIndex = sortedIndex(1:(mu-lambda));
        eliteChrom = Chrom(eliteIndex,:);
        eliteObjV = ObjV(eliteIndex);

        %the others are part of the pool!
        othersIndex = sortedIndex((mu-lambda+1):end);
        poolChrom = [Chrom(othersIndex,:) ; poolChrom];
        poolObjV = [ObjV(othersIndex) ; poolObjV];
    else %in the event no elitism is used, so mu <= lambda
        poolChrom = [Chrom ; poolChrom];
        poolObjV = [ObjV ; poolObjV];
    end

    %select the best lambda number of candidates out of the pool
    [~, sortedIndex] = sort(poolObjV);
    bestlambdaIndex = sortedIndex(1:lambda);

    %add elites as well
    NewChrom = [eliteChrom ; poolChrom(bestlambdaIndex,:)];
    NewObjV = [eliteObjV ; poolObjV(bestlambdaIndex)];

end


