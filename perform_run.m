function [ Path, BestFitness, BestFV, MeanFV, WorstFV ] = perform_run( ...
    Repre, ...
    ParSel, ...
    SurvSel, ...
    ...
    NCities, ...
    Distances, ...
    NInd, ...
    MaxGen, ...
    PXover, ...
    PMut, ...   
    Elite, ...
    StopPercentage)

%PERFORM_RUN Perform a run of GA (based on run_ga.m)
%   Input args:
%   * Repre - used encoding (adj or path)
%   * ParSel - parent selection method; might be in three forms:
%       [0, _selection pressure_, _ranking method_] - ranking method
%           selection pressure is a value in range (1; 2>
%           ranking method can be: 0 - linear; or 1 - non-linear.
%       [1, _tournament size_] - tournament method
%       [2] - proportional method
%   * SurvSel - survivor selection (insertion) method. Scalar of:
%       0 - elitism only, 1 - elitism and mu+lambda selection
%   * NCities - number of cities
%   * Distances - calculated distance matrix
%   * NInd - size of the population
%   * MaxGen - max number of generations
%   * PXover - probability of cross-over
%   * PMut - probability of mutation
%   * Elite - percentage of elite
%   * StopPercetage - percentage of equal fitness (stop criterium)
%
%   Output args:
%   * Path - the best path found
%   * BestFitness - fitness of the found path
%   * BestFV - matrix of how the best fitness evolved
%   * MeanFV - matrix of how the mean fitness evolved
%   * WorstFV - matrix of how the worst fitness evolved

    global TourSize;
    global Dist;
    Dist = Distances;

    Chrom = zeros(NInd,NCities);
    
    fitness_fun = 'tspfun';
    xover_fun = 'xalt_edges';
    mut_fun = 'inversion';
    
    switch Repre
        case 'adj'
            Repre = 1;
            
            for row = 1:NInd
                Chrom(row, :) = path2adj(randperm(NCities));
            end
            
        case 'path'
            Repre = 2;
            fitness_fun = 'path_fitness';
            xover_fun = 'perform_scx';
            mut_fun = 'reciprocal_exchange';
            
            for row = 1:NInd
                Chrom(row, :) = randperm(NCities);
            end
    end
    
    GGAP = 1 - Elite;
    
    BestFV = zeros(1, MaxGen);
    MeanFV = zeros(1, MaxGen);
    WorstFV = zeros(1, MaxGen);
    
    % number of individuals of equal fitness needed to stop
    stopN = ceil(StopPercentage * NInd);
    
    % evaluate initial population
    ObjV = feval(fitness_fun, Chrom, Dist);
    
    % generational loop
    gen = 0;
    while gen < MaxGen - 1
        sObjV = sort(ObjV);
        BestFV(gen + 1) = sObjV(1); %min(ObjV);
        BestFitness = BestFV(gen + 1);
        MeanFV(gen + 1) = mean(ObjV);
        WorstFV(gen + 1) = sObjV(end); %max(ObjV);
        
        for t = 1:size(ObjV, 1)
            if (ObjV(t) == BestFitness)
                break;
            end
        end
        
        Path = Chrom(t, :);

        if (sObjV(stopN) - sObjV(1) <= 1e-15)
            break;
        end

        switch ParSel(1)
            case 0
                % assign fitness values to entire population
                FitnV = ranking(ObjV, ParSel(2:end));
                
                % select individuals for breeding
                SelCh = select('sus', Chrom, FitnV, GGAP);

            case 1
                TourSize = ParSel(2);
                SelCh = select('tournament', Chrom, ObjV, GGAP);

            case 2
                SelCh = select('rws', Chrom, ObjV, GGAP);
                
        end

        % recombine individuals (crossover)
        SelCh = recombin(xover_fun, SelCh, PXover);
        SelCh = mutateTSP(mut_fun, SelCh, PMut, Repre);
        
        % evaluate offsprings, call objective function
        ObjVSel = feval(fitness_fun, SelCh, Dist);
        
        % reinsert offsprings into population
        switch SurvSel
            case 0 %elitism
                [Chrom, ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
            case 1 %mu+lambda with elitism
                [Chrom, ObjV]=perform_mu_lambda_survival(Chrom,SelCh,ObjV,ObjVSel);
        end

        % increment generation counter
        gen = gen + 1;
    end
    
    %If the run was finished because of Maxgen: the last generation is not 
    %even used (just thrown away), so there is actually Maxgen+1 with 
    %gen < Maxgen condition => use gen < Maxgen-1 instead!
    %Hence the last evaluation is here:
    if (gen == MaxGen-1)
        BestFV(gen + 1) = min(ObjV);
        BestFitness = BestFV(gen + 1);
        MeanFV(gen + 1) = mean(ObjV);
        WorstFV(gen + 1) = max(ObjV);

        for t = 1:size(ObjV, 1)
            if (ObjV(t) == BestFitness)
                break;
            end
        end

        Path = Chrom(t, :);
    end
    
    gen = gen + 1;
    
    BestFV = BestFV(1:gen);
    MeanFV = MeanFV(1:gen);
    WorstFV = WorstFV(1:gen);

end
