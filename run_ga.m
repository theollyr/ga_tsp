function run_ga(x, y, NIND, MAXGEN, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, CROSSOVER, LOCALLOOP, ah1, ah2, ah3)
% usage: run_ga(x, y,
%               NIND, MAXGEN, NVAR,
%               ELITIST, STOP_PERCENTAGE,
%               PR_CROSS, PR_MUT, CROSSOVER,
%               ah1, ah2, ah3)
%
%
% x, y: coordinates of the cities
% NIND: number of individuals
% MAXGEN: maximal number of generations
% ELITIST: percentage of elite population
% STOP_PERCENTAGE: percentage of equal fitness (stop criterium)
% PR_CROSS: probability for crossover
% PR_MUT: probability for mutation
% CROSSOVER: the crossover operator
% calculate distance matrix between each pair of cities
% ah1, ah2, ah3: axes handles to visualise tsp
% {NIND MAXGEN NVAR ELITIST STOP_PERCENTAGE PR_CROSS PR_MUT CROSSOVER LOCALLOOP}


        global TourSize;
        TourSize = 10;

        parent_selection = 'ranking';
        fitness_fun = 'tspfun';
        survivor_selection = 'elitism';
        mut_fun = 'inversion';
        ind_repre = 2;

        if ind_repre == 2
            fitness_fun = 'path_fitness';
            CROSSOVER = 'perform_scx';
        end

        GGAP = 1 - ELITIST;
        mean_fits=zeros(1,MAXGEN+1);
        worst=zeros(1,MAXGEN+1);

        global Dist
        Dist=zeros(NVAR,NVAR);

        for i=1:size(x,1)
            for j=1:size(y,1)
                Dist(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            end
        end
        % initialize population
        Chrom=zeros(NIND,NVAR);

        if ind_repre == 1
            for row=1:NIND
                Chrom(row,:)=path2adj(randperm(NVAR));
            end
        else
            for row=1:NIND
                Chrom(row,:)=randperm(NVAR);
            end
        end

        gen=0;
        % number of individuals of equal fitness needed to stop
        stopN=ceil(STOP_PERCENTAGE*NIND);
        % evaluate initial population
        ObjV = feval(fitness_fun, Chrom, Dist);
        best=zeros(1,MAXGEN);
        % generational loop
        while gen<MAXGEN
            sObjV=sort(ObjV);
          	best(gen+1)=min(ObjV);
        	minimum=best(gen+1);
            mean_fits(gen+1)=mean(ObjV);
            worst(gen+1)=max(ObjV);
            for t=1:size(ObjV,1)
                if (ObjV(t)==minimum)
                    break;
                end
            end

            if ind_repre == 1
                visualizeTSP(x,y,adj2path(Chrom(t,:)), minimum, ah1, gen, best, mean_fits, worst, ah2, ObjV, NIND, ah3);
            else
                visualizeTSP(x,y,Chrom(t,:), minimum, ah1, gen, best, mean_fits, worst, ah2, ObjV, NIND, ah3);
            end

            if (sObjV(stopN)-sObjV(1) <= 1e-15)
                  break;
            end

            switch parent_selection
                case 'ranking'
                    % assign fitness values to entire population
                    FitnV=ranking(ObjV);
                    %select individuals for breeding
                    SelCh=select('sus', Chrom, FitnV, GGAP);

                case 'tournament'
                    SelCh = select('tournament', Chrom, ObjV, GGAP);

                case 'proportional'
                    SelCh = select('rws', Chrom, ObjV, GGAP);
            end

        	%recombine individuals (crossover)
            SelCh = recombin(CROSSOVER,SelCh,PR_CROSS);
            %SelCh=mutateTSP('inversion',SelCh,PR_MUT,ind_repre);
            SelCh=mutateTSP(mut_fun,SelCh,PR_MUT,ind_repre);
            %evaluate offsprings, call objective function
        	ObjVSel = feval(fitness_fun, SelCh, Dist);
            
            %reinsert offsprings into population
            switch survivor_selection
                case 'elitism'
                   [Chrom, ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
                case 'mu_lambda'
                   [Chrom, ObjV]=perform_mu_lambda_survival(Chrom,SelCh,ObjV,ObjVSel);
            end
            

            Chrom = tsp_ImprovePopulation(NIND, NVAR, Chrom,LOCALLOOP,Dist);
        	%increment generation counter
        	gen=gen+1;
        end
end