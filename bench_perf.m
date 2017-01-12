function bench_perf( ...
    OutDir, ...
    Datasets, ...
    NRuns, ...
    ...
    Repres, ...
    ParSels, ...
    SurvSels, ...
    ...
    NInds, ...
    MaxGens, ...
    PXovers, ...
    PMuts, ...
    Elites)
%BENCH_PERF Benchmark performance
%   Runs benchmark for multiple datasets with various algorithm setups
%   and parameters.
%
%   Input args:
%   * OutDir - directory to save the bench results
%   * Datasets - filepaths to datasets (array of strings)
%   * NRuns - number of runs per one setup
%   * Repres - array of strings with encodings to try out
%   * ParSels - matrix with parent selection setups to use (setup per row)
%   * SurvSels - array with survival selection method
%   * NInds - array with generation sizes
%   * MaxGens - array with maximal number of generations to use
%   * PXovers - array with values of xover probability
%   * PMuts - array with values of mutation probability
%   * Elites - array with values of elite percentage
%
%   The function creates the following file structure
%   - OutDir/
%       - bench_report.txt - file with full report of the benchmark
%       - Dataset1
%           - map.pdf - map of the given dataset
%           - algo0
%               - run0
%                   - path.pdf - the best path found
%                   - fitness.pdf - a graph showing evolution of fitness
%                     values
%                   - result.txt - summary for the given run
%               - run1 ...
%           - algo1
%               -run0 ...
%           - algo2 ...
%       - Dataset2 ...

    tic

    report_dir = OutDir;
    mkdir(report_dir)
    report_file = fopen(fullfile(report_dir, 'bench_report.txt'), 'w');

    fprintf(report_file, 'Benchmark run at %s\n\n', datestr(now));
    fprintf(report_file, 'Using these datasets: %s\n', strjoin(Datasets, ', '));
    fprintf(report_file, 'Number of runs per setup: %d\n', NRuns);
    fprintf(report_file, '\nAlgorithms:\n');

    algon = 0;
    for repre = Repres
        for parSelIdx = 1:size(ParSels, 1)
            parSel = ParSels(parSelIdx, :);

            switch parSel(1)
                case 0
                    if parSel(3) == 0
                        rankMethod = 'linear';
                    else
                        rankMethod = 'non-linear';
                    end

                    parDesc = sprintf('parent selection method: ranking (selection pressure %.4f, %s)', parSel(2), rankMethod);

                case 1
                    parDesc = sprintf('parent selection method: tournament (size %d)', parSel(2));

                case 2
                    parDesc = 'parent selection method: proportional';
            end

            for survSel = SurvSels

                if survSel == 0
                    survKind = 'elitism';
                else
                    survKind = 'mu+lambda with elitism';
                end

                fprintf(report_file, '%d - %s encoding, %s, survivor selection type: %s\n', ...
                    algon, repre, parDesc, survKind);

                algon = algon + 1;
            end
        end
    end

    fprintf(report_file, '\nRuns:\n');

    for datasetpath = Datasets
        dataset = load(datasetpath);

        [~, datasetname, ~] = fileparts(char(datasetpath));

        fprintf(report_file, '\nDataset %s:\n', datasetname);

        dataset_dir = fullfile(report_dir, datasetname);
        mkdir(dataset_dir);

        mx = max([dataset(:, 1); dataset(:, 2)]);
        x = dataset(:, 1) / mx;
        y = dataset(:, 2) / mx;

        figure('visible', 'off');
        plot(x, y, 'ko');
        saveas(gcf, fullfile(dataset_dir, 'map.pdf'), 'pdf');
        close;

        NCities = size(x, 1);

        Distances = zeros(NCities);

        for i=1:NCities
            for j=1:NCities
                Distances(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            end
        end

        algon = 0;
        for repre = Repres
            for parSelIdx = 1:size(ParSels, 1)
                parSel = ParSels(parSelIdx, :);

                for survSel = SurvSels

                    fprintf(report_file, '\nAlgorithm %d runs:\n', algon);
                    algo_name = sprintf('algo%d', algon);
                    algo_dir = fullfile(dataset_dir, algo_name);
                    algon = algon + 1;

                    mkdir(algo_dir);

                    AlgoStatF = [0.0 Inf Inf];
                    AlgoStatN = zeros(1, 3);

                    runn = 0;
                    for NInd = NInds
                        for MaxGen = MaxGens
                            for PXover = PXovers
                                for PMut = PMuts
                                    for Elite = Elites

                                        BestBestFitness = Inf;
                                        WorstBestFitness = 0.0;
                                        BestFitnesses = 0.0;

                                        run_name = sprintf('run%d', runn);
                                        run_dir = fullfile(algo_dir, run_name);
                                        mkdir(run_dir);

                                        for r = 1:NRuns
                                            [PathTmp, BestFTmp, BestFVTmp, MeanFVTmp, WorstFVTmp] = perform_run(...
                                                repre, parSel, survSel, ...
                                                NCities, Distances, NInd, MaxGen, PXover, PMut, Elite, .95);

                                            BestFitnesses = BestFitnesses + BestFTmp;

                                            if BestFTmp < BestBestFitness
                                                Path = PathTmp;
                                                BestFV = BestFVTmp;
                                                MeanFV = MeanFVTmp;
                                                WorstFV = WorstFVTmp;

                                                BestBestFitness = BestFTmp;
                                            end

                                            if BestFTmp > WorstBestFitness
                                                WorstBestFitness = BestFTmp;
                                            end
                                        end

                                        MeanBestFitness = BestFitnesses / NRuns;

                                        if WorstBestFitness > AlgoStatF(1)
                                            AlgoStatF(1) = WorstBestFitness;
                                            AlgoStatN(1) = runn;
                                        end

                                        if MeanBestFitness < AlgoStatF(2)
                                            AlgoStatF(2) = MeanBestFitness;
                                            AlgoStatN(2) = runn;
                                        end

                                        if BestBestFitness < AlgoStatF(3)
                                            AlgoStatF(3) = BestBestFitness;
                                            AlgoStatN(3) = runn;
                                        end

                                        fprintf(report_file, ['%3d - NInd: %3d, MaxGen: %3d, PXover: %.2f, PMut: %.2f, Elite: %.2f; ' ...
                                            'BestBF: %.4f, MeanBF: %.4f, WorstBF: %.4f\n'], ...
                                            runn, NInd, MaxGen, PXover, PMut, Elite, BestBestFitness, MeanBestFitness, WorstBestFitness);

                                        runn = runn + 1;

                                        result_file = fullfile(run_dir, 'result.txt');
                                        path_file = fullfile(run_dir, 'path.pdf');
                                        fitness_file = fullfile(run_dir, 'fitness.pdf');

                                        result_file = fopen(result_file, 'w');

                                        fprintf(result_file, 'Parameters:\n\n');
                                        fprintf(result_file, 'Generation size: %d\n', NInd);
                                        fprintf(result_file, 'Maximum generations: %d\n', MaxGen);
                                        fprintf(result_file, 'Xover probability: %.2f\n', PXover);
                                        fprintf(result_file, 'Mutation probability: %.2f\n', PMut);
                                        fprintf(result_file, 'Elite %%: %.2f\n', Elite);

                                        fprintf(result_file, '\nResults:\n\n');
                                        fprintf(result_file, 'Best BestFitness found: %.4f\n', BestBestFitness);
                                        fprintf(result_file, 'Mean BestFitness: %.4f\n', MeanBestFitness);
                                        fprintf(result_file, 'Worst BestFitness found: %.4f\n', WorstBestFitness);

                                        fclose(result_file);

                                        figure('visible', 'off');
                                        plot(x(Path), y(Path), 'ko-', 'MarkerFaceColor', 'Black');
                                        hold on;
                                        plot([x(Path(end)) x(Path(1))], [y(Path(end)) y(Path(1))], 'ko-');
                                        hold off;
                                        saveas(gcf, path_file, 'pdf');
                                        close;

                                        figure('visible', 'off');
                                        plot(BestFV, 'r');
                                        hold on;
                                        plot(MeanFV, 'b');
                                        plot(WorstFV, 'g');
                                        hold off;
                                        xlabel('Generation');
                                        ylabel('Distance (Min. - Mean - Max.)');
                                        saveas(gcf, fitness_file, 'pdf');
                                        close;

                                    end
                                end
                            end
                        end
                    end

                    best_run_file = fopen(fullfile(algo_dir, 'best_run.txt'), 'w');
                    fprintf(best_run_file, 'Worst best fitness in run: %d\n', AlgoStatN(1));
                    fprintf(best_run_file, 'Best MBF in run: %d\n', AlgoStatN(2));
                    fprintf(best_run_file, 'Best best fitness in run: %d\n', AlgoStatN(3));
                    fclose(best_run_file);

                    fprintf(report_file, 'Worst best fitness in run: %d; ', AlgoStatN(1));
                    fprintf(report_file, 'Best MBF in run: %d; ', AlgoStatN(2));
                    fprintf(report_file, 'Best best fitness in run: %d\n', AlgoStatN(3));
                end
            end
        end
    end

    elapsedTime = toc;
    fprintf(report_file, '\n\nBenchmark elapsed time: %.4f seconds\n', elapsedTime);

    fclose(report_file);
end

