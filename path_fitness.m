function [Fitness] = path_fitness(Gen, Dist)
%PATH_FITNESS Calculates fitness of path format phenotype
%   Gen matrix of individuals in path format
%   Dist matrix with precalculated distances

    dist_size = size(Dist);

    Fitness = Dist(sub2ind(dist_size, Gen(:,1), Gen(:,2)));
    for t = 3:size(Gen, 2)
        Fitness = Fitness + Dist(sub2ind(dist_size, Gen(:, t-1), Gen(:, t)));
    end

end

