function SelV = tournament(FitV, NSel)
%TOURNAMENT Performs tournament selection
%   FitV vector with fitness values for the individuals
%   NSel number of individuals to select
%
%   Uses global variable TourSize - tournament size
%
%   Returns vector SelV with indices for selected individuals

    global TourSize;

    PoolSize = size(FitV, 1);
    
    SelV = zeros(NSel, 1);
    
    for i = 1:NSel
        
        RandSel = randperm(PoolSize, TourSize);
        RandSel = [RandSel' FitV(RandSel)];
    
        [~, midx] = min(RandSel(:,2));
        
        SelV(i) = RandSel(midx);
        
    end

end
