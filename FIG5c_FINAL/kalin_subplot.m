% subplot with smaller margins
function ax = kalin_subplot(M,N, k, space, label_offset)

    if (nargin < 4) || isempty(space)
        spaceX = 0.1;
        spaceY = 0.1;
    else
        if length(space) == 2
            spaceX = space(1);
            spaceY = space(2);
        else
            spaceX = space(1);
            spaceY = spaceX;
        end
    end
    
    if (nargin < 5) || isempty(label_offset)
        label_offset = [0.05 0.05];
    end
    
    
    [i, j] = ind2sub([M N], k);
    Lx = (1-label_offset(1))/N; Ly = (1-label_offset(2))/M;
    ax = axes('position', [label_offset(1)+spaceX*Lx/2 + (j-1)*Lx, (1-label_offset(2)) - i*Ly + spaceY*Ly/2, Lx*(1-spaceX), Ly*(1-spaceY)]);
    
end