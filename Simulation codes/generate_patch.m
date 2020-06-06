% [p,Nactual,a] = generate_patch(N)
% p - data structure vector
% Nactual - size of p
% a - distance between sites (neighboring HA trimers)
% N - approximately how many particles do you want?

function [p,Nactual,a] = generate_patch(N)
       
    a = 1;
    etasafety = 1.1;   % 1.1 is safety factor

    diam = 2*sqrt(N/pi)*a;
    NX = ceil(diam/a * etasafety /2)*2;
    NY = ceil(diam/(sqrt(3)/2*a) * etasafety /4)*4;
    xvec = [-NX/2:NX/2] * a;
    yvec = [-NY/2:NY/2] * a*sqrt(3)/2;
    
    [xx,yy] = meshgrid(xvec,yvec);        % Create 2D matrix
    xx(2:2:end,:) = xx(2:2:end,:) + a/2;  % Shifted pixels to create hexagonal lattice
    
    idx = find((xx.^2 + yy.^2) <= (diam/2)^2);  % Find only points inside circle;
    xfinal = xx(idx); yfinal = yy(idx);
    Nactual = length(xfinal);

    % Plot Contact Patch
    figure; plot(xfinal,yfinal,'bx',xx(:),yy(:),'yo','MarkerSize',15); axis image;
    xlim(2*[-1.2 1.2]*diam*etasafety/2); ylim(2*[-1.2 1.2]*diam*etasafety/2);
    title(sprintf('Actual Patch Size = %d', Nactual));


    % Build data structure of circular network
    for kk = 1:length(xfinal)
        p.id(kk,1) = kk;                          % Site id
        p.xy(kk,:) = [xfinal(kk) yfinal(kk)];     % Site coordinates
        p.flipped(kk,1) = NaN;                    % Time flipped
        p.neighbors{kk} = [];                     % Not decided yet
    end
                
    for kk = 1:length(xfinal)
        for jj = 1:length(xfinal)
            if (jj ~= kk) && (sum((p.xy(jj,:) - p.xy(kk,:)).^2) < (1.01*a)^2)
                p.neighbors{kk}(end+1) = jj;
            end
        end
    end
    
end