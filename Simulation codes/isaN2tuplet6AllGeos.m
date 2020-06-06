function [TF] = isaN2tuplet6AllGeos(p, n, N2, a)
%isaN2tuplet determines whether nth element of p completed a set of active N2 neighbors that flipped
    
TF = 0;
[ffn] = findFlippedNeighbors(p,n,n); % will go through the neighbors of p.id(n) and will include their indices within ffn if their p.flipped < p.flipped(n)
Lffn = length(ffn);

switch N2
    case 2
        if Lffn>0
        TF = 1;
        end
    case 3
        if Lffn>3
            TF = 1;
        elseif Lffn==2 || Lffn==3   % they will form a triplet if any two of the neighbors are connected
            fl = 0;
            for ii = 1 : (Lffn-1)
                for jj = (ii+1): Lffn
                    if (sum((p.xy(ffn(ii),:) - p.xy(ffn(jj),:)).^2) < (1.01*a)^2)
                    fl = 1;
                    break
                    end
                end
                if fl == 1
                TF = 1;
                break
                end
            end          
        end
    case 4
        if Lffn>4
            TF = 1;
        elseif Lffn == 2     % they will form a quadruplet if they have a shared neighbor
            neighbors1=findFlippedNeighbors(p,ffn(1),n);
            neighbors2=findFlippedNeighbors(p,ffn(2),n);
            if ~isempty(intersect(neighbors1,neighbors2));
                TF = 1;
            end
        elseif Lffn == 3    % they will form a quadrluplet if 1) three of them connected or 2) if two of them connected and they have an additional shared neighbor
            fl = 0;
            for ii = 1 : (Lffn-1)
                for jj = (ii+1): Lffn
                    if (sum((p.xy(ffn(ii),:) - p.xy(ffn(jj),:)).^2) < (1.01*a)^2)
                        fl = fl+1;
                        if fl == 1
                            idx=[ffn(ii) ffn(jj)];
                        end
                    end
                end
            end
            if fl>1
                TF = 1;
            elseif fl == 1; %check if the neighboring pair has a shared neighbor
                neighbors1=findFlippedNeighbors(p,idx(1),n);
                neighbors2=findFlippedNeighbors(p,idx(2),n);
                if ~isempty(intersect(neighbors1,neighbors2))
                    TF = 1;
                end
            end
        elseif Lffn == 4    % they will form a quadruplet if 1) all of them are connected, 2) three of them are connected, 3) if there are two pairs of connected ones either one of which share a neighbor
            fl = 0;
            idx=[];
            for ii = 1 : (Lffn-1)
                for jj = (ii+1): Lffn
                    if (sum((p.xy(ffn(ii),:) - p.xy(ffn(jj),:)).^2) < (1.01*a)^2)
                        fl = fl+1;
                        idx = [idx; ffn(ii) ffn(jj)];
                    end
                end
            end
            if fl == 3
                TF = 1;
            elseif fl == 2
                pp = idx(:);
                ol = 0;
                for ii = 1:(length(pp)-1)
                    for jj=(ii+1):length(pp)
                        if pp(ii)==pp(jj)
                            ol=1;
                            break
                        end
                    end
                end
                if ol==0
                    for ii=1:fl
                        neighbors1=findFlippedNeighbors(p,idx(ii,1),n);
                        neighbors2=findFlippedNeighbors(p,idx(ii,2),n);
                        if ~isempty(intersect(neighbors1,neighbors2))
                            ol=1;
                            break
                        end
                    end
                end                               
                if ol==1
                    TF = 1;
                end
            end
        end                                   
    case 5
        if Lffn>4
            TF = 1;
        elseif Lffn == 2 % the two neighbors have to have a shared neighbor and one of them at least has to have a shared neighbor with their shared neighbor that is not one of the orginal ones
            if (sum((p.xy(ffn(1),:) - p.xy(ffn(2),:)).^2) < (1.01*a)^2)
                neighbors1=findFlippedNeighbors(p,ffn(1),n);
                neighbors2=findFlippedNeighbors(p,ffn(2),n);
                sn=intersect(neighbors1,neighbors2);
                if ~isempty(sn)
                    neighbors3=findFlippedNeighbors(p,sn,n);
                    neighbors3new=[];
                    for ii=1:length(neighbors3)
                        if isempty(intersect(neighbors3(ii),ffn))
                            neighbors3new=[neighbors3new,neighbors3(ii)];
                        end
                    end
                    if ~isempty(intersect(neighbors1,neighbors3new)) || ~isempty(intersect(neighbors2,neighbors3new))
                        TF = 1;
                    end
                end
            end
        elseif Lffn == 3 % if it has 3 connected neighbors, either pair of connected neighbors has to have a shared neighbor; if it has 2 connected neighbors, it is like Lffn == 2 case above; if it has no connected neighbors keep going
            fl = 0;
            idx=[];
            for ii = 1 : (Lffn-1)
                for jj = (ii+1): Lffn
                    if (sum((p.xy(ffn(ii),:) - p.xy(ffn(jj),:)).^2) < (1.01*a)^2)
                        fl = fl+1;
                        idx = [idx; ffn(ii) ffn(jj)];
                    end
                end
            end
            if fl == 2
                for ii = 1:length(idx)
                    neighbors1=findFlippedNeighbors(p,idx(ii,1),n);
                    neighbors2=findFlippedNeighbors(p,idx(ii,2),n);
                    if ~isempty(intersect(neighbors1,neighbors2))
                        TF = 1;
                        break
                    end
                end
            elseif fl == 1
                neighbors1=findFlippedNeighbors(p,idx(1),n);
                neighbors2=findFlippedNeighbors(p,idx(2),n);
                sn=intersect(neighbors1,neighbors2);
                if ~isempty(sn)
                    neighbors3=findFlippedNeighbors(p,sn,n);
                    neighbors3new=[];
                    for ii=1:length(neighbors3)
                        if isempty(intersect(neighbors3(ii),idx(:)))
                            neighbors3new=[neighbors3new,neighbors3(ii)];
                        end
                    end
                    if ~isempty(intersect(neighbors1,neighbors3new)) || ~isempty(intersect(neighbors2,neighbors3new))
                    TF = 1;
                    end
                end
            end
        elseif Lffn == 4 % if all 4 connected, done; 3 connected plus 1 - does either of the connected pairs have a shared neighbor; 2 and 2 connected - run the Lffn == 2 case for each pair
            fl = 0;
            idx=[];
            for ii = 1 : (Lffn-1)
                for jj = (ii+1): Lffn
                    if (sum((p.xy(ffn(ii),:) - p.xy(ffn(jj),:)).^2) < (1.01*a)^2)
                        fl = fl+1;
                        idx = [idx; ffn(ii) ffn(jj)];
                    end
                end
            end   
            if fl == 3
                TF = 1;
            elseif fl == 2
                pp = idx(:);
                ol = 0;
                for ii = 1:(length(pp)-1)
                    for jj=(ii+1):length(pp)
                        if pp(ii)==pp(jj)
                            ol=1;
                            break
                        end
                    end
                end
                if ol==1
                    for ii=1:length(idx)
                        neighbors1=findFlippedNeighbors(p,idx(ii,1),n);
                        neighbors2=findFlippedNeighbors(p,idx(ii,2),n);
                        if ~isempty(intersect(neighbors1,neighbors2))
                            TF = 1;
                            break
                        end
                    end
                else
                    for ii=1:length(idx)
                        neighbors1=findFlippedNeighbors(p,idx(ii,1),n);
                        neighbors2=findFlippedNeighbors(p,idx(ii,2),n);
                        sn=intersect(neighbors1,neighbors2);
                        if ~isempty(sn)
                            neighbors3=findFlippedNeighbors(p,sn,n);
                            neighbors3new=[];
                            for jj=1:length(neighbors3)
                                if isempty(intersect(neighbors3(jj),idx(ii,:)))
                                    neighbors3new=[neighbors3new,neighbors3(jj)];
                                end
                            end
                    
                            if ~isempty(intersect(neighbors1,neighbors3new)) || ~isempty(intersect(neighbors2,neighbors3new))
                            TF = 1;
                            end
                        end
                    end
                end
            end
        end
    case 6 %case six-mer 1 to 3 (Methods Figure 1)
        if Lffn > 4
            TF = 1;
        elseif Lffn == 2 % either of the two adjacent neighbors had 4 neighbors before the last one flipped or the two adjacent neighbors have a shared neighbor and either both of them have an additional neighbor that's shared with their shared neighbor or only one of them has an additional neighbor that's shared with their shared neighbor but it and the shared neighbor share a neighbor
            if (sum((p.xy(ffn(1),:) - p.xy(ffn(2),:)).^2) < (1.01*a)^2)
                neighbors1=findFlippedNeighbors(p,ffn(1),n);
                neighbors2=findFlippedNeighbors(p,ffn(2),n);
                if length(neighbors1)>3 || length(neighbors2)>3
                    TF=1;
                elseif ~isempty(intersect(neighbors1,neighbors2))
                    sn=intersect(neighbors1,neighbors2);
                    neighbors3=findFlippedNeighbors(p,sn,n);
                    neighbors3new=[]; % will include neighbors of the shared neighbor in addition to the ffn
                    for ii=1:length(neighbors3)
                        if isempty(intersect(neighbors3(ii),ffn))
                            neighbors3new=[neighbors3new,neighbors3(ii)];
                        end
                    end
                    if ~isempty(intersect(neighbors1,neighbors3new)) && ~isempty(intersect(neighbors2,neighbors3new))
                        TF = 1;
                    elseif ~isempty(intersect(neighbors1,neighbors3new)) || ~isempty(intersect(neighbors2,neighbors3new))
                        if ~isempty(intersect(neighbors1,neighbors3new))
                            sn2=intersect(neighbors1,neighbors3new);
                        elseif ~isempty(intersect(neighbors2,neighbors3new))
                            sn2=intersect(neighbors2,neighbors3new);
                        end
                        neighbors4=findFlippedNeighbors(p,sn2,n);
                        neighbors4new=[];
                        for ii=1:length(neighbors4)
                            if isempty(intersect(neighbors4(ii),ffn)) && isempty(intersect(neighbors4(ii),sn))
                            neighbors4new=[neighbors4new,neighbors4(ii)];
                            end
                        end
                        if ~isempty(intersect(neighbors4new,neighbors3new))
                            TF=1;
                        end
                    end
                end
            end
        elseif Lffn == 3 % if none connected keep going, if 2 connected as in Lffn == 2, if 3 connected neighbors and each pair of connected neighbors shares a neighbor, TF=1; or if 3 connected neighbors and either pair of connected neighbors has a shared neighbor, which also shares a neighbor with either` ffn
            fl = 0;
            idx = [];
            for ii = 1 : (Lffn-1)
                for jj = (ii+1): Lffn
                    if (sum((p.xy(ffn(ii),:) - p.xy(ffn(jj),:)).^2) < (1.01*a)^2)
                        fl = fl+1;
                        idx=[idx, ffn(ii), ffn(jj)];
                    end
                end
            end
            if  fl == 1; %check if the neighboring pair has a shared neighbor
                neighbors1=findFlippedNeighbors(p,idx(1),n);
                neighbors2=findFlippedNeighbors(p,idx(2),n);
                if length(neighbors1)>3 || length(neighbors2)>3
                    TF=1;
                elseif ~isempty(intersect(neighbors1,neighbors2))
                    sn=intersect(neighbors1,neighbors2);
                    neighbors3=findFlippedNeighbors(p,sn,n);
                    neighbors3new=[];
                    for ii=1:length(neighbors3)
                        if isempty(intersect(neighbors3(ii),ffn))
                            neighbors3new=[neighbors3new,neighbors3(ii)];
                        end
                    end
                    if ~isempty(intersect(neighbors1,neighbors3new)) && ~isempty(intersect(neighbors2,neighbors3new))
                        TF = 1;
                    elseif ~isempty(intersect(neighbors1,neighbors3new)) || ~isempty(intersect(neighbors2,neighbors3new))
                        if ~isempty(intersect(neighbors1,neighbors3new))
                            sn2=intersect(neighbors1,neighbors3new);
                        elseif ~isempty(intersect(neighbors2,neighbors3new))
                            sn2=intersect(neighbors2,neighbors3new);
                        end
                        neighbors4=findFlippedNeighbors(p,sn2,n);
                        neighbors4new=[];
                        for ii=1:length(neighbors4)
                            if isempty(intersect(neighbors4(ii),ffn)) && isempty(intersect(neighbors4(ii),sn))
                            neighbors4new=[neighbors4new,neighbors4(ii)];
                            end
                        end
                        if ~isempty(intersect(neighbors4new,neighbors3new))
                            TF=1;
                        end
                    end
                end
            elseif fl == 2;
                Lidx = length(idx);
                for ii = 1 : (Lidx-1)
                    for jj = (ii+1) : Lidx
                        if idx(ii) == idx(jj)
                            middle = idx(ii);
                        end
                    end
                end
                side = [];
                for ii = 1 : Lidx
                    if ~(idx(ii) == middle)
                        side = [side, idx(ii)];
                    end
                end
                neighbors1=findFlippedNeighbors(p,middle,n);
                neighbors2=findFlippedNeighbors(p,side(1),n);
                neighbors3=findFlippedNeighbors(p,side(2),n);
                if length(neighbors1)>3 || length(neighbors2)>3 || length(neighbors3)>3
                    TF=1;
                elseif length(neighbors2) == 3 && length(neighbors3) == 3
                    fl = 0;
                    for ii = 1 : 2
                        for jj = (ii+1) : 3
                            if (sum((p.xy(neighbors2(ii),:) - p.xy(neighbors2(jj),:)).^2) < (1.01*a)^2)
                                fl = fl+1;
                            end
                        end
                    end
                        if fl > 1
                            TF = 1;
                        else fl = 0;
                            for ii = 1 : 2
                                for jj = (ii+1) : 3
                                    if (sum((p.xy(neighbors3(ii),:) - p.xy(neighbors3(jj),:)).^2) < (1.01*a)^2)
                                    fl = fl+1;
                                    end
                                end
                            end
                        end
                        if fl > 1
                            TF = 1;
                        end
                elseif length(neighbors2) == 3
                    fl = 0;
                    for ii = 1 : 2
                        for jj = (ii+1) : 3
                            if (sum((p.xy(neighbors2(ii),:) - p.xy(neighbors2(jj),:)).^2) < (1.01*a)^2)
                                fl = fl+1;
                            end
                        end
                    end
                    if fl > 1
                        TF = 1;
                    end
                elseif length(neighbors3) == 3
                    fl = 0;
                    for ii = 1 : 2
                        for jj = (ii+1) : 3
                            if (sum((p.xy(neighbors3(ii),:) - p.xy(neighbors3(jj),:)).^2) < (1.01*a)^2)
                                fl = fl+1;
                            end
                        end
                    end
                    if fl > 1
                        TF = 1;
                    end
                end
            end
        elseif Lffn == 4 % if 4 connected, than any pair with shared neighbor will yield TF=1, if 3+1 then the case of 3 connected, if 2+2 then two cases of 2 either can yield TF=1
            fl = 0;
            idx = [];
            for ii = 1 : (Lffn-1)
                for jj = (ii+1): Lffn
                    if (sum((p.xy(ffn(ii),:) - p.xy(ffn(jj),:)).^2) < (1.01*a)^2)
                        fl = fl+1;
                        idx=[idx, ffn(ii), ffn(jj)];
                    end
                end
            end
            Lidx = length(idx);
            if fl>2
                for ii = 1 : 2 : (Lidx-1) 
                    jj = ii+1;
                    neighbors1=findFlippedNeighbors(p,idx(ii),n);
                    neighbors2=findFlippedNeighbors(p,idx(jj),n);
                    if ~isempty(intersect(neighbors1,neighbors2))
                        TF=1;
                    end
                end
            elseif fl == 2
                middle = [];
                for ii = 1 : (Lidx-1)
                    for jj = (ii+1) : Lidx
                        if idx(ii) == idx(jj)
                            middle = idx(ii);
                        end
                    end
                end
                if ~isempty(middle)
                    side = [];
                    for ii = 1 : Lidx
                        if ~(idx(ii) == middle)
                            side = [side, idx(ii)];
                        end
                    end
                    neighbors1=findFlippedNeighbors(p,middle,n);
                    neighbors2=findFlippedNeighbors(p,side(1),n);
                    neighbors3=findFlippedNeighbors(p,side(2),n);
                    if length(neighbors1)>3 || length(neighbors2)>3 || length(neighbors3)>3
                        TF=1;
                    elseif length(neighbors2) == 3 && length(neighbors3) == 3
                        fl = 0;
                        for ii = 1 : 2
                            for jj = (ii+1) : 3
                                if (sum((p.xy(neighbors2(ii),:) - p.xy(neighbors2(jj),:)).^2) < (1.01*a)^2)
                                    fl = fl+1;
                                end
                            end
                        end
                            if fl > 1
                                TF = 1;
                            else fl = 0;
                                for ii = 1 : 2
                                    for jj = (ii+1) : 3
                                        if (sum((p.xy(neighbors3(ii),:) - p.xy(neighbors3(jj),:)).^2) < (1.01*a)^2)
                                            fl = fl+1;
                                        end
                                    end
                                end
                            end
                            if fl > 1
                                TF = 1;
                            end
                    elseif length(neighbors2) == 3
                        fl = 0;
                        for ii = 1 : 2
                            for jj = (ii+1) : 3
                                if (sum((p.xy(neighbors2(ii),:) - p.xy(neighbors2(jj),:)).^2) < (1.01*a)^2)
                                    fl = fl+1;
                                end
                            end
                        end
                        if fl > 1
                            TF = 1;
                        end
                    elseif length(neighbors3) == 3
                        fl = 0;
                        for ii = 1 : 2
                            for jj = (ii+1) : 3
                                if (sum((p.xy(neighbors3(ii),:) - p.xy(neighbors3(jj),:)).^2) < (1.01*a)^2)
                                    fl = fl+1;
                                end
                            end
                        end
                        if fl > 1
                            TF = 1;
                        end
                    end
                else for ii = 1 : 2 : Lidx-1
                        for jj = (ii+1) : 2 : Lidx
                            neighbors1=findFlippedNeighbors(p,idx(ii),n);
                            neighbors2=findFlippedNeighbors(p,idx(jj),n);
                            if length(neighbors1)>3 || length(neighbors2)>3
                                TF=1;
                            elseif ~isempty(intersect(neighbors1,neighbors2))
                                sn=intersect(neighbors1,neighbors2);
                                neighbors3=findFlippedNeighbors(p,sn,n);
                                neighbors3new=[]; % will include neighbors of the shared neighbor in addition to the ffn
                                for ii=1:length(neighbors3)
                                    if isempty(intersect(neighbors3(ii),idx))
                                    neighbors3new=[neighbors3new,neighbors3(ii)];
                                    end
                                end
                                if ~isempty(intersect(neighbors1,neighbors3new)) && ~isempty(intersect(neighbors2,neighbors3new))
                                    TF = 1;
                                elseif ~isempty(intersect(neighbors1,neighbors3new)) || ~isempty(intersect(neighbors2,neighbors3new))
                                    if ~isempty(intersect(neighbors1,neighbors3new))
                                        sn2=intersect(neighbors1,neighbors3new);
                                    elseif ~isempty(intersect(neighbors2,neighbors3new))
                                        sn2=intersect(neighbors2,neighbors3new);
                                    end
                                    neighbors4=findFlippedNeighbors(p,sn2,n);
                                    neighbors4new=[];
                                    for kk=1:length(neighbors4)
                                        if isempty(intersect(neighbors4(kk),idx)) && isempty(intersect(neighbors4(kk),sn))
                                        neighbors4new=[neighbors4new,neighbors4(kk)];
                                        end
                                    end
                                    if ~isempty(intersect(neighbors4new,neighbors3new))
                                        TF=1;
                                    end
                                end
                            end
                        end
                    end
                    
                            
                end
                    
                
                     
        
        
        
        
        end
        
                
                        
                        
        
end
end            
                  
        
            
                
        
    
        
        
       



 