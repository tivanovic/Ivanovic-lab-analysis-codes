function [ffn] = findFlippedNeighbors(p,n,k)
%goes through the neighbors of the nth element (the one that just flipped)
%and returns the ones that have flipped before it into ffn (find flipped neighbors) 

ffn=[];
idx=p.neighbors{n};

for ii = 1:length(idx)
    if (p.flipped(idx(ii))<p.flipped(k)) && (p.inactive(idx(ii))==0) && (p.unproductive(idx(ii))==0)
       ffn=[ffn;idx(ii)];
    end
end
end

