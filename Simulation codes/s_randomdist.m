function [tvec] = s_randomdist(k,N)
% Samples randomly from an exponentially decaying function

tvec = random('exp',1/k,N,1);

end
