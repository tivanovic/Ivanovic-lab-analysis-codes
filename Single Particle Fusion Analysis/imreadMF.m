function A = imreadMF(movielengths,filename,filetype,nframe,varargin)
% Run imread on multi-TIF movie file.

cumulMovieFrameCount = [0 cumsum(movielengths)];
nmovie = find(nframe <= cumulMovieFrameCount(2:end)); nmovie = nmovie(1);
nframelocal = nframe - cumulMovieFrameCount(nmovie);
A = imread(filename{nmovie},filetype,nframelocal,varargin{:});
