%%% Adapted from Ivanovic et al. "Influenza-virus membrane fusion: cooperative fold-back of stochastically induced hemaggutinin intermediates"
%%% Revised for "Distinct functional determinants of influenza hemagglutinin-mediated membrane fusion" to include Nh=6 and
%%% the unproductive-HA population. Also, here, hemifusion time is measured from
%%% the start instead of from the arrest intermediate.
%%%
%%% Main Code Text for MATLAB version R2012a
%%%
%%% 1)  first get lag times for all trimers in a pre-defined 'patch' of size N sampled randomly from an exponentially decaying function
%%%
%%%     HAprefusion -> HAextended (lag times for this process will be described by an exponentially decaying function; 
%%%                                pre-define rate, k (k_sim in the paper), and initially get lag times for ALL trimers in a patch 
%%%                                (as if this was an entire reaction (all HAs are allowed to extend and hemifusion does not end the process)
%%% 2)  sort particles by their lag times
%%% 3)  find arrest and hemifusion times
%%%         arrest time: the time at which the N1th active HA flipped
%%%         (arrest was not studied in the current paper except to generate plots and results shown in Figure 2B and the model figures (first event delay: arrest time when N1=1 and inactive=unproductive=0);
%%%         hemifusion time: the time at which that HA flipped, which completed the first set of N2 active 'neighbors'                   
%%%         note that in Ivanovic et al., 2013, hemifusion time:(hemifusion time-arrest time)
%%% 4)  repeat 1 to 3 for Nvirions (number of virions analyzed in a single experiment) 
%%%     and collect arrest and hemifusion times in 2 vectors: arrest and hemi

Nvirions = input('How many virions would you like to analyze in this experiment : ');
N = input('Approximately how many HA trimers are in the contact patch (enter 50,100,250 for n-actual = 55,121,295, respectively) : '); %note actual number (Nactual) is bigger than yout input estimate; 
N1 = input('How many extended HA trimers lead to arrest : '); % N1 represents the number of HA trimers that need to extend and insert for arrest to take place;
N2 = input('How many extended HA neighbors lead to hemifusion (enter number between 2 and 6) : '); % 2<=N2<=6; N2 represents the number of neighboring HA trimers that need to extend and insert for hemifusion to take place
sf = input('What is the fraction of inactive HAs (enter number between 0 and 1): ');
np = input('What is the fraction unproductive HAs (enter number between 0 and 1): ');

Vplot = input('Would you like the code to plot final states for each analyzed virion [y/N] :', 's');
if(isempty(Vplot) || upper(Vplot)=='N')
    YesPlot=0;
else
    YesPlot=1;
end

k = 0.015;
%k = 0.02;
%k = 0.035;
%k = 0.0025;
Ndeadvirions=0;

tic
[p,Nactual,a] = generate_patch(N);   % circular patch of N elements arranged in a hexagonal lattice

for ii = 1:Nvirions
    NotDeadFlag=0;
    fl = 0; nins = 0;
    if YesPlot
        figure; plot(p.xy (:,1),p.xy(:,2),'mx', 'MarkerSize',20); axis image; set(gcf,'Color',[1 1 1]); hold on;
    end
    p.flipped = s_randomdist(k, Nactual);  % samples randomly from an exponentially decaying funcion with rate k; gets lag times for particles in order from 1st to Nactual-th
    
    [y,idx] = sortrows([p.id p.xy p.flipped],4); % sorts data according to lag times 
    p.id = y(:,1); p.xy = y(:,2:3); p.flipped = y(:,4);
    p.neighbors = p.neighbors(idx);   % Re-sort cell array.
    for kk = 1:length(p.neighbors)
        for jj = 1:length(p.neighbors{kk})
            newidx = find(p.neighbors{kk}(jj) == idx);
            p.neighbors{kk}(jj) = newidx;  
        end
    end
    
    p.inactive = zeros(Nactual,1);
    p.unproductive = zeros(Nactual,1);
        
    if sf~=0 && np==0
        r=randperm(Nactual);
        for jj = 1:ceil((sf * (Nactual)))
            p.inactive(r(jj))=1; 
        end
    elseif sf==0 && np~=0
        r=randperm(Nactual);
        for jj = 1:ceil((np * (Nactual)))
            p.unproductive(r(jj))=1; 
        end
    elseif sf~=0 && np~=0
        r=randperm(Nactual);
        ia = ceil(sf * (Nactual));
        rm = ceil(np * ((Nactual)-ia));
        for jj = 1:ia
            p.inactive(r(jj))=1; 
        end
        for jj=(ia+1):(ia+rm)
            p.unproductive(r(jj))=1;
        end
    end
                
    for jj = 1: Nactual
         if p.inactive(jj) == 1
            if YesPlot
                plot(p.xy(jj,1),p.xy(jj,2),'kx', 'MarkerSize',20); hold on;
                plot(p.xy(jj,1),p.xy(jj,2),'ks', 'MarkerFaceColor', 'k', 'MarkerSize',7); hold on;
            end
         end
         if p.unproductive(jj) == 1
            if YesPlot
                plot(p.xy(jj,1),p.xy(jj,2),'ks', 'MarkerFaceColor', 'c', 'MarkerSize',7); hold on;
            end
         end
    end
    
    for jj = 1: Nactual
        if p.inactive(jj) == 0 && p.unproductive(jj) == 0
            fl=fl+1;
        end
          
        if fl==N1;
            tN1 = p.flipped(jj);
            break
        end
    end       
     
    for n = 1:Nactual   % go through all HAs and mark active and unproductive events until hemifusion; starting with N2th active HA, check for neighborhoods
        if p.inactive(n)==0 && p.unproductive(n)==0
            nins=nins+1;
            if YesPlot
                plot(p.xy(n,1), p.xy(n,2), 'mo', 'MarkerSize',20); hold on;
                if nins<N1
                    plot(p.xy(n,1), p.xy(n,2), 'mo', 'MarkerSize',10); hold on;
                elseif nins==N1
                    plot(p.xy(n,1), p.xy(n,2), 'mo', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm','MarkerSize',10); hold on;
                end
            end
            
            if nins>=N2 
                if isaN2tuplet6AllGeos(p, n, N2, a)
                    NotDeadFlag=1;
                    tN2=p.flipped(n);
                    if YesPlot
                        plot(p.xy(n,1), p.xy(n,2),'mo', 'MarkerFaceColor','m', 'MarkerSize',20);
                    end
                    break;
                end
            end
        elseif YesPlot && p.unproductive(n) == 1
            plot(p.xy(n,1), p.xy(n,2), 'k*', 'MarkerSize',20); hold on;
            plot(p.xy(n,1), p.xy(n,2),'ks', 'MarkerFaceColor', 'c', 'MarkerSize',7); hold on;
        end   
    end
    
    allarrest(ii) = tN1;
    if NotDeadFlag==0
        Ndeadvirions=Ndeadvirions+1;
        allhemi(ii) = NaN;
    else
%         dt = tN2-tN1; % in Ivanovic et al., 2013 study, hemifusion time was defined as the lag time between arrest and hemifusion
 %       allhemi(ii) = dt;
        allhemi(ii) = tN2;
    end    
end
toc

arrest1=allarrest(allhemi>allarrest); %analyze those virions that stop moving before hemifusing
arrest2=allarrest(isnan(allhemi)); %analyze those virions that arrest but that do not hemifuse
arrest=[arrest1,arrest2];
hemi=allhemi(allhemi>0); %include only hemifusion competent virions
% hemi=allhemi(allhemi>allarrest); %includes hemifusion competent virions that arrest before hemifusing

   
