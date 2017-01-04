function [ dxArray, dxiArray, dxArray5, dxiArray5, tcellcycle, Delta_xi ]...
    = histStepCell( dirname )
% get steps size from the tracks

% define the number of bins in x and time
xbin = 9;
tbin = 3;

% large step size
dt5 = 10;

% Delta xi is the starting point and ending point for a track
Delta_xi = [];

% cell array of step size
dxArray  = cell(xbin,tbin);
dxiArray = cell(xbin,tbin);

% cell array of step size for large step
dxArray5  = cell(xbin,tbin);
dxiArray5 = cell(xbin,tbin);


% list of the cell cycle length
tcellcycle = [];


dirname_cell = [dirname,'/xy1/cell/'];
CONST = load( [dirname,'/CONST.mat'] );

% get the cell files
contents = dir( [dirname_cell,'Cell*.mat'] );

numCell = numel( contents );

for ii = 1:numCell
    
    disp('This many left');
    disp(ii/numCell);
    
    celld = load(  [dirname_cell,contents(ii).name] );
    
    %try
        data = getLocusTracksDev( celld, true, CONST);
        drawnow;
        
        ss = size(data.track1.xtrack)
        
        % Orient cells old pole versus new pole
        data.track1.xtrack = celld.CellA{1}.pole.op_ori * data.track1.xtrack;
        
        % calculate relative position in the cell
        xi = data.track1.xtrack./(data.track1.lx*ones(1,ss(2)))+.5;
   
        % calculate the bin in which this location belongs
        nn   = ceil(xi(1:end-1,:)*xbin);
        nn(nn>xbin) = xbin;
        nn(nn<1)    = 1;
        
        % calculate the step size
        % relative step size
        dxi = xi(2:end,:)-xi(1:end-1,:);
        % absolute stepsize
        dx = data.track1.xtrack(2:end,:)-data.track1.xtrack(1:end-1,:);
   
        % calculate the bin in which this location belongs
        nn5   = ceil(xi(1:(end-dt5),:)*xbin);
        nn5(nn5>xbin) = xbin;
        nn5(nn5<1)    = 1;
        
        % calculate the step size
        % relative step size
        dxi5 = xi((1+dt5):end,:)-xi(1:(end-dt5),:);
        % absolute stepsize
        dx5 = data.track1.xtrack((1+dt5):end,:)-data.track1.xtrack(1:(end-dt5),:);
        
        
        
        % time bin size
        DT = (ss(1)-1)/tbin;
        
        % cut off track score is set here to eliminate short tracks
        sumS_cut = 100;
        
        
        % ss(2) is the number of different tracks
        % loop through the tracks (jj)
        for jj = 1:ss(2);
            
            % check to make sure that the score for this track is greater
            % than the cut.
            if data.track1.sumS(jj) > sumS_cut
                
                % add the current cell cycle list to the vector of cell
                % cycle lengths
                tcellcycle = [tcellcycle, numel(celld.CellA)];
                
                % loop through the time bins
                for kk = 1:tbin
                    
                    % this make a vector of of the time steps (tt)
                    % corresponding to this bin
                    tt = floor((kk-1)*DT+1):floor((kk)*DT);
                    
                    % this make a vector of of the time steps (tt)
                    % corresponding to this bin
                    ss5 = size( dx5 );
                    tt5 = (round((kk-1)/tbin*ss5(1))+1):round((kk)/tbin*ss5(1));
                    
                    % loop through xbins
                    for ll = 1:xbin
                        
                        %% Do the short steps here
                        % get the dx steps and put them in the dxArray
                        tmp = dx(nn(tt,jj)==ll,jj);
                        tmp = tmp(~isnan(tmp));
                        dxArray{ll,kk}  = [dxArray{ll,kk}; tmp];
                        
                        % get the dxi (relative) steps and put them in 
                        % the dxArray
                        tmp = dxi(nn(tt,jj)==ll,jj);
                        tmp = tmp(~isnan(tmp));
                        dxiArray{ll,kk}  = [dxiArray{ll,kk}; tmp];
                        
                        %% do the long steps here
                        % get the dx steps and put them in the dxArray
                        tmp = dx5(nn5(tt5,jj)==ll,jj);
                        tmp = tmp(~isnan(tmp));
                        dxArray5{ll,kk}  = [dxArray5{ll,kk}; tmp];
                        
                        % get the dxi (relative) steps and put them in 
                        % the dxArray
                        tmp = dxi5(nn5(tt5,jj)==ll,jj);
                        tmp = tmp(~isnan(tmp));
                        dxiArray5{ll,kk}  = [dxiArray5{ll,kk}; tmp];                        
                        
                        
                   end
                end
            end
            
            % record starting and ending of the track
            tmp = xi(:,jj);
            tmp = tmp(~isnan(tmp));
            Delta_xi = [Delta_xi; tmp(1), tmp(end)];
        end
    %catch
    %    disp( 'Error' );
    %end
    
    
end

end

