function [auto] = getVelVelCorr( dirname )

%% THIS FUNTION IS MEANT TO MATCH THE VELOCITY-VELOCITY CALCUALTION MADE IN
% STEPH WEBER'S PRL
%
% NJK. MAY 2011.

%% INITIALIZE CONSTANTS

global CONST

if isempty( CONST )
    if exist('loadConstantsMine','file');
        loadConstantsMine
    else
        loadConstantsMine
    end
end

%% GET TRAJECTORIES FOR EACH CELL IN DIRECTORY dirname

[ tr,track_mean, variance ] = getTrackMean( dirname );
 
%% CALCUALTE VELOCITY FUNTION FOR EACH TRAJECTORY
%
% NOTE: THIS IS ONLY INTERESTED IN THE ORIGIN TRACKS: ri, rh, rl

ri = tr.ri' ;
rh = tr.rh' ;
rl = tr.rl' ;

ssrh = size(rh)


for jj = 1:10;
    
    tau = jj;
    
    ri_aut(jj) = getVelVel( ri, tau ) ;
    rh_aut(jj) = getVelVel( rh(1:20,:), tau ) ;
    rl_aut(jj) = getVelVel( rl(1:20,:), tau ) ;
    
end

auto.ri = ri_aut./ri_aut(1) ;
auto.rh = rh_aut./rh_aut(1)  ;
auto.rl = rl_aut./rl_aut(1) ;


end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function velvel = getVelVel( mat, tau )

% NOTE: THIS IS FOR DELTA = 1 MIN

ss = size(mat) ; 

for ii = 1:ss(2)
   
    v0 = mat(2:ss(1)) - mat(1:(ss(1)-1));
    
    ss0 = length(v0) - tau +1 ;
    
    VV = v0(tau:end).*v0(1:ss0);

    VV_mean(ii) = mean(VV(~isnan(VV)));
    
end

    velvel = mean(VV_mean(~isnan(VV_mean)));

end