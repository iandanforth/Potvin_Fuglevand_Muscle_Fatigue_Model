% Motor Unit Based Muscle Fatigue Model by Jim Potvin & Andrew Fuglevand
% front end (rested size-principle) based on Fuglevand, Winter & Patla (1993)
% last updated 2017-05-28 by Jim Potvin
% 
% Partial vectorization by Ian Danforth

% fthtime = 10
%   - original matlab code:     Run time: 1372.633170 seconds
%   - vectorized matlab code:   Run time: 19.214544 seconds ~ ~70x
%   - python code:              Run time: 3.785 seconds - ~360x


%clear all
clc
start_time = cputime;

%% Model input parameters
nu = 120;           % number of neurons (ie. motor units) in the modeled pool ("n") 

samprate = 10;      % sample rate (10 Hz is suggested)
res = 100;          % resolution of activations (set = 10 for 0.1 activation resolution, 100 for 0.01)
hop = 20;           % allows for hopping through activations to more rapidly find that which meets the threshold (10 means every 1/10th of maxact)
r = 50;             % range of recruitment thresholds (30 or 50)

fat = 180;          % range of fatigue rates across the motor units (300 best) (fat)
FatFac = 0.0225;    % fatigue factor (FF/S) percent of peak force of MU per second (FatFac)

tau = 22;           % 22 based on Revill & Fuglevand (2011)
adaptSF = 0.67;     % 0.67 from Revill & Fuglevand (2011)
ctSF = 0.379;       % 0.379 based on Shields et al (1997)
 
mthr = 1;           % minimum recruitment threshold
a = 1;              % recruitment gain paramter (1)
minfr = 8;          % minimum firing rate (8)
pfr1 = 35;          % peak firing rate of first recruited motor unit (35)
pfrL = 25;          % peak firing rate of last recruited motor unit (25)
mir = 1;            % slope of the firing rate increase vs excitation (1)
rp = 100;           % range of twitch tensions (100)
rt = 3;             % range of contraction times (3 fold)
tL = 90;            % longest contraction time (90)

%% Various methods to create, or read in, force (%MVC)time-histories 

%% Create isotonic data -----------------------------------


% Excitation Level
fthscale = 0.5;             % sets %MVC level for the trial duration (100% MVC is 1.00)
con = '0.50';               % for output file names
fthtime = 10;              % duration to run trial (seconds)

fthsamp = fthtime * samprate; 
excitations = zeros(1, fthsamp);
for z = 1:fthsamp
    excitations(z) = fthscale;
end

%% Calculations from the Fuglevand, Winter & Patla (1993) Model 
ns = 1:1:fthsamp;   % array of samples for fth
excitations = excitations(ns);

% Recruitment Threshold Excitation (thr)
thr = zeros(1, nu);
n = 1:1:nu;
b = log(r + (1-mthr)) / (nu-1);             % this was modified from Fuglevand et al (1993) RTE(i) equation (1)
for i = 1:nu                                % as that did not create the exact range of RTEs (ie. 'r') entered
    thr(i) = a * exp((i-1) * b) - (1 - mthr);
end

% Peak Firing Rate (frp)
% modified from Fuglevand et al (1993) PFR equation (5) to remove thr(1) before ratio
frdiff = pfr1 - pfrL;
frp = pfr1 - frdiff * ((thr(n) - thr(1))/(r - thr(1)));
maxex = thr(nu) + (pfrL - minfr)/mir;       % maximum excitation
maxact = round(maxex * res);                % max excitation x resolution
ulr = 100 * thr(nu)/maxex;                  % recruitment threshold (%) of last motor unit


% Calculation of the rested motor unit twitch properties (these will change with fatigue)

% PRE-CALCULATE FIRING RATES FOR EACH MOTOR UNIT BY EXCITATION LEVEL

% Firing Rates for each MU with increased excitation (act)
mufr = zeros(nu, maxact);
for act = 1:maxact
    acti = act / res; % Each step is actually a step of <resolution> size.
    above_thresh_indices = acti >= thr;
    new_vals = mir * (acti - thr) + minfr;
    above_peak_indices = new_vals > frp;
    new_vals(above_peak_indices) = frp(above_peak_indices);
    mufr(above_thresh_indices, act) = new_vals(above_thresh_indices);
end

% UNVECTORIZED CODE
% for mu = 1:nu
%     for act = 1:maxact
%         acti = act / res; % Each step is actually a step of <resolution> size.
%         if acti >= thr(mu) % If the activation level is equal to or above threshold
%             mufr(mu, act) = mir * (acti - thr(mu)) + minfr;  
%             if mufr(mu, act) > frp(mu)
%                 mufr(mu, act) = frp(mu);
%             end
%         elseif acti < thr(mu)
%             mufr(mu, act) = 0;
%         end
%     end
% end


k = 1:1:maxact;  %range of excitation levels

% Twitch peak force (P)
b = log(rp)/(nu-1);                 % this was modified from Fuglevand et al (1993) P(i) equation (13)
P = exp(b * (n-1));                 % as that didn't create the exact range of twitch tensions (ie. 'rp') entered

% Twitch contraction time (ct) 
c = log(rp)/log(rt);                % scale factor
ct = zeros(1,nu);        
for mu = 1:nu                       % assigns contraction times to each motor unit (moved into loop)
   ct(mu) = tL * (1/P(mu))^(1/c);
end

% Normalized motor unit firing rates (nmufr) with increased excitation (act)
nmufr = zeros(nu, maxact);
for act = 1:maxact
    nmufr(:, act) = ct' .* (mufr(:, act) / 1000);  % same as CT / ISI  
end

% UNVECTORIZED CODE
% for mu = 1:nu                    
%     for act = 1:maxact
%         nmufr(mu, act) = ct(mu) * (mufr(mu, act) / 1000);  % same as CT / ISI  
%     end
% end

% Normalized Force
% Motor unit force, relative to full fusion (Pr) with increasing excitation
% based on Figure 2 of Fuglevand et al (1993)
sPr = 1 - exp(-2 * (0.4^3)); % NOTE - THIS IS JUST SIMPLIFIED TO 0.3 IN THE PAPER
Pr = zeros(nu, maxact);
for act = 1:maxact
    below_thresh_indices = nmufr(:, act) <= 0.4;
    Pr(below_thresh_indices, act) = nmufr(below_thresh_indices,act)/0.4 * sPr;
    above_thresh_indices = nmufr(:, act) > 0.4;
    Pr(above_thresh_indices,act) = 1 - exp(-2 * (nmufr(above_thresh_indices,act).^3));
end

% UNVECTORIZED CODE
% for mu = 1:nu                       
%     for act = 1:maxact
%         if nmufr(mu,act) <=0.4                      %linear portion of curve
%             Pr(mu,act) = nmufr(mu,act)/0.4 * sPr;   %Pr = MU force relative to rest 100% max excitation of 67
%         end
%         if nmufr(mu,act) > 0.4
%             Pr(mu,act) = 1 - exp(-2 * (nmufr(mu,act)^3));
%         end
%     end
% end

% Motor unit force (muP) with increased excitation
muP = zeros(nu, maxact);
for act = 1:maxact 
   muP(:,act) = Pr(:,act) .* P';
end
% UNVECTORIZED CODE
% for mu = 1:nu
%     for act = 1:maxact 
%        muP(mu,act) = Pr(mu,act) * P(mu);
%     end
% end

totalP = sum(muP,1);                                % sum of forces across MUs for each excitation (dim 1)
maxP = totalP(maxact); % CHECKED


% Total Force across all motor units when rested
Pnow = zeros(nu, fthsamp);
Pnow(:,1) = P(:);


%% Calculation of Fatigue Parameters (recovery currently set to zero in this version)
    
% note, if rp = 100 & fatigue_range = 180, there will be a 100 x 180 = 1800-fold difference in
% the absolute fatigue of the highest threshold vs the lowest threshold.
% The highest threshold MU will only achieve ~57% of its maximum (at 25 Hz), so the actual range of fatigue
% rates is 1800 x 0.57 = 1026 

% fatigue rate for each motor unit  (note: "log" means "ln" in Matlab)
b2 = log(fat)/(nu-1);       
mufatrate = exp(b2 * (n-1)); % Curve from 1 to 180

% NOMINAL FATIGUE RATES
fatigue = zeros(1,nu); % VERIFIED
for mu = 1:nu
     fatigue(mu) = mufatrate(mu) * (FatFac / fat) * P(mu);    
        % the full fatigue rate is mufatrate(mu) * [max_fatigue_rate / fatigue_range] * Pr(mu,act) * P(mu)
        % the only variable is the relative force: Pr(mu,act), so this part is calculated once here
end

% Establishing the rested excitation required for each target load level (if 0.1% resolution, then 0.1% to 100%)    
startact = zeros(1, 100); % VERIFIED
percent_total_force = totalP / maxP;
for force = 1:100
    startact(force) = 0;    % excitation will never be lower than that needed at rest for a given force
    indices = find((percent_total_force * 100) >= force);
    first_act_level = indices(1);
    startact(force) = first_act_level - 2;
    % UNVECTORIZED CODE
    % for act = 1:maxact      % so it speeds the search up by starting at this excitation
    %     if (totalP(act)/maxP * 100) < force
    %         startact(force) = act - 1;
    %     end        
    % end
end

% Pchangecurves = zeros(nu, maxact);
% for act = 1:maxact
%     for mu = 1:nu
%         Pchangecurves(mu,act) = fatigue(mu) * Pr(mu, act) * P(mu);  % just used for graphical display
%     end
% end

mes = 'start of fatigue analysis'


%% Moving through force time-history and determing the excitation required to meet the target force at each time

TmuPinstant = zeros(nu, fthsamp); 
m = zeros(1, fthsamp);
mufrFAT = zeros(nu, fthsamp);
ctFAT = zeros(nu, fthsamp);
ctREL = zeros(nu, fthsamp);
nmufrFAT = zeros(nu, fthsamp);
PrFAT = zeros(nu, fthsamp);
muPt = zeros(nu, fthsamp);
TPt = zeros(nu, fthsamp);
Ptarget = zeros(1, fthsamp);
Tact = zeros(1, fthsamp);
Pchange = zeros(nu, fthsamp); 
TPtMAX = zeros(1, fthsamp);
muPtMAX = zeros(nu, fthsamp); 
muON = zeros(1, nu);
adaptFR = zeros(nu,fthsamp);
Rdur = zeros(1,nu);
% acttemp = zeros(fthsamp, maxact);
muPna = zeros(nu, fthsamp);
muForceCapacityRel = zeros(nu,fthsamp);
PrMAX = zeros(1, nu);
timer = 0;

fthsamp
   
for i = 1:fthsamp    
    if i == (timer + 1) * samprate * 60         % shows a timer value every 15 seconds
        timer = timer + 1;
        current = i / samprate
    end
    
    force = round(excitations(i) * 100) + 1;      % used to start at the minimum possible excitation (lowest it can be is 1)
    if force > 100                                % so start with excitation needed for fth(i) when rested (won't be lower than this)
        force = 100;
    end
    s = startact(force) - (5 * res);            % starts a little below the minimum (0.95)
    if s < 1
        s = 1;
    end

    acthop = round(maxact / hop);               % resets 'acthop' to larger value for new sample
    act = s;                                    % start at lowest value then start jumping by 'acthop'
    while true                                  % this starts at the mimimum (s) then searches for excitation required to meet the target
        % For this timestep, activation level, store our starting activation (in the first loop)
        % acttemp(i, a) = act; % never used even in original code

        % Update our active durations for each motor unit
        above = muON > 0;
        Rdur(above) = (i - muON(above) + 1) / samprate;
        % Zero out durations if we go below zero for some reason
        Rdur(Rdur < 0) = 0;

        % Firing rate adaptation curve
        % q(i)(Eq. 13)
        ratio = (thr - 1) / (thr(nu) - 1);
        q = adaptSF * (mufr(:, act) - minfr + 2)' .* ratio;

        % (Eq. 12)
        intrinsic = 1 - exp(-1 * Rdur / tau);
        adaptFR(:,i) = q .* intrinsic;
        adaptFR(adaptFR(:,i) < 0,i) = 0;

        % Motor unit firing rates with fatigue
        mufrFAT(:, i) = mufr(:, act) - adaptFR(:, i); % adapted motor unit firing rate: based on time since recruitment
        mufrMAX = mufr(:, maxact) - adaptFR(:, i); % adapted FR at max excitation

        % Update contraction times
        ctFAT(:, i) = ct .* (1 + ctSF * (1 - (Pnow(:, i) ./ P')))';  % corrected contraction time: based on MU fatigue
        % ctREL(:, i) = ctFAT(:, i) ./ ct'; % unused

        % Calculate normalized firing rate
        nmufrFAT(:, i) = ctFAT(:, i) .* (mufrFAT(:, i) / 1000);  % adapted normalized Stimulus Rate (CT * FR)
        nmufrMAX = ctFAT(:, i) .* (mufrMAX / 1000);  % normalized FR at max excitation

        % Calculate normalized forces
        below = nmufrFAT(:, i) <= 0.4;
        PrFAT(below, i) = nmufrFAT(below, i) / 0.4 * sPr;  % fusion level at adapted firing rate (linear portion)
        above = nmufrFAT(:, i) > 0.4;
        PrFAT(above, i) = 1 - exp(-2 * (nmufrFAT(above, i).^3));  % Non linear portion of curve

        % Calculate current force
        muPt(:, i) = PrFAT(:, i) .* Pnow(:, i); % MU force at the current time (muPt): based on adapted postion on fusion curve

        % fusion force at 100% maximum excitation
        below = nmufrMAX <= 0.4;
        PrMAX(below) = nmufrMAX(below) / 0.4 * sPr;
        above = nmufrMAX > 0.4;
        PrMAX(above) = 1 - exp(-2 * (nmufrMAX(above).^3));

        muPtMAX(:, i) = PrMAX .* Pnow(:, i)';  % Max MU force capacity at the current time

        % UNVECTORIZED CODE
        % for mu = 1:nu  
        
        %     MU firing rate adaptation modified from Revill & Fuglevand (2011)
        %     this was modified to directly calculate the firing rate adaption, as 1 unit change in excitation causes 1 unit change in firing rate
        %     scaled to the mu threshold (highest adaptation for hightest threshold mu)

        %     % Update our duration if on
        %     if muON(mu) > 0
        %         Rdur(mu) = (i - muON(mu) + 1)/samprate;                             % duration since mu was recruited at muON
        %     end
        %     % Zero out durations if we go below zero for some reason
        %     if Rdur(mu) < 0     
        %         Rdur(mu) = 0;
        %     end      
            
        %     adaptFR(mu,i) = (thr(mu)-1)/(thr(nu)-1) * adaptSF * (mufr(mu,act) - minfr + 2) * (1 - exp(-1 * Rdur(mu) / tau )); 
        %     if adaptFR(mu,i) < 0                                                % firing rate adaptation
        %         adaptFR(mu,i) = 0;
        %     end

        %     mufrFAT(mu,i) = mufr(mu,act) - adaptFR(mu,i);                           % adapted motor unit firing rate: based on time since recruitment
        %     mufrMAX = mufr(mu,maxact) - adaptFR(mu,i);                              % adapted FR at max excitation

        %     ctFAT(mu,i) = ct(mu) * (1 + ctSF * (1 - Pnow(mu,i)/P(mu)));             % corrected contraction time: based on MU fatigue
        %     ctREL(mu,i) = ctFAT(mu,i)/ct(mu);

        %     nmufrFAT(mu,i) = ctFAT(mu,i) * (mufrFAT(mu, i) / 1000);                 % adapted normalized Stimulus Rate (CT * FR)
        %     nmufrMAX = ctFAT(mu,i) * (mufrMAX / 1000);                              % normalized FR at max excitation

        %     if nmufrFAT(mu,i) <=0.4                                               % fusion level at adapted firing rate 
        %         PrFAT(mu,i) = nmufrFAT(mu,i) / 0.4 * sPr;                         %   linear portion of curve
        %     end
        %     if nmufrFAT(mu,i) > 0.4                                                 %   non-linear portion of curve
        %           PrFAT(mu,i) = 1 - exp(-2 * (nmufrFAT(mu,i)^3));
        %     end

        %     muPt(mu, i) = PrFAT(mu, i) * Pnow(mu, i);                               % MU force at the current time (muPt): based on adapted postion on fusion curve
            
        %     if nmufrMAX <=0.4                                                   % fusion force at 100% maximum excitation
        %         PrMAX = nmufrMAX / 0.4 * sPr;                                  
        %     end
        %     if nmufrMAX > 0.4
        %         PrMAX = 1 - exp(-2 * (nmufrMAX^3));
        %     end
        %     muPtMAX(mu, i) = PrMAX * Pnow(mu, i);                               % Max MU force capacity at the current time 
        % end % next motor unit (mu)
        
        TPt(i) = sum(muPt(:,i))/maxP;                                               % total sum of MU forces at the current time (TPt)
        TPtMAX(i) = sum(muPtMAX(:,i))/maxP;   

        % used to speed up the search for the right excitation to meet the current target
        if TPt(i) < excitations(i) && act == maxact
            break;
        end
        if TPt(i) < excitations(i)
            act = act + acthop;
            if act > maxact
                act = maxact;
            end
        end
        if TPt(i) >= excitations(i) && acthop == 1     
            break;  % stop searching as the correct excitation is found
        end            
        if TPt(i) >= excitations(i) && acthop > 1
            act = act - (acthop - 1);  % if the last large jump was too much, it goes back and starts increasing by 1
            if act < 1
                act = 1;
            end
            acthop = 1;
        end
            
    end % next excitation (act)
    
    % We bypass this in the python code by checking firing rates against thresholds directly
    just_recruited = muON == 0 & ((act / res) >= thr);
    muON(just_recruited) = i;

    % UNVECTORIZED CODE
    % for mu = 1:nu
    %     if muON(mu) == 0 && (act/res) >= thr(mu)         % can be modified to reset if the MU turns off
    %         muON(mu) = i;                                       % time of onset of mu recruitment (s)
    %     end                                                     
    % end

    Ptarget(i) = TPt(i);        % modeled force level ?? do I need to do this, or can I just use TPt(i)
    Tact(i) = act;              % descending (not adapted) excitation required to meet the target force at the current time
    
    % Calculating the fatigue (force loss) for each motor unit

    % Note this assumes that mufrFAT is >= 0 for all (see orig code below)
    Pchange(:, i) = -1 * (fatigue / samprate) .* PrFAT(:, i)'; % fatigue - nominal fatigue rates

    if i < 2
        Pnow(:, i+1) = P;
    elseif i >= 2
        Pnow(:, i+1) = Pnow(:, i) + Pchange(:, i);
    end

    above_peak = Pnow(:, i+1) >= P'; % Can't happen without recovery
    Pnow(above_peak, i+1) = P(above_peak);
    below_zero = Pnow(:, i+1) < 0;
    Pnow(below_zero, i+1) = 0;  % does not let it fatigue below zero

    % UNVECTORIZED CODE
    % for mu = 1:nu
    %     if  mufrFAT(mu, i) >= 0                                                 % Force loss of each MU for each interval
    %         Pchange(mu,i) = -1 * (fatigue(mu) / samprate) * PrFAT(mu, i);       % based on % of MU fusion force
    %     % NOTE: This line never gets called in this version. recminfr is not defined.
    %     elseif mufrFAT(mu, i) < recminfr                                        
    %         Pchange(mu,i) = recovery(mu) / samprate;   
    %     end

    %     if i < 2
    %         Pnow(mu, i+1) = P(mu);
    %            % Pnow(mu, i+1) = 0;  % Use this to start the muscle fully exhausted
    %     elseif i >=2
    %         Pnow(mu, i+1) = Pnow(mu, i) + Pchange(mu,i);                            % instantaneous strength of MU
    %     end                                                                         % right now without adaptation

    %     if Pnow(mu, i+1) >= P(mu)
    %        Pnow(mu, i+1) = P(mu);                                                   % does not let it increase past rested strength
    %     end

    %     if Pnow(mu, i+1) < 0
    %        Pnow(mu, i+1) = 0;                                                       % does not let it fatigue below zero
           
    %     end
        
    % end % next motor unit
      
end % next fthsamp

Tstrength = zeros(1, fthsamp);  
for i = 1:fthsamp
    for mu = 1:nu
        muPna(mu,i) = Pnow(mu,i) * muP(mu,maxact) / P(mu);                          % non-adapted MU max force at 100% excitation (muPna)            
    end
    Tstrength(i) = sum(muPna(:,i)) / maxP;                                          % Current total strength without adaptation relative to max rested capacity
end
for i = 1:fthsamp                                   
    endurtime = i / samprate; 
    if TPtMAX(i) < excitations(i)
        break;            
    end        
end
    
                
%% Output
                
EndStrength = (TPtMAX(fthsamp) * 100);

endurtime
            
for mu = 0:nu
    if mu == 0
        mu = 1;
    end
    muForceCapacityRel(mu,ns) = Pnow(mu,ns)*100/P(mu);  % for outputs below
end

hold off;
   combo = [ns(:)/samprate, excitations(:), Tact(:)/res/maxex * 100,Tstrength(:) * 100 ,Ptarget(:) * 100,TPtMAX(:)* 100];
   dlmwrite(strcat(con,' A - Target - Act - Strength (no adapt) - Force - Strength (w adapt).csv'), combo)

   dlmwrite(strcat(con,' B - Firing Rate.csv'),transpose(mufrFAT(:,:)))
   dlmwrite(strcat(con,' C - Individual MU Force Time-History.csv'), transpose(muPt(:,:)))
   dlmwrite(strcat(con,' D - MU Capacity - relative.csv'),transpose(muForceCapacityRel(:,:)))
 
beep;

printf('Total cpu time: %f seconds\n', cputime-start_time);
  


