function S = init_state(P)

S.k = 0;

% Queues (bits)
S.Q.urlLC_bits = 0;
S.Q.eMBB_bits  = 0;

% Serving link state
% 0 = TN, 1 = NTN
S.servingLink = 0;   % ξεκινάμε από TN


% Handover state
S.HO.candidate = -1;        % -1 σημαίνει "κανένας υποψήφιος"
S.HO.timer = 0;             % πόσα slots κρατάει η HO συνθήκη
S.HO.active = false;        % αν εκτελείται handover τώρα
S.HO.interruptTimer = 0;    % πόσα slots διακοπής απομένουν


% Filtered metrics
S.filtRateTN  = 0;
S.filtRateNTN = 0;

end



