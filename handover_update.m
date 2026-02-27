function [S, hoTN, hoNTN] = handover_update(S, P, L)
% Συνάρτηση handover (HO) ελεγχόμενου από το δίκτυο για:
%   - TN (Terrestrial BS handover)
%   - NTN (Satellite handover)
%
% Κύρια χαρακτηριστικά αλγορίθμου:
%   1) Φιλτραρισμένο SINR (EWMA) για αποφυγή αστάθειας λόγω fast fading
%   2) Hysteresis margin (σε linear) για αποφυγή ping-pong
%   3) Time-To-Trigger (TTT): η συνθήκη πρέπει να ισχύει συνεχόμενα για κάποιο χρόνο
%   4) Execution timer: περίοδος interruption όπου το link θεωρείται μη διαθέσιμο
%
% Έξοδοι:
%   - S: ενημερωμένο state (serving node, HO timers/flags)
%   - hoTN: true αν ολοκληρώθηκε TN handover στο slot
%   - hoNTN: true αν ολοκληρώθηκε NTN handover στο slot

% =====================================================================
%                     Αρχικοποίηση HO event flags
% =====================================================================
hoTN  = false;
hoNTN = false;
%----------------------------------------------------------------------

% =====================================================================
% 0) Update filtered metrics (RSRP/RSRQ-like filtering)
% =====================================================================
% Χρησιμοποιούμε EWMA φίλτρο:
%   filt = alpha * filt_prev + (1-alpha) * instant
% όπου alpha κοντά στο 1 => ισχυρή εξομάλυνση (λιγότερη ευαισθησία στο fading)

S.HO_TN.filtSINR  = P.HO.alpha * S.HO_TN.filtSINR  + (1 - P.HO.alpha) * L.TN.SINR_DL;
S.HO_NTN.filtSINR = P.HO.alpha * S.HO_NTN.filtSINR + (1 - P.HO.alpha) * L.NTN.SINR_DL;
%----------------------------------------------------------------------

% --------- Φιλτραρισμένες μετρικές για αποφάσεις HO -------------------
TNm  = S.HO_TN.filtSINR;   % Filtered SINR προς κάθε BS candidate
NTNm = S.HO_NTN.filtSINR;  % Filtered SINR προς κάθε SAT candidate
%----------------------------------------------------------------------

% =====================================================================
% 1) TN HO (BS selection & execution)
% =====================================================================
if ~S.HO_TN.active
    % --------------------------------------------------------------
    % Decision phase (δεν είμαστε σε HO execution)
    % --------------------------------------------------------------

    % --------- Επιλογή καλύτερου candidate BS ----------------------
    [bestMetric, bestBS] = max(TNm);     % Καλύτερος BS βάσει filtered SINR
    servMetric = TNm(S.servBS);         % Metric του τρέχοντος serving BS
    %---------------------------------------------------------------

    % --------- Έλεγχος HO condition (hysteresis + διαφορετικός BS) --
    % HO trigger αν:
    %   - bestBS != serving
    %   - bestMetric > servingMetric * hyst (σε linear κλίμακα)
    if bestBS ~= S.servBS && bestMetric > servMetric * P.HO.hyst

        % Μετράμε TTT όσο η συνθήκη ισχύει συνεχόμενα
        S.HO_TN.ttt_s = S.HO_TN.ttt_s + P.dt;
        S.HO_TN.cand  = bestBS;

        % Αν συμπληρωθεί το TTT, ξεκινά HO execution
        if S.HO_TN.ttt_s >= P.HO.TTT_TN_s
            S.HO_TN.active  = true;           % Μπαίνουμε σε interruption phase
            S.HO_TN.timer_s = P.HO.exec_TN_s; % Χρόνος execution (κομμένο link)
            S.HO_TN.ttt_s   = 0;              % Reset TTT counter
        end

    else
        % Αν η συνθήκη δεν ισχύει συνεχόμενα:
        % - Reset TTT (αποφυγή ψευδο-trigger)
        % - Candidate επιστρέφει στο serving
        S.HO_TN.ttt_s = 0;
        S.HO_TN.cand  = S.servBS;
    end

else
    % --------------------------------------------------------------
    % Execution phase (interruption)
    % --------------------------------------------------------------
    % Κατά την execution phase, θεωρούμε ότι το TN link δεν μεταφέρει δεδομένα
    % (αυτό εφαρμόζεται στην apply_action_and_serve με capTN=0).
    S.HO_TN.timer_s = S.HO_TN.timer_s - P.dt;

    % Όταν τελειώσει ο execution timer, ολοκληρώνεται το HO:
    if S.HO_TN.timer_s <= 0
        S.servBS = S.HO_TN.cand;   % Νέος serving BS
        S.HO_TN.active  = false;   % Έξοδος από interruption
        S.HO_TN.timer_s = 0;       % Reset timer
        hoTN = true;               % Καταγραφή HO event στο slot
    end
end
%----------------------------------------------------------------------

% =====================================================================
% 2) NTN HO (Satellite selection & execution)
% =====================================================================
if ~S.HO_NTN.active
    % --------------------------------------------------------------
    % Decision phase για δορυφόρους
    % --------------------------------------------------------------

    % --------- Επιλογή καλύτερου candidate SAT ----------------------
    [bestMetric, bestSAT] = max(NTNm);   % Καλύτερος SAT βάσει filtered SINR
    servMetric = NTNm(S.servSAT);       % Metric του τρέχοντος serving SAT
    %---------------------------------------------------------------

    % --------- Έλεγχος HO condition (hysteresis + διαφορετικός SAT) --
    if bestSAT ~= S.servSAT && bestMetric > servMetric * P.HO.hyst

        % Μετράμε TTT όσο η συνθήκη ισχύει συνεχόμενα
        S.HO_NTN.ttt_s = S.HO_NTN.ttt_s + P.dt;
        S.HO_NTN.cand  = bestSAT;

        % Αν συμπληρωθεί το TTT, ξεκινά HO execution
        if S.HO_NTN.ttt_s >= P.HO.TTT_NTN_s
            S.HO_NTN.active  = true;            % interruption phase
            S.HO_NTN.timer_s = P.HO.exec_NTN_s; % execution time
            S.HO_NTN.ttt_s   = 0;               % reset TTT
        end

    else
        % Reset TTT αν η συνθήκη δεν κρατάει συνεχόμενα
        S.HO_NTN.ttt_s = 0;
        S.HO_NTN.cand  = S.servSAT;
    end

else
    % --------------------------------------------------------------
    % Execution phase (interruption) για NTN
    % --------------------------------------------------------------
    S.HO_NTN.timer_s = S.HO_NTN.timer_s - P.dt;

    % Ολοκλήρωση NTN HO όταν μηδενίσει ο timer
    if S.HO_NTN.timer_s <= 0
        S.servSAT = S.HO_NTN.cand;  % Νέος serving SAT
        S.HO_NTN.active  = false;
        S.HO_NTN.timer_s = 0;
        hoNTN = true;              % Καταγραφή HO event στο slot
    end
end
%----------------------------------------------------------------------

end