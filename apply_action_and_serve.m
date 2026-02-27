function [S, stepKPI] = apply_action_and_serve(S, P, L, action)
% Συνάρτηση εξυπηρέτησης (serving) των υπηρεσιών URLLC και eMBB στο τρέχον slot,
% με βάση την επιλεγμένη ενέργεια (action) και λαμβάνοντας υπόψη HO interruptions.
%
% Modes:
%   URLLC modeU: 0 TN, 1 NTN, 2 SPLIT, 3 DUP
%   eMBB  modeE: 0 TN, 1 NTN, 2 SPLIT (aggregation)
%
% Έξοδοι:
%   - S: ενημερωμένο state (queue flags, counters, head pointer)
%   - stepKPI: KPI του slot (URLLC delivered/success/delay, eMBB bits)

% =====================================================================
%                         Αρχικοποίηση KPI slot
% =====================================================================
stepKPI = struct('URLLC_delivered',0,'URLLC_delay',NaN,'URLLC_success',0,...
                 'eMBB_bitsDL',0,'eMBB_bitsUL',0);
%----------------------------------------------------------------------

% =====================================================================
%                 Υπολογισμός χωρητικοτήτων ανά slot (bits/slot)
% =====================================================================
% Μετατρέπουμε τους στιγμιαίους ρυθμούς (bits/sec) σε χωρητικότητα ανά slot:
capTNdl  = L.R_TN_DL  * P.dt;   % TN downlink bits/slot
capTNul  = L.R_TN_UL  * P.dt;   % TN uplink   bits/slot
capNTNdl = L.R_NTN_DL * P.dt;   % NTN downlink bits/slot
capNTNul = L.R_NTN_UL * P.dt;   % NTN uplink   bits/slot
%----------------------------------------------------------------------

% =====================================================================
%                 HO interruptions: μηδενισμός capacity
% =====================================================================
% Αν εκτελείται handover (execution phase), το link θεωρείται “κομμένο”
% και άρα η χωρητικότητα του slot είναι 0.
if S.HO_TN.active
    capTNdl = 0; 
    capTNul = 0;
end
if S.HO_NTN.active
    capNTNdl = 0; 
    capNTNul = 0;
end
%----------------------------------------------------------------------

% =====================================================================
%                         URLLC (DL, packet-level)
% =====================================================================
% Μοντέλο:
% - FIFO ουρά πακέτων
% - Εξυπηρετεί το πολύ 1 πακέτο ανά slot (head-of-line)
% - Αν δεν “χωράει” στο slot (capacity < pktBits) δεν γίνεται partial send
%----------------------------------------------------------------------

% --------- Safety: ύπαρξη head/tail (για robustness) ------------------
if ~isfield(S.URLLC,'tail')
    S.URLLC.tail = numel(S.URLLC.queue);
end
if ~isfield(S.URLLC,'head')
    S.URLLC.head = 1;
end
%----------------------------------------------------------------------

% --------- Εντοπισμός head-of-line πακέτου ---------------------------
h = S.URLLC.head;

% Προχωράμε το head μέχρι να βρούμε πακέτο που δεν είναι delivered/expired
while (h <= S.URLLC.tail) && (S.URLLC.queue(h).delivered || S.URLLC.queue(h).expired)
    h = h + 1;
end
S.URLLC.head = h;
%----------------------------------------------------------------------

% Αν υπάρχει ενεργό πακέτο, το idx δείχνει στο head-of-line
idx = [];
if h <= S.URLLC.tail
    idx = h;
end
%----------------------------------------------------------------------

% --------- Εξυπηρέτηση 1 πακέτου max ανά slot -------------------------
if ~isempty(idx)

    pktBits = S.URLLC.pktBits; % Μέγεθος URLLC πακέτου σε bits

    % Επιλογή serving mode για URLLC σύμφωνα με action.modeU
    switch action.modeU

        % ==============================================================
        % URLLC Mode 0: TN μόνο
        % ==============================================================
        case 0 % TN
            if capTNdl >= pktBits
                % Delay = αναμονή στην ουρά + propagation delay του TN link
                delay = (S.t - S.URLLC.queue(idx).genTime) + L.Dprop_TN;
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end

        % ==============================================================
        % URLLC Mode 1: NTN μόνο
        % ==============================================================
        case 1 % NTN
            if capNTNdl >= pktBits
                % Delay = αναμονή + propagation delay του NTN link
                delay = (S.t - S.URLLC.queue(idx).genTime) + L.Dprop_NTN;
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end

        % ==============================================================
        % URLLC Mode 2: SPLIT (συντηρητικό split σε 2 links)
        % ==============================================================
        % Εδώ υποθέτουμε ότι το πακέτο “σπάει” σε δύο ίσα μισά και για να
        % ολοκληρωθεί χρειάζονται ΚΑΙ τα δύο μισά μέσα στο ίδιο slot.
        case 2 % SPLIT
            if (capTNdl/2 >= pktBits/2) && (capNTNdl/2 >= pktBits/2)
                % Propagation delay λαμβάνει το χειρότερο link (max)
                delay = (S.t - S.URLLC.queue(idx).genTime) + max(L.Dprop_TN, L.Dprop_NTN);
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end

        % ==============================================================
        % URLLC Mode 3: DUP (διπλή μετάδοση)
        % ==============================================================
        % Αρκεί να μπορέσει να περάσει το πακέτο από ΕΝΑ από τα δύο links.
        % Delay θεωρείται του πρώτου που θα φτάσει (min propagation delay).
        case 3 % DUP
            if (capTNdl >= pktBits) || (capNTNdl >= pktBits)
                delay = (S.t - S.URLLC.queue(idx).genTime) + min(L.Dprop_TN, L.Dprop_NTN);
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end
    end
end
%----------------------------------------------------------------------

% =====================================================================
%               Έλεγχος λήξης deadline (FIFO expiry sweep)
% =====================================================================
% Κάνουμε αποδοτικό έλεγχο μόνο από το head προς τα εμπρός:
% - Αν το head έχει ήδη delivered/expired -> προχωράμε
% - Αν το head έχει ξεπεράσει deadline -> το μαρκάρουμε expired και fail++
% - Μόλις βρούμε πακέτο που δεν έχει λήξει, σταματάμε (FIFO ιδιότητα)
while S.URLLC.head <= S.URLLC.tail

    if S.URLLC.queue(S.URLLC.head).delivered || S.URLLC.queue(S.URLLC.head).expired
        S.URLLC.head = S.URLLC.head + 1;
        continue;
    end

    if (S.t - S.URLLC.queue(S.URLLC.head).genTime) > P.URLLC.deadline
        S.URLLC.queue(S.URLLC.head).expired = true;
        S.URLLC.failCount = S.URLLC.failCount + 1;
        S.URLLC.head = S.URLLC.head + 1;
    else
        break;
    end
end
%----------------------------------------------------------------------

% =====================================================================
%                      eMBB (DL/UL throughput)
% =====================================================================
% eMBB δεν είναι packet-level εδώ. Μετράμε bits ανά slot ανάλογα με modeE.
switch action.modeE

    % eMBB Mode 0: TN μόνο
    case 0 % TN
        stepKPI.eMBB_bitsDL = capTNdl;
        stepKPI.eMBB_bitsUL = capTNul;

    % eMBB Mode 1: NTN μόνο
    case 1 % NTN
        stepKPI.eMBB_bitsDL = capNTNdl;
        stepKPI.eMBB_bitsUL = capNTNul;

    % eMBB Mode 2: SPLIT / Aggregation
    % Υποθέτουμε ότι μπορούμε να αθροίσουμε throughput από TN+NTN
    case 2 % SPLIT aggregate
        stepKPI.eMBB_bitsDL = capTNdl + capNTNdl;
        stepKPI.eMBB_bitsUL = capTNul + capNTNul;
end
%----------------------------------------------------------------------

end


% =====================================================================
%                    Παράδοση URLLC πακέτου
% =====================================================================
function [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI)
% Ενημερώνει την ουρά URLLC όταν ένα πακέτο παραδοθεί:
% - θέτει delivered=true
% - αποθηκεύει delay
% - ενημερώνει success/fail counters ανάλογα με deadline

% --------- Ενημέρωση πακέτου στην ουρά --------------------------------
S.URLLC.queue(idx).delivered = true;
S.URLLC.queue(idx).delay     = delay;
%----------------------------------------------------------------------

% --------- Έλεγχος deadline ------------------------------------------
ok = (delay <= P.URLLC.deadline);
if ok
    S.URLLC.succCount = S.URLLC.succCount + 1;
else
    S.URLLC.failCount = S.URLLC.failCount + 1;
end
%----------------------------------------------------------------------

% --------- Ενημέρωση KPI slot ----------------------------------------
stepKPI.URLLC_delivered = 1;
stepKPI.URLLC_delay     = delay;
stepKPI.URLLC_success   = ok;
%----------------------------------------------------------------------

end