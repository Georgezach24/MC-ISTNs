function C = thesisColors()
%THESISCOLORS Κοινή παλέτα (διακρίνεται και με αχρωματοψία).
hex = @(h) sscanf(h(2:end), '%2x%2x%2x', [1 3]) / 255;
C = struct();
C.cat     = [hex('#2a78d6'); hex('#eb6834'); hex('#1baf7a')];
C.tn      = C.cat(1,:);           % επίγειο (TN)
C.ntn     = C.cat(2,:);           % δορυφορικό (NTN)
C.third   = C.cat(3,:);
C.neutral = hex('#8a8984');       % outage / ουδέτερο
C.text    = hex('#1f1f1f');
C.muted   = hex('#52514e');
C.band    = hex('#b7d3f6');       % ζώνη εκατοστημορίων
end
