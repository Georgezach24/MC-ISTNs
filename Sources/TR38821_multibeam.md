# Διάταξη πολλαπλών δεσμών — αποσπάσματα πηγών

Για την παρατήρηση #4. Αντιγραμμένα από το πρωτότυπο· η εξαγωγή έγινε
απευθείας από το δυαδικό `.doc`, οπότε τα κελιά των πινάκων εμφανίζονται
χωρισμένα με `|`.

- `38821-g10.doc` → TR 38.821 **V16.1.0** (Rel-16, 2021)
  <https://www.3gpp.org/ftp/Specs/archive/38_series/38.821/38821-g10.zip>

---

## 1. Πίνακας 6.1.1.1-4 — Beam layout definition for single satellite simulation

> **Beam layout definition** | Baseline: Hexagonal mapping of the beam bore
> sight directions on UV plane defined in the satellite reference frame.
> Only the 3dB beam width parameters should be used. The beam diameter and
> beam spacing values can be computed directly from the 3 dB beam width
> assumptions and should be considered as informative.
>
> **Number of beams** | Baseline: **19-beam layout** considering wrap-around
> mechanism (i.e. 18 beams surrounding the central beam and allocated on 2
> distinct "tiers")
>
> **UV plane convention** | U axis is defined as the perpendicular line to the
> satellite-earth line on the orbital plane […] The straight line being
> orthogonal to UV plane is pointing towards the Earth centre. UV coordinates
> of the nadir of the reference satellite is (0,0)
>
> **Adjacent beam spacing on UV plane** | Baseline: Adjacent beam spacing
> computation based on 3dB beam width of the satellite antenna pattern:
>
>     ABS = sqrt(3) x sin(HPBW/2 [rad])
>
> **Central beam bore sight direction definition** | Baseline:
> Case 1: Central beam center is considered at nadir point
> Case 2: Central beam boresight direction computed based on elevation […]

Και για τον μηχανισμό αναδίπλωσης:

> For FRF = 1, two additional tiers of beams are considered in the simulation
> surrounding the 19-beam layout. For FRF > 1, four additional tiers of beams
> are considered […]

## 2. Πίνακας 6.1.1.1-5 — οι σχετικές γραμμές

> **Satellite antenna pattern** | See section 6.4.1 in [2]: Bessel function
> **Satellite polarization configuration** | Circular
> **Beam layout definition** | For single satellite simulation: See Table 6.1.1.1-4
> **Frequency re-use factor** | Option 1: 1 · Option 2: 3 · Option 3: 2 if polarization re-use is enabled
> **Polarization re-use** | Option 1: Disable · Option 2: Enable
> **UE distribution** | Base-line for calibration: at least X=10 UEs per beam with uniform distribution in all the Voronoi cell area associated to each beam.
> **Handover Margin** | 0 dB
> **UE attachment** | RSRP

Και από τον Πίν. 6.1.1.1-7 (παράμετροι αξιολόγησης επίδοσης, σε αντιδιαστολή
με τη βαθμονόμηση):

> **Handover Margin** | To be reported by the companies

*Σημείωση που αφορά την ήδη υλοποιημένη υστέρηση:* το περιθώριο 0 dB είναι της
**βαθμονόμησης**. Για αξιολόγηση επίδοσης το πρότυπο δεν δίνει τιμή αλλά ζητά
να δηλωθεί. Η επιλογή των 3 dB με χρόνο επιβεβαίωσης 2560 ms είναι επομένως
συμβατή με το πλαίσιο, αρκεί να αναφέρεται — που γίνεται.

---

## 3. Τι σημαίνει αριθμητικά για το σενάριό μας

Χρησιμοποιείται το διάγραμμα του §6.4.1 του TR 38.811 με `ka = 41.98`, τιμή
που προκύπτει από το ημιεύρος δέσμης 4,4127° του Πίν. 6.1.1.1-1 (η παραγωγή
είναι στο CLAUDE.md, παρατήρηση #12).

Από τον τύπο του πίνακα, `ABS = sqrt(3)·sin(HPBW/2) = 0.066681` στο επίπεδο UV,
δηλαδή γειτονική δέσμη **3,82° εκτός άξονα** — όχι 4,41° όπως θα υπέθετε
κανείς ταυτίζοντας την απόσταση με το ημιεύρος. Σημειώνεται ότι το γινόμενο
`ka · ABS = 2.7995` είναι **ανεξάρτητο του ka**, άρα η εξασθένηση προς τη
γειτονική δέσμη είναι σταθερή ιδιότητα της διάταξης και όχι της κεραίας.

Αθροίζοντας το εξαγωνικό πλέγμα για χρήστη στον άξονα της κεντρικής δέσμης:

| | ανά δέσμη | συνολικό I/C | οροφή SINR |
|---|---|---:|---:|
| **FRF = 1** | 6 δέσμες στις 3,82° → −10,7 dB καθεμία | **−1,3 dB** | **1,3 dB** |
| **FRF = 3** | 6 δέσμες στις 6,63° → −18,0 dB καθεμία | **−9,0 dB** | **9,0 dB** |

Για σύγκριση, το σημερινό δορυφορικό SINR της προσομοίωσης έχει **διάμεσο
7,91 dB** και μέσο 7,12 dB.

**Συνέπεια.** Με FRF = 1 η δορυφορική ζεύξη καταρρέει: η οροφή των 1,3 dB είναι
κάτω από το σημερινό σημείο λειτουργίας. Με FRF = 3 η οροφή των 9,0 dB κόβει
μόνο τις καλύτερες ζεύξεις και αφήνει το μεγαλύτερο μέρος της κατανομής. Ο
συντελεστής επαναχρησιμοποίησης δεν είναι λεπτομέρεια: **είναι η παράμετρος που
καθορίζει αν το δορυφορικό σκέλος παραμένει ανταγωνιστικό**, και το πρότυπο τον
δίνει ρητά ως επιλογή και όχι ως τιμή.

Τρεις επιφυλάξεις στον παραπάνω υπολογισμό, όλες προς την ίδια κατεύθυνση
(υπερεκτίμηση της παρεμβολής):
1. Οι πλευρικοί λοβοί του ιδανικού κυκλικού ανοίγματος είναι ψηλότεροι από
   πραγματικό ανακλαστήρα με βαθμιαίο φωτισμό.
2. Υποτίθεται ότι όλες οι ομοδιαυλικές δέσμες εκπέμπουν συνεχώς σε πλήρη ισχύ
   (πλήρης απασχόληση, συνεπές με ITU-R M.2412-0 ΠΙΝ. 5 β).
3. Ο υπολογισμός γίνεται για χρήστη ακριβώς στον άξονα της δέσμης του. Επειδή
   όμως όλη η περιοχή μελέτης υποτείνει ≤0,48° (παρατήρηση #12) και υιοθετείται
   στόχευση Case 2, όλοι οι χρήστες είναι πρακτικά στον άξονα — άρα εδώ η
   παραδοχή είναι ακριβής, όχι αισιόδοξη.

## 4. Τι θα απαιτούσε η υλοποίηση

- Στόχευση Case 2 για την κεντρική δέσμη (ήδη δηλωμένη στη μεθοδολογία).
- Θέσεις των 18 περιφερειακών δεσμών από το εξαγωνικό πλέγμα με βήμα ABS.
- Γωνία εκτός άξονα ανά χρήστη προς **κάθε** δέσμη, και το διάγραμμα του
  §6.4.1 — δηλαδή ακριβώς ο υπολογισμός που η #12 απέδειξε ότι ήταν περιττός
  για **μία** δέσμη και γίνεται απαραίτητος για πολλές.
- Επιλογή FRF (1 ή 3 ή 2 με επαναχρησιμοποίηση πόλωσης) **με αιτιολόγηση**.
- Η παραδοχή ότι οι υπόλοιπες δέσμες είναι φορτωμένες, δηλωμένη όπως και η
  αντίστοιχη για τους σταθμούς βάσης.

Αλλάζει το δορυφορικό SINR, άρα όλα τα κατάντη αποτελέσματα: μπαίνει στην ίδια
επανεκτέλεση με τις υπόλοιπες εκκρεμείς αλλαγές, όχι σε δική της.
