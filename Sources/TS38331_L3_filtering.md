# Φιλτράρισμα μέτρησης πριν τον A3 — αποσπάσματα πηγών

Αντιγραμμένα αυτολεξεί από τα πρωτότυπα, για την παρατήρηση #1 (ρυθμός
μεταπομπών αντί για άνω φράγμα). Τα έγγραφα κατέβηκαν από το 3GPP
(το ETSI μπλοκάρει τις αυτοματοποιημένες λήψεις).

- `38331-g10.docx` → TS 38.331 **V16.1.0** (Rel-16, 2020-07)
  <https://www.3gpp.org/ftp/Specs/archive/38_series/38.331/38331-g10.zip>
- `38133-gt0.docx` → TS 38.133 **V16.29.0** (Rel-16)
  <https://www.3gpp.org/ftp/Specs/archive/38_series/38.133/38133-gt0.zip>

---

## 1. TS 38.331 §5.5.3.2 — Layer 3 filtering

Στο PDF του ETSI η ενότητα είναι στη σελίδα 141.

> The UE shall:
>
> **1>** for each cell measurement quantity, each beam measurement quantity,
> each sidelink measurement quantity as needed in sub-clause 5.8.10, and for
> each CLI measurement quantity that the UE performs measurements according to
> 5.5.3.1:
>
> **2>** filter the measured result, before using for evaluation of reporting
> criteria or for measurement reporting, by the following formula:
>
>     Fn = (1 – a)*Fn-1 + a*Mn
>
> where
>
> *Mn* is the latest received measurement result from the physical layer;
>
> *Fn* is the updated filtered measurement result, that is used for evaluation
> of reporting criteria or for measurement reporting;
>
> *Fn-1* is the old filtered measurement result, where *F0* is set to *M1* when
> the first measurement result from the physical layer is received; and for
> MeasObjectNR, a = 1/2^(ki/4), where *ki* is the filterCoefficient for the
> corresponding measurement quantity of the i:th QuantityConfigNR in
> quantityConfigNR-List, and i is indicated by quantityConfigIndex in
> MeasObjectNR; for other measurements, a = 1/2^(k/4), where *k* is the
> filterCoefficient for the corresponding measurement quantity received by the
> quantityConfig; for UTRA-FDD, a = 1/2^(k/4), where *k* is the
> filterCoefficient for the corresponding measurement quantity received by
> quantityConfigUTRA-FDD in the QuantityConfig;
>
> **2>** adapt the filter such that the time characteristics of the filter are
> preserved at different input rates, observing that the filterCoefficient *k*
> assumes a sample rate equal to **X** ms; The value of X is equivalent to one
> intra-frequency L1 measurement period as defined in TS 38.133 [14] assuming
> non-DRX operation, and depends on frequency range.
>
> NOTE 1: If k is set to 0, no layer 3 filtering is applicable.
>
> NOTE 2: The filtering is performed in the same domain as used for evaluation
> of reporting criteria or for measurement reporting, i.e., logarithmic
> filtering for logarithmic measurements.
>
> NOTE 3: The filter input rate is implementation dependent, to fulfil the
> performance requirements set in TS 38.133 [14]. For further details about the
> physical layer measurements, see TS 38.133 [14].
>
> NOTE 4: For CLI-RSSI measurement, it is up to UE implementation whether to
> reset filtering upon BWP switch.

**Τι χρησιμοποιείται από αυτό**

| Σημείο | Συνέπεια για τον κώδικα |
|---|---|
| «before using for evaluation of reporting criteria» | Ο A3 (§5.5.4.4) τροφοδοτείται με το `Fn`, όχι με το `Mn`. Αυτό έλειπε. |
| NOTE 2, λογαριθμικό φιλτράρισμα | Το φίλτρο εφαρμόζεται στο SINR σε dB, όπως ήδη το έχουμε. |
| NOTE 1, k = 0 | Το «χωρίς φίλτρο» είναι επιτρεπτή ρύθμιση του ίδιου του προτύπου, άρα η παλιά συμπεριφορά παραμένει αναπαραγώγιμη. |
| «adapt the filter … at different input rates» | Με βήμα Δt ≠ X χρειάζεται μετατροπή του `a`. |
| `F0 = M1` | Αρχικοποίηση στην πρώτη μέτρηση κάθε υποψηφίου, όχι στο μηδέν. |

**Προσαρμογή ρυθμού.** Διατήρηση των χρονικών χαρακτηριστικών σημαίνει ίδια
σταθερά χρόνου, δηλαδή `(1 − a′)^(1/Δt) = (1 − a)^(1/X)`, άρα

    a′ = 1 − (1 − a)^(Δt/X),      a = 2^(−k/4)

Για Δt = 1 s και X = 200 ms:

| k | a (στα 200 ms) | σταθερά χρόνου τ = −X/ln(1−a) | a′ (στο 1 s) |
|---|---|---|---|
| 0  | — | — | χωρίς φίλτρο (NOTE 1) |
| 4 (προεπιλογή) | 0,5000 | 289 ms | 0,969 |
| 8  | 0,2500 | 695 ms | 0,762 |
| 9  | 0,2102 | 848 ms | 0,692 |
| 11 | 0,1487 | 1,24 s | 0,551 |
| 13 | 0,1051 | 1,80 s | 0,426 |
| 15 | 0,0743 | 2,59 s | 0,321 |
| 17 | 0,0526 | 3,70 s | 0,238 |
| 19 | 0,0372 | 5,28 s | 0,174 |

Με την προεπιλογή `fc4` η σταθερά χρόνου (289 ms) είναι μικρότερη από το βήμα
του 1 s, οπότε το φίλτρο L3 είναι σχεδόν διαφανές. Το βάρος το σηκώνει η μέση
τιμή της ίδιας της μέτρησης (ενότητα 3 παρακάτω).

---

## 2. TS 38.331 — επιτρεπτές τιμές του filterCoefficient

ASN.1, ενότητα 6.3.2:

```
FilterCoefficient ::= ENUMERATED {
    fc0, fc1, fc2, fc3, fc4, fc5, fc6, fc7, fc8, fc9,
    fc11, fc13, fc15, fc17, fc19, spare1, ... }
```

Από το `QuantityConfigNR`, όπου ορίζεται χωριστός συντελεστής ανά μέγεθος:

```
filterCoefficientRSRP      FilterCoefficient   DEFAULT fc4,
filterCoefficientRSRQ      FilterCoefficient   DEFAULT fc4,
filterCoefficientRS-SINR   FilterCoefficient   DEFAULT fc4
```

Το `filterCoefficientRS-SINR` είναι το σχετικό: το μέγεθος απόφασης στην
προσομοίωση είναι SINR.

---

## 3. TS 38.133 Πίν. 9.2.5.2-1 — η περίοδος μέτρησης X

Ενότητα 9.2.5.2 «Measurement period», για ενδοσυχνοτικές μετρήσεις χωρίς
διάκενα μέτρησης:

> **Table 9.2.5.2-1: Measurement period for intra-frequency measurements
> without gaps (FR1)**
>
> | DRX cycle | T_SSB_measurement_period_intra |
> |---|---|
> | No DRX | **max(200 ms**, ceil(5 × Kp) × SMTC period) × CSSF_intra |
> | DRX cycle ≤ 320 ms | max(200 ms, ceil(1,5 × 5 × Kp) × max(SMTC period, DRX cycle)) × CSSF_intra |
> | DRX cycle > 320 ms | ceil(5 × Kp) × DRX cycle × CSSF_intra |

Η φέρουσα της προσομοίωσης είναι 3,5 GHz, δηλαδή **FR1**, και δεν
μοντελοποιείται DRX. Άρα **X = 200 ms**.

*(Για FR2 ο αντίστοιχος Πίν. 9.2.5.2-2 δίνει κατώτατο όριο 400 ms. Δεν μας
αφορά.)*

**Η δεύτερη, σημαντικότερη συνέπεια.** Ο πίνακας λέει ότι το `Mn` του
§5.5.3.2 δεν είναι στιγμιαίο δείγμα: είναι το αποτέλεσμα περιόδου μέτρησης
τουλάχιστον 200 ms. Στα 3,5 GHz και με πεζό χρήστη στα 3 km/h ο χρόνος
συνοχής είναι

    λ   = c/f    = 0,0857 m
    f_D = v/λ    = 9,73 Hz
    T_c ≈ 0,423/f_D = 43 ms

δηλαδή μέσα σε μια περίοδο μέτρησης χωρούν ≈5 ανεξάρτητες πραγματώσεις των
γρήγορων διαλείψεων, και μέσα σε ένα βήμα του 1 s χωρούν ≈23. Η προσομοίωση
τροφοδοτούσε την απόφαση με **μία** από αυτές.

---

## 4. Τι άλλαξε στον κώδικα

Δύο σκέλη, και τα δύο στο `PROD/simulateScenario.m`:

1. **Μέση τιμή μέτρησης.** Οι γρήγορες διαλείψεις υπολογίζονται ως μέση τιμή
   σε γραμμική ισχύ πάνω σε `NAvg` ανεξάρτητες πραγματώσεις ανά βήμα, αντί
   για μία. Στο δορυφορικό σκέλος η συνιστώσα σκίασης κατά Nakagami-m
   κληρώνεται **μία φορά** ανά βήμα και μέσος όρος παίρνεται μόνο στη
   σκεδαζόμενη συνιστώσα, γιατί η σκίαση δεν είναι γρήγορο φαινόμενο.
2. **Φίλτρο L3** κατά §5.5.3.2, με `a′` από τη ρήτρα προσαρμογής ρυθμού,
   εφαρμοσμένο στο SINR κάθε υποψηφίου πριν τη σύγκριση A3. Η τιμή που
   καταγράφεται στο σύνολο δεδομένων και η χωρητικότητα παραμένουν οι
   πραγματικές, όχι οι φιλτραρισμένες: το φίλτρο αφορά τη μέτρηση που βλέπει
   ο μηχανισμός απόφασης, όχι την ίδια τη ζεύξη.

Με `NAvg = 1` και `k = 0` (NOTE 1) η συμπεριφορά ανάγεται ακριβώς στην
προηγούμενη.
