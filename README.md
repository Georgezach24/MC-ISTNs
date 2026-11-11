# MC-ISTNs — Ολοκληρωμένα Δορυφορικά-Επίγεια Δίκτυα

Προσομοίωση σε MATLAB ενός δικτύου 5G NR (UMa) μαζί με δορυφόρους LEO, με
κινητικότητα κατά 3GPP (Event A3, RLF, MLB) και απόφαση handover με ML στον
σταθμό, που μειώνει τα ping-pong handovers.

## Σενάριο

| Στοιχείο | Τιμή | Πηγή |
|---|---|---|
| Επίγειο δίκτυο | 19 sites × 3 τομείς, ISD 500 m, 25 m, wrap-around | TR 38.901 Πίν. 7.2-1, 7.8-1 |
| Κανάλι επίγειο | UMa, χωρικά συνεπές LOS (soft), σκίαση, Rician/Rayleigh | TR 38.901 §7.4, §7.6.3.3 |
| Δορυφόροι | 3 × LEO-600, S-band, μάσκα 20°, FRF 1 | TR 38.821 Πίν. 6.1.1.1-1 |
| Κανάλι δορυφορικό | FSPL + αέρια, Shadowed Rician | TR 38.811 §6.6.4, Abdi et al. 2003 |
| Χρήστες | 10–30 ανά διέλευση, 3 km/h, τυχαία κατεύθυνση | ITU-R M.2412-0 §8.4 |
| Βήμα | 0.2 s (περίοδος L1→L3) | TR 36.839 Πίν. 5.3.2.1, TS 38.133 |
| Μέτρηση | μέσος όρος 200 ms + φίλτρο L3 (K = 4) | TS 38.331 §5.5.3.2 |
| Handover | Event A3, Hys 3 dB, TTT 480 ms | TS 38.331 §5.5.4.4, TR 36.839 Set 1 |
| RLF | Qout −8 dB, Qin −6 dB, T310 1 s | TS 38.331 §5.3.10, TR 36.839 |
| Ping-pong | επιστροφή με παραμονή < MTS 1 s | TR 36.839 Ορισμός 5 |
| MLB | Ocn = 6 dB × διαφορά εικονικού φόρτου | TS 38.300 §15.5.1, Nasri & Altman, Awada et al. |
| Χωρητικότητα | Shannon με όριο MCS 28, ίσος διαμοιρασμός | TS 38.214 |
| Ενέργεια | EARTH (BS), γραμμικός ενισχυτής (δορυφόρος) | Auer et al. 2011 |
| Διάρκεια | 300 διελεύσεις × 900 s, 95% CI | ITU-R M.2412-0 §7.1 |

## Πολιτικές

Όλες τρέχουν στην ίδια πραγματοποίηση καναλιού, οπότε οι συγκρίσεις γίνονται
κατά ζεύγη ανά διέλευση.

| Πολιτική | Περιγραφή |
|---|---|
| `Actual` | επίγειοι και δορυφορικοί υποψήφιοι, A3 |
| `TerrestrialOnly` / `SatelliteOnly` | μόνο ένα σκέλος |
| `ActualMLB` / `TerrestrialOnlyMLB` | A3 με εξισορρόπηση φόρτου |
| `ActualML` / `ActualMLBML` | ίδιες με απόφαση ML: αναμονή αν P(ping-pong) > κατώφλι |

## Απόφαση handover με ML

- **Είσοδοι.** Μόνο ό,τι ξέρει ο σταθμός τη στιγμή της αναφοράς A3
  (`handoverFeatures.m`): μετρήσεις L3 και τάσεις τους (Aggarwal et al.),
  ιστορικό παραμονής (TS 38.413), φόρτος (TS 38.423) και γεωμετρία.
- **Ετικέτα.** Ping-pong κατά TR 36.839 Ορισμός 5. Χρησιμοποιείται μόνο ως
  ετικέτα, ποτέ ως είσοδος.
- **Εκπαίδευση.** Random Forest σε ξεχωριστούς σπόρους (seed 1000), με CV ανά
  διέλευση· κατώφλι μέγιστου F1 στις out-of-fold προβλέψεις.
- **Αξιολόγηση.** Σε κλειστό βρόχο (seed 42).
- **Αποτελέσματα.** Συμπεράσματα και αριθμοί στο [THESIS_NOTES.md](THESIS_NOTES.md).

## Εκτέλεση

Απαιτήσεις:

- MATLAB R2026a ή νεότερο (δοκιμάστηκε σε R2026a και R2026b) με 5G, Mapping, Phased Array, Communications και Aerospace Toolbox (`ecef2lla`).
- Parallel Computing Toolbox (προαιρετικό).
- Python 3 με `Model/requirements.txt`.

```matlab
addpath('PROD');
% 1) Δεδομένα εκπαίδευσης (χωρίς ML)
d = fullfile(pwd,'Dataset','ml_train','seed1000');
runSimulation(struct('rngSeed',1000,'numPasses',150,'mlGatePath','', ...
    'outputDir',d,'datasetPath',fullfile(d,'dataset.csv'), ...
    'handoverEventsPath',fullfile(d,'handover_events.mat'),'label','mltrain_seed1000'));
% 2) Εκπαίδευση
system('python Model/train_handover_gate.py');
% 3) Αξιολόγηση όλων των πολιτικών
runSimulation();
% Επαληθεύσεις
linkBudgetValidation(); geometryValidation(); energyModelValidation();
```

- **Διάρκεια** σε laptop με 6 workers: ~1.5 h για το βήμα 1, ~2 λεπτά για το
  βήμα 2, ~4 h για το βήμα 3.
- **Διακοπή.** Μια διακοπείσα εκτέλεση συνεχίζει από το σημείο ελέγχου.
- **Αναπαραγωγιμότητα.** Κάθε διέλευση έχει δικό της σπόρο
  (`rngSeed + passIdx`), οπότε η παράλληλη εκτέλεση δίνει τα ίδια αρχεία με
  τη σειριακή.

## Δομή

```
PROD/
  runSimulation.m          κύρια εκτέλεση: διελεύσεις, πολιτικές, στατιστική, γραφήματα
  simulateScenario.m       ένα βήμα: κανάλι, SINR, A3/RLF/MLB/ML, κατανομή πόρων, ενέργεια
  handoverFeatures.m       χαρακτηριστικά της απόφασης handover (κοινά για εκπαίδευση και εφαρμογή)
  handoverGatePredict.m    πρόβλεψη του Random Forest μέσα στο MATLAB
  bsSectorGainDb.m         κέρδος κεραίας τομέα (TR 38.901 Πίν. 7.3-1)
  plotNetworkSnapshot.m    γράφημα τοπολογίας και ζεύξεων
  saveRunVersion.m         αρχειοθέτηση εκτέλεσης στο Results/runs/
  *Validation.m            επαληθεύσεις: ισοζύγιο ζεύξης, γεωμετρία, ενεργειακό μοντέλο
Model/
  train_handover_gate.py   εκπαίδευση της απόφασης handover
  results_handover_gate/   μοντέλο, metrics.json, γραφήματα
Results/                   γραφήματα και simulation_results.mat της τελικής εκτέλεσης
Results/runs/              παράμετροι και έξοδοι κάθε εκτέλεσης
Dataset/                   dataset.csv, handover_events.mat (παράγονται, εκτός git)
```
