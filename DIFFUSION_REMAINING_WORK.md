# Что ещё не реализовано в diffusion backend

Крупные незавершённые блоки из текущих планов DVBFixer. Документ предназначен
для рабочего использования и не заменяет детальные планы в `docs/plans/`.

## 1. Генерация N/C-концов

**Статус: реализовано для Protpardelle 2026-10-06.** One-anchor contract,
детерминированная нумерация, materialization, validation, localized refinement и
provenance поддерживают N/C regions. Одновременный N+internal+C 500-step MPS
smoke прошёл все gates. Другие profiles остаются fail closed; расширенный frozen
hardware cohort и сложные multi-chain numbering fixtures ещё нужны.

Сейчас реализован только `--no-terminal`: target обрезается снаружи первого и
последнего наблюдаемого остатка, после чего моделируются внутренние gap’ы.

Остаётся для promotion:

- frozen MPS/CUDA cohort для N-only, C-only и N+C случаев;
- multi-chain tight-numbering и insertion-code fixtures;
- принятие one-anchor protocol для CUDA/Protenix.

Детальный план находится в
`docs/plans/diffusion-cli-integration.md`, раздел **A.8.1**.

## 2. Настоящее межцепочечное conditioning

**Статус: локальный protein-partner context реализован 2026-10-06.** Для каждой
gap-bearing chain Protpardelle получает отдельные contiguous crops соседних
protein chains, если fixed heavy atoms находятся в пределах 12 Å от anchor.
Partner residues входят в реальный multi-chain sampler axis как fixed motif;
trace записывает точное отображение target chain → partner chains/atoms. Общий
лимит axis остаётся 512 residues. Gap-bearing chains по-прежнему запускаются
независимо — это не whole-complex joint sampling.

Сейчас Protpardelle:

- отдельно запускается для каждой gap-bearing chain;
- затем результаты объединяются;
- остальные цепи сохраняются;
- полная структура используется при OpenMM refinement и clash validation.

Остаётся:

- совместное представление нескольких цепей без превышения лимита 512;
- адаптивный spatial crop относительно предсказанной траектории gap, а не только
  исходных anchor coordinates;
- whole-complex joint sampling нескольких gap-bearing chains;
- отдельный frozen interface benchmark и causal ablation partner-on/off;
- специализированные межцепочечные contact features сверх нативного
  Protpardelle multi-chain representation.

Это особенно важно для интерфейсов антител, олигомеров и protein–protein
contacts.

## 3. Heterogen-aware diffusion

**Статус: узкий isolated-SMILES geometric-context slice реализован 2026-10-06.**
Protpardelle получает fixed HETATM coordinates в differentiable repulsion guidance;
request/trace/provenance фиксируют exact atom mapping и authoritative graph digest.
500-step MPS smoke прошёл validation и сохранил HETATM coordinates точно. Это не
general learned chemical interaction model и не heterogen generation.

Сейчас есть только явное `--strip-heterogens`, создающее приватный protein-only
input. Это не ligand conditioning.

Остаётся:

- learned graph/charge/bond-order interaction features вместо geometric exclusion;
- metals и coordination state;
- glycans, PTMs, cofactors и waters;
- covalent protein–heterogen links;
- heterogen-aware refinement и clash/contact validation;
- расширенный benchmark strata для малых лигандов, металлов, glycans, covalent
  ligands, cofactors и waters.

Основной документ:
`docs/plans/diffusion-heterogen-gap-reconstruction.md`.

## 4. Generated-only diffusion state и per-step fixation

Текущий Protpardelle sampler обрабатывает выбранную protein chain, используя
motif conditioning. После sampling DVBFixer точно восстанавливает fixed
coordinates.

Не реализовано полностью:

- stochastic state только для generated residues;
- полная target topology до sampling;
- обязательная проекция fixed atoms после каждого denoising step в production
  MPS profile;
- доказательство identity mapping на каждом шаге;
- per-step fixed-coordinate error trace для production Apple runner.

Apple production profile использует template/motif conditioning и окончательное
exact restoration. Это отражается в trace, но слабее исходной архитектурной
цели.

## 5. Несколько candidates и seeds

CLI принимает повторяемый `--diffusion-seed`, но production Protpardelle adapter
пока требует ровно один seed и один candidate.

Не реализовано:

- последовательный запуск нескольких seeds;
- независимые traces для каждого кандидата;
- refinement каждого кандидата;
- ranking только прошедших candidates;
- публикация смеси accepted/rejected candidates;
- top-N selection.

## 6. Более длинные и сложные gap’ы

Текущий заявленный диапазон — 3–12 residues.

Не реализовано:

- отдельный research stratum 13–25 residues;
- adjacent или overlapping gaps;
- gap’ы без достаточного numbering space;
- структурно неоднозначные sequence placements;
- большие disordered linkers;
- несколько gap’ов, расстояние между которыми требует совместного sampling;
- adaptive local context/crop;
- gap-specific candidate count и sampling budget.

Несколько обычных внутренних gap’ов уже поддерживаются, но моделируются в
рамках независимых chain invocation’ов.

## 7. Полная геометрическая validation

Уже проверяются:

- atom identities;
- fixed coordinates;
- canonical heavy-atom completeness;
- peptide junctions;
- backbone breaks;
- bond lengths и angles;
- peptide planarity;
- общая Ramachandran quality;
- pooled χ1/χ2;
- clashes;
- Cα chirality.

Остаётся:

- отдельные GLY, PRO и pre-PRO Ramachandran distributions;
- χ1-only residues;
- χ3–χ5;
- residue-specific rotamer libraries;
- backbone-dependent rotamers;
- более полная side-chain geometry;
- явно записываемый `fix_ca_chirality` recovery path вместо простого rejection;
- проверка geometry после пользовательской full-system minimization.

## 8. Refinement failure как rejected output

**Статус: реализовано 2026-10-06.** Ожидаемый сбой localized OpenMM refinement
сохраняет неизменённый raw sampler candidate и публикует его только как
`validation_failed` bundle с hard gate
`localized-openmm-boundary-refinement-failed`. Ошибки artifact integrity,
containment и неожиданные исключения остаются фатальными.

Реализованное поведение:

- raw sampler candidate сохраняется без изменения координат и digest;
- стадия ошибки записывается в contract/provenance без внутренних путей;
- independent validation принудительно отклоняет такой candidate;
- CLI рекомендует другой seed или full-system minimization;
- protocol, digest, containment и неожиданные runtime errors не маскируются под
  scientific rejection.

## 9. Linux/NVIDIA и Protenix production acceptance

Не завершено:

- финальный lock Linux environment;
- NVIDIA/CUDA-compatible image;
- immutable container digest;
- checkpoint-backed production smoke;
- frozen 231-case acceptance;
- hard-gate pass-rate comparison;
- VRAM/RAM profiling;
- production launcher acceptance;
- self-hosted/manual CUDA CI lane;
- сравнение Protenix с Apple Protpardelle на одинаковых strata.

Linux profile зарегистрирован архитектурно, но не имеет завершённой production
acceptance.

## 10. Жёсткие resource limits

Уже ограничены timeout, stdout/stderr, artifact sizes, workspace и
process-group cleanup.

Не реализовано полностью:

- hard RAM limit;
- hard VRAM limit;
- cgroup/container enforcement;
- GPU admission по доступной VRAM;
- multi-process GPU coordination;
- защита от общей memory pressure на macOS;
- scheduler для нескольких hosts.

## 11. Batch, GUI и managed jobs

Diffusion доступен в single-input CLI и в Homology Model workflow GUI. Generated
GUI schema знает backend options; Homology API запускает процесс через общий
bounded runner и регистрирует bundle artifacts. Общего Model-panel workflow и
полноценного bundle viewer пока нет.

Не реализовано:

- batch directory output вида `<stem>_model_diffusion/`;
- backend-dependent output mode в command registry;
- общий `model --backend diffusion` workflow в GUI вне Homology;
- durable managed-job retry/resume для diffusion;
- специализированное отображение bundle/candidates/provenance в GUI;
- retry с другим seed из GUI;
- просмотр validation failures;
- архивирование directory bundles;
- backend capability discovery на GUI стороне.

## 12. Homology/mosaic diffusion

**Статус: начальный mosaic-first slice реализован 2026-10-06.**
`homology --backend diffusion` и GUI Model backend используют authoritative
`selected_template_mosaic.pdb`; companion coverage сохраняет zero-based
half-open masks, precedence и template ownership. Covered matching atoms fixed,
uncovered regions generated через общий contract/validation/publication path.
One-chain five-residue internal insertion прошёл 500-step MPS acceptance без
fixed-coordinate drift. MODELLER остаётся default, fallback отсутствует.

Реализовано:

- использовать `selected_template_mosaic.pdb` как authoritative frame;
- передавать template coverage metadata;
- фиксировать покрытые template atoms;
- генерировать uncovered insertion/terminal regions;
- сохранять маски template ownership в provenance;
- сохранять distinct antibody H/L PDB chain IDs;
- запускать и отображать bundle через Homology GUI workflow;
- добавить отдельный `homology --backend diffusion`.

Остаётся:

- bounded-window generation для template-covered substitutions;
- frozen multi-template и multi-chain mosaic cohorts;
- causal mosaic-adherence benchmark и сравнение с MODELLER на одинаковом plan;
- richer candidate/provenance viewer и retry seed в GUI.

## 13. Scientific benchmarking и production promotion

Apple MPS прошёл descriptive cohort, но diffusion всё ещё experimental.

Остаётся:

- окончательный license audit кода, weights и auxiliary artifacts;
- training-membership/leakage review;
- immutable dataset manifests;
- сравнение с MODELLER по всем заявленным strata;
- lDDT, GDT-HA или TM-score;
- cluster-bootstrap confidence intervals;
- natural-gap shadow set;
- отдельные показатели для multi-chain split mode;
- отдельный benchmark N/C generation;
- отдельный benchmark heterogen conditioning;
- формальное решение о production promotion.

## 14. Собственная компактная gap-модель

`docs/plans/small-diffusion-apple-silicon-experiments.md` описывает отдельную
исследовательскую программу, почти полностью не реализованную:

- dataset preparation;
- sequence-cluster split;
- temporal holdout;
- обучение моделей 5M/15M/30M/50M;
- frame/torsion/coordinate architectures;
- MPS/MLX comparison;
- checkpoint provenance;
- сравнение с MODELLER и Protenix;
- публикация собственного checkpoint.

Это не обязательный шаг для текущего Protpardelle backend, а возможная
долгосрочная замена.

## 15. Долг документации планов

Некоторые checkbox’ы в старых планах уже устарели и всё ещё выглядят
незакрытыми, хотя код реализован. Например:

- public CLI dispatch;
- Apple MPS smoke;
- multi-chain internal gaps;
- generated-reference updates;
- request-builder tests;
- bounded process cleanup;
- documentation updates.

В `diffusion-cli-integration.md` также могут оставаться старые утверждения, что
validation failure не публикует output, хотя теперь публикуется
`validation_failed` bundle.

Нужен отдельный cleanup plan-документов: не новая функциональность, а
синхронизация checklist’ов с текущим состоянием.

## Рекомендуемый практический порядок

1. Refinement-failed rejected bundles.
2. Несколько seeds/candidates.
3. N/C-terminal generation.
4. Локальное межцепочечное conditioning.
5. Linux/Protenix production acceptance.
6. Heterogen-aware conditioning.
7. Batch/GUI.
8. Homology mosaic diffusion.
