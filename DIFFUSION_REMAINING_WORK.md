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

Сейчас Protpardelle:

- отдельно запускается для каждой gap-bearing chain;
- затем результаты объединяются;
- остальные цепи сохраняются;
- полная структура используется при OpenMM refinement и clash validation.

Но соседние цепи не попадают в denoiser.

Не реализовано:

- совместное представление нескольких цепей без превышения лимита 512;
- локальный crop вокруг gap с partner-chain atoms;
- conditioning gap одной цепи атомами других цепей;
- проверяемое отображение chain/residue/atom identities в sampler state;
- межцепочечные attention/contact features;
- доказательство того, что partner atoms действительно влияют на denoising.

Это особенно важно для интерфейсов антител, олигомеров и protein–protein
contacts.

## 3. Heterogen-aware diffusion

Сейчас есть только явное `--strip-heterogens`, создающее приватный protein-only
input. Это не ligand conditioning.

Не реализовано:

- сохранение ligand/glycan/metal context в sampler input;
- authoritative molecular graph для heterogen;
- bond order, formal charge, stereochemistry и coordination state;
- точное engine atom/token mapping;
- ligand SMILES как источник химической истины;
- conditioning protein gap по ligand/glycan/metal atoms;
- сохранение covalent protein–heterogen links;
- heterogen-aware refinement и clash/contact validation;
- отдельные benchmark strata для малых лигандов, металлов, glycans, covalent
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

Diffusion сейчас доступен только для single-input CLI.

Не реализовано:

- batch directory output вида `<stem>_model_diffusion/`;
- backend-dependent output mode в command registry;
- diffusion в GUI schema;
- запуск через managed-job API;
- отображение bundle/candidates/provenance в GUI;
- retry с другим seed из GUI;
- просмотр validation failures;
- архивирование directory bundles;
- backend capability discovery на GUI стороне.

## 12. Homology/mosaic diffusion

Практически отдельная большая фаза, пока не реализованная.

Планируется:

- использовать `selected_template_mosaic.pdb` как authoritative frame;
- передавать template coverage metadata;
- фиксировать покрытые template atoms;
- генерировать только uncovered regions;
- сохранять маски template ownership;
- поддерживать multi-template и multi-chain mosaics;
- корректно обрабатывать antibody H/L;
- сравнивать с MODELLER на одинаковом template plan;
- добавить отдельный `homology --backend diffusion`.

Публично включать это до mosaic-adherence acceptance не планируется.

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
