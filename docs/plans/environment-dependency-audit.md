# План B. Ревизия модулей и состава окружения

**Одобрено пользователем 2026-10-28:** реализовать Standard/Desktop и общий
установочный helper. **PySCF и RESP на его основе остаются в Standard и Desktop**,
как и MDAnalysis/Plotly, MODELLER, химический и антительный стек. Отдельного
дополнения PySCF не будет. Актуальная точка входа и порядок реализации:
[handoff профилей установки](installation-profiles-handoff.md).
Реализация начата: composition/specs и базовый helper проверены на Linux,
включая PySCF и защищённое лицензирование. Актуальные выполненные и оставшиеся
этапы: [profile acceptance](../installation-profile-acceptance.md).
Эта редакция имеет приоритет над прежними предложениями о дроблении extras.

**Обновление 2026-10-28:** исследование Probe/xTB завершено; по
[ADR 0012](../adr/0012-remove-xtb-refinement-and-repair-probe.md) xTB удалён из
поддерживаемых команд и основного окружения, интеграция опционального Probe
исправлена. Упоминания будущего xTB-профиля ниже — исходный план, не актуальный
список возможностей. Это не завершает остальные пункты ревизии окружения.

**Статус на 2026-10-20: частично подготовлен, миграция packaging не выполнена.**
Отдельные ML runner environments и `doctor` уже существуют; это не означает,
что выполнены описанные ниже optional installation profiles, полная ревизия
transitive/executable dependencies или clean-install matrix. В `pyproject.toml`
пока определён только optional extra `dev`. Исправление eager Torch import в
Protenix helper не добавляет Torch в core и не завершает весь аудит.
Сводный статус: [индекс планов](README.md).

**Дополнение на 2026-10-27:** принято направление Standard/Desktop (B.4),
проверки xTB/Probe до миграции (B.12–B.13) и helper на conda/mamba/micromamba
без uv/Poetry (B.14). Эти решения имеют приоритет над исходными предложениями
о дроблении extras и переносе MDAnalysis/Plotly ниже. Packaging ещё не изменён.

## B.1. Цель и правила аудита

- Провести отдельную ревизию:
  - `environment.yml`;
  - обязательных dependencies в `pyproject.toml`;
  - будущих optional dependency groups;
  - внешних executable;
  - GUI/dev dependencies;
  - platform-specific diffusion environments.
- Не считать отсутствие direct import доказательством ненужности.
- Для каждого пакета проверить:
  - eager и lazy imports;
  - command entrypoints;
  - subprocess executable lookups;
  - optional backend paths;
  - tests/fixtures;
  - docs/install instructions;
  - GUI/build scripts;
  - возможную transitive/runtime роль.
- Категоризировать каждую dependency:
  - core runtime;
  - optional command/backend;
  - dev/test/build;
  - external executable;
  - platform diffusion runner;
  - transitive/redundant;
  - unused pending proof.
- Не удалять пакет до:
  - negative-evidence audit;
  - clean-environment smoke;
  - focused regression tests;
  - синхронизации packaging/docs/doctor/CI.

## B.2. Зафиксированная evidence-based классификация

Ниже сохранена исходная классификация потребителей для аудита. Предложения
о пользовательских extras/core/full не являются актуальной комплектацией:
для реализации использовать одобренный состав B.4 и handoff. В частности,
PySCF и MDAnalysis/Plotly не выносятся из Standard.

- **Оставить в минимальном core:**
  - `python >=3.11,<3.14`:
    - interpreter constraint в repository policy и package metadata;
  - `numpy <2.5`:
    - базовые scientific paths и совместимость;
  - `scipy`:
    - steric diagnostics и соответствующие tests;
    - может быть вынесен только вместе с формальным переводом diagnose feature в optional extra;
  - `openmm >=8.1`:
    - preparation, minimization, topology/force-field paths;
  - `pdbfixer >=1.9`:
    - default legacy preparation backend;
  - `biopython`:
    - sequence/alignment paths;
  - `gemmi >=0.7.5`:
    - advertised mmCIF/PDBx detection/conversion на CLI boundary;
  - `pip`:
    - editable install/bootstrap flow;
  - `propka >=3.5`:
    - standalone protonation и integrated preparation variant pass.
- **Оставить в full environment, но выделить licensed modeling extra:**
  - `modeller`:
    - production-supported `model` и `homology`;
    - отдельная license setup;
    - отсутствует в обычных pip dependencies.
- **Выделить ligand/chemistry extra:**
  - `openmmforcefields`:
    - `GAFFTemplateGenerator` в ligand parameterization;
  - `openff-toolkit`:
    - direct use `Molecule` и AmberTools wrapper;
    - отсутствует в `pyproject.toml`;
  - `rdkit`:
    - glycan/heterogen hydrogen path;
    - tests уже допускают отсутствие;
  - `openbabel`/`openbabel-wheel`:
    - small-molecule CIF conversion;
    - ligand extraction;
    - heterogen hydrogen processing;
    - optional minimization utilities;
    - оставить в full environment, пока поддерживается broad-input chemistry.
- **Выделить analysis/trajectory extra:**
  - `MDAnalysis`:
    - `cluster.py`, `align.py`, trajectory/GRO workflows;
  - `plotly`:
    - optional HTML plotting;
    - код уже умеет предупреждать и пропускать plot.
- **Выделить parameterization/topology extra:**
  - `parmed`:
    - parameterization и GROMACS export;
  - `ambertools >=23`:
    - `sqm`, `antechamber`, `parmchk2`, `tleap`;
    - нужен ligand parameterization и `tleap-reduce` backend;
  - `reduce`:
    - standalone executable `reduce -build -nuclear -quiet`;
    - сохранить отдельно до platform package proof;
  - `acpype`:
    - только `top --acpype`/topology export.
- **Выделить antibody extra:**
  - `anarci`:
    - antibody numbering/CDR paths в `homology.py` и `antibody.py`;
  - `hmmer`:
    - прямой DVBFixer consumer не найден;
    - статус `unused pending proof`;
    - не удалять до проверки package metadata и isolated ANARCI runtime smoke.
- **Оставить один default MSA executable, альтернативы вынести:**
  - решение принято пользователем 2026-10-27; изменение packaging ещё не выполнено;
  - `mafft`:
    - `msa.py --auto` выбирает MAFFT первым;
    - docs указывают, что одного MAFFT достаточно;
  - `muscle >=5`:
    - supported explicit alternative;
    - version floor важен из-за CLI syntax;
  - `clustalo`:
    - supported explicit alternative.
  - ANARCI/HMMER сохранить единым антительным набором. Добавить preflight
    Python-пакета, HMMER executable и профильных данных, включая `doctor`.
    Отсутствие/ошибку инструмента отличать от отсутствия антительного домена:
    явно запрошенная антительная схема не должна незаметно заменяться SEQRES.
  - Обычный `salign --engine biopython` также требует один MSA executable;
    отсутствие MODELLER не означает отсутствие внешних зависимостей.
- **Выделить refinement/QM extras:**
  - `xtb`:
    - только explicit `--refine xtb`;
    - на карандаше: пока сохранить, до решения о профиле выполнить проверки B.12;
  - `pyscf`:
    - только RESP `--qm-engine pyscf`.
- **Перенести в GUI/dev environment:**
  - `nodejs >=22`:
    - нужен GUI build/typecheck/test;
    - не нужен Python CLI runtime.
- **Хранить только в platform runner environments:**
  - PyTorch;
  - CUDA runtime/tooling;
  - MPS-compatible engine dependencies;
  - Protenix packages/patches;
  - Protpardelle packages/patches;
  - ML caches/checkpoint configuration.

## B.3. Противоречия, которые надо устранить

- `pyproject.toml` сейчас делает обязательными command-specific dependencies:
  - `MDAnalysis`;
  - `openbabel-wheel`;
  - `parmed`;
  - `acpype`;
  - `plotly`.
- `pyproject.toml` одновременно не описывает direct optional dependencies:
  - OpenFF Toolkit;
  - RDKit;
  - MODELLER;
  - ANARCI;
  - PySCF.
- Conda route использует `openbabel`, pip route — `openbabel-wheel`:
  - определить одну поддерживаемую binding route для каждого installation type;
  - не устанавливать обе реализации в одно environment без compatibility evidence;
  - добавить smoke import/feature test для выбранной route.
- `environment.yml` comment упоминает `pdb4amber`, но consumer не найден:
  - удалить/исправить комментарий;
  - не удалять AmberTools, поскольку остальные executables доказанно используются.
- `prep_backend.py` говорит, что AmberTools предоставляет Reduce, но:
  - environment отдельно устанавливает `reduce`;
  - код ищет самостоятельный executable;
  - conda package content может различаться по платформам;
  - сохранить обе записи до metadata verification на Linux и macOS.
- Три MSA engines избыточны для default install:
  - оставить MAFFT в full/default feature set;
  - MUSCLE и Clustal Omega оформить как alternatives extra/environment.
- HMMER не имеет прямого DVBFixer consumer, но локальные conda-метаданные
  ANARCI подтверждают runtime dependency `hmmer >=3.3.2`. Решение пользователя:
  сохранить ANARCI/HMMER вместе; требуется smoke всей цепочки, а не удаление
  HMMER по отсутствию прямого вызова в DVBFixer.

## B.4. Новая модель installation profiles

### Одобренная комплектация и старт реализации (2026-10-28)

Пользователь выбрал вариант с двумя готовыми комплектациями:

- **Standard**: вся обычная научная функциональность — MODELLER/Homology,
  подготовка белков и поддерживаемых комплексов, протонирование, минимизация,
  автоматическая параметризация лигандов, топологии/ACPYPE, антительная
  нумерация, анализ траекторий с MDAnalysis/Plotly, MAFFT и **PySCF/RESP**.
- **Desktop**: Standard плюс GUI. Научный состав не дублируется отдельным
  независимым списком зависимостей.
- Отдельные дополнения: платформенные diffusion runtimes с PyTorch,
  исходниками/патчами движка, refinement, checkpoint verification и регистрацией;
  отдельный PSI4/psiresp runtime; внешняя Gaussian-интеграция;
  MUSCLE 5/Clustal Omega. Это не требование установить все дополнения сразу.
- Лицензия MODELLER и предоставление diffusion checkpoints — явные этапы
  готовности установки, не автоматически выполненные действия.
- xTB уже удалён по ADR 0012 и не входит ни в один профиль. Probe остаётся
  опциональным диагностическим компонентом с проверкой поддерживаемой версии.
- PySCF уже установлен текущим `environment.yml`; сохраняется в обоих основных
  профилях. Его перенос в отдельный optional extra не согласован и не выполняется.

**Порядок работ:** исследования xTB/Probe завершены, решение принято и включено
в `main` с релизом 0.10.1. Пользователь одобрил реализацию плана с сохранением
PySCF в Standard/Desktop. После обновления сессии продолжить в
`feature/installation-profiles` по handoff: состав/контракты → спецификации →
helper/preflight → чистые установки/smoke → документация и CI → перенос
GUI-only/альтернативных MSA компонентов после подтверждённой замены.
Не требовать повторного одобрения уже согласованной комплектации.

### Историческое предложение: не принятая матрица профилей

Следующий детальный конструктор core/full и функциональных extras сохранён
как исходный материал аудита, а не как актуальный пользовательский интерфейс.
Принятое направление Standard/Desktop выше имеет приоритет.

- Сначала определить стабильные profiles, затем менять dependencies.
- Предлагаемые profiles:
  - `core`:
    - baseline structure preparation;
    - CIF boundary support;
    - OpenMM/PDBFixer;
    - core diagnostics/protonation;
  - `full`:
    - backward-compatible conda environment со всеми production-supported non-ML features;
  - `modeling`:
    - MODELLER и license instructions;
  - `chemistry`:
    - OpenFF, openmmforcefields, RDKit, Open Babel;
  - `analysis`:
    - MDAnalysis и Plotly;
  - `parameterization`:
    - ParmEd, AmberTools, Reduce, ACPYPE;
  - `antibody`:
    - ANARCI и подтверждённые runtime requirements;
  - `msa-alternatives`:
    - MUSCLE 5 и Clustal Omega;
  - `refinement`:
    - xTB;
  - `qm`:
    - PySCF;
  - `gui-dev`:
    - Node.js/npm toolchain;
  - `diffusion-protenix-cuda`:
    - отдельный reproducible Linux runner environment;
  - `diffusion-protpardelle-mps`:
    - отдельный reproducible Apple runner environment.
- Не обязательно публиковать все extras одним релизом:
  - сначала зафиксировать matrix и user-facing feature guarantees;
  - затем мигрировать по одному coherent profile;
  - сохранять `full` как compatibility path до завершения smoke coverage.

## B.5. Изменения файлов

Это перечень затрагиваемых областей из исходного аудита. Пользовательские
комплектации и границы переноса определяются одобренным B.4 и handoff;
исторические названия core/full не вводят дополнительные продукты.

- `environment.yml`:
  - сохранить full/backward-compatible назначение либо переименовать/документировать его явно;
  - убрать GUI-only и alternative tools только после появления отдельных environment specs;
  - исправить misleading comments;
  - не добавлять ML dependencies.
- `pyproject.toml`:
  - оставить истинный core в `dependencies`;
  - добавить optional dependency groups там, где packages доступны через pip и поддержаны;
  - не обещать через pip packages с отдельной лицензией/conda-only route без корректных instructions;
  - согласовать naming extras с docs и `doctor`.
- Добавить отдельные environment specs/lock inputs:
  - GUI/dev;
  - Linux Protenix CUDA runner;
  - Apple Protpardelle MPS runner;
  - при необходимости conda feature profiles, если pip extras недостаточны для executables.
- `docs/installation.md`:
  - описать core/full/optional profiles;
  - указать feature-to-extra mapping;
  - дать actionable install instructions;
  - документировать MODELLER license;
  - документировать Open Babel route;
  - документировать external executable discovery.
- `docs/known-issues.md`:
  - platform/package conflicts;
  - Reduce/AmberTools differences;
  - Open Babel binding ambiguity;
  - MPS/CUDA constraints.
- `doctor.py`:
  - не объявлять optional packages общими failures;
  - проверять requirements выбранной feature/profile;
  - выдавать actionable missing-extra/missing-executable diagnostics.
- CI:
  - minimal core lane;
  - full environment lane;
  - focused optional-profile smoke lanes;
  - GUI lane;
  - diffusion runner acceptance остаётся hardware/manual/self-hosted.

## B.6. Порядок миграции dependencies

- Этап 1 — inventory freeze:
  - зафиксировать package, version constraint, consumer, install channel и category;
  - сверить `environment.yml`, `pyproject.toml`, imports, executables, docs и CI.
- Этап 2 — contract definition:
  - зафиксировать научные гарантии Standard, включая PySCF/RESP;
  - определить наследование Desktop и дополнительные GUI-проверки;
  - определить ошибки/отложенную готовность выбранных дополнений.
- Этап 3 — packaging alignment:
  - добавить optional groups/specs;
  - исправить direct/optional dependency declarations;
  - сохранить совместимый ручной маршрут до готовности подтверждённой замены.
- Этап 4 — clean-environment proof:
  - создать чистый Standard и проверить заявленные научные команды;
  - создать Desktop с тем же научным составом и проверить GUI;
  - создать каждое выбранное дополнение с нужной изоляцией;
  - выполнить representative command smoke, не ограничиваться imports.
- Этап 5 — slim full/default при подтверждённой политике:
  - перенести Node.js в GUI/dev;
  - перенести alternative MSA engines;
  - xTB уже удалён; PySCF остаётся в Standard/Desktop;
  - сохранить analysis/parameterization/chemistry features в Standard; изменения
    каналов установки/metadata не должны незаметно исключать их пакеты;
  - не удалять HMMER/Reduce до отдельных proof gates.
- Этап 6 — удалить obsolete declarations/comments:
  - только после успешных smoke lanes;
  - синхронно обновить docs, doctor и CI.

## B.7. Проверки для каждого optional profile

Ниже перечислены функциональные группы проверок из аудита, не дополнительные
пользовательские presets. Обычные научные группы и RESP/PySCF входят в
Standard/Desktop; отдельно устанавливаются только согласованные дополнения.

- `core`:
  - package import;
  - CLI startup/help;
  - PDB preparation legacy backend;
  - CIF boundary conversion;
  - protonation;
  - core diagnostics;
  - representative no-solvent minimization.
- `modeling`:
  - MODELLER import/license preflight;
  - `model` default path;
  - `homology` smoke.
- `chemistry`:
  - Open Babel import/executable path;
  - small-molecule CIF conversion;
  - OpenFF molecule creation;
  - GAFF template path;
  - RDKit heterogen/glycan H path.
- `analysis`:
  - trajectory/GRO input;
  - cluster/align;
  - Plotly present и absent behavior.
- `parameterization`:
  - discovery `antechamber`, `parmchk2`, `sqm`, `tleap`, `reduce`;
  - ParmEd import;
  - ACPYPE opt-in export;
  - `tleap-reduce` preflight.
- `antibody`:
  - ANARCI numbering/CDR smoke;
  - HMMER requirement experiment;
  - actionable failure без profile.
- `msa-alternatives`:
  - MAFFT auto/default;
  - MUSCLE 5 syntax;
  - Clustal Omega explicit selection;
  - stable missing-executable errors.
- `refinement`:
  - xTB уже удалён: проверять отклонение старых CLI-флагов и отсутствие API;
  - существующий OpenBabel post-pass остаётся отдельным явным запросом и не
    получает новую гарантию качества за счёт этой миграции.
- RESP/PySCF в составе Standard/Desktop:
  - PySCF RESP route;
  - no PySCF import on unrelated commands.
- `gui-dev`:
  - `npm run typecheck`;
  - `npm test -- --run`;
  - schema generation.

## B.8. Критерии удаления или переноса

- Dependency можно убрать из core/default только если:
  - нет core import/subprocess consumer;
  - optional command имеет отдельный install profile;
  - отсутствие dependency даёт раннюю actionable ошибку, а не traceback;
  - clean core smoke проходит;
  - focused feature smoke проходит после установки extra;
  - docs, doctor, CI и package metadata синхронизированы.
- Dependency нельзя удалять, если:
  - роль только кажется transitive;
  - external executable используется runtime;
  - package нужен advertised broad-input support;
  - platform distributions расходятся;
  - отсутствие покрывается лишь skipped tests;
  - нет clean-environment proof.
- Одобренные границы переноса:
  - `nodejs` → GUI/dev;
  - `muscle`, `clustalo` → alternative MSA;
  - xTB уже удалён, не создавать для него refinement profile;
  - PySCF, Plotly, MDAnalysis, ACPYPE и ParmEd сохраняются в Standard/Desktop.
- Не считать безопасными кандидатами на удаление без дополнительной проверки:
  - `hmmer`;
  - standalone `reduce`;
  - `openbabel`;
  - `ambertools`;
  - `modeller`;
  - `rdkit` при сохранении текущей broad chemistry behavior.

## B.9. Проверки репозитория после ревизии

- Запустить focused tests каждой dependency-группы.
- Запустить обязательные общие checks:
  - `pytest -m 'not slow' -q`;
  - `ruff check src/dvbfixer`;
  - `mypy src/dvbfixer/cli.py src/dvbfixer/ffutils src/dvbfixer/pdbutils src/dvbfixer/align.py`;
  - `python scripts/check_agent_docs.py`;
  - `python scripts/gen_cli_reference.py --check`;
  - `python scripts/gen_gui_spec.py --check`;
  - `git diff --check`.
- Для GUI/environment changes:
  - `npm run typecheck`;
  - `npm test -- --run`.
- Проверить установки минимум двумя независимыми путями:
  - clean core pip/venv route, если она заявлена;
  - full conda/micromamba route.
- На macOS Docker/VirtioFS не размещать micromamba root на host bind mount; использовать native container overlay согласно repository instructions.

## B.10. Результаты ревизии

- Подготовить отдельный dependency audit artifact в документации с полями:
  - dependency;
  - version/channel;
  - category;
  - concrete consumers;
  - profile;
  - keep/move/remove decision;
  - uncertainty;
  - verification command.
- Обновить installation matrix и `doctor` output.
- Не смешивать environment cleanup с diffusion scientific claims.
- Не удалять пакеты только ради уменьшения YAML.
- Итогом первого change set должны быть:
  - согласованная dependency model;
  - optional profiles/specs;
  - исправленные contradictions;
  - clean-environment evidence;
  - только затем фактическое slimming.

## B.11. Stop/go критерии

- **STOP:** dependency считается unused только по отсутствию import.
- **STOP:** pip и conda после изменения обещают разные возможности без документации.
- **STOP:** optional command падает поздним `ImportError`/`FileNotFoundError` вместо preflight diagnostic.
- **STOP:** удаление ломает advertised PDB/PDBx, preparation, MODELLER или chemistry coverage.
- **STOP:** HMMER/Reduce/Open Babel удаляются без platform/runtime proof.
- **STOP:** diffusion ML stack попадает в core environment.
- **GO:** dependency перемещена только после profile definition, clean install и focused smoke.
- **GO:** удаление означает отсутствие доказанного consumer и успешно подтверждённую работу всех заявленных profiles.

## B.12. xTB на карандаше: проверить до решения об установочном профиле

**Исторический раздел.** Исследование завершено; по
[ADR 0012](../adr/0012-remove-xtb-refinement-and-repair-probe.md) интеграция и
управляемая зависимость xTB удалены. Перечень ниже сохраняет исходный протокол,
а не блокирует одобренную реализацию профилей и не разрешает возвращать xTB.

**Решение на 2026-10-27:** xTB пока сохраняем. Удаление, обновление пакета и
изменение алгоритма не согласованы. Следующие пункты — план проверок, не
результаты выполненных экспериментов.

### Текущая интеграция и upstream

- `src/dvbfixer/minimize/refine.py` запускает `xtb --opt --gfnff`: это
  оптимизация геометрии GFN-FF, а не полуэмпирическая QM-оптимизация GFN2-xTB.
  Проход включается явно после OpenMM (`minimize --xtb-refine`,
  `zbs --refine xtb`), не заменяет параметризацию и не создаёт MD-топологию.
- В режиме гетерогенов извлекаются все гетерогены и ковалентно связанные
  остатки-якоря. Полное белковое окружение кармана не включается. Фиксируются
  выбранные атомы якорей, а не обязательно все атомы таких остатков.
- При обратном объединении исходные координаты фиксированных атомов
  сохраняются, но дрейф в сыром результате xTB не измеряется. Это не доказывает
  корректность оптимизированной относительно якоря геометрии.
- Свыше 5000 атомов код автоматически переключается на гетерогены. При
  отсутствии executable или сбое проход может пропускаться с предупреждением.
- На дату проверки локальные conda-метаданные показывают xTB 6.7.1; GitHub
  latest release и последний тег также указывают на v6.7.1.
- Upstream [PR #1429](https://github.com/grimme-lab/xtb/pull/1429) принят
  2026-08-09, merge commit `a45f05d397268fd683b5f9f824ed59912c5490b5`.
  Исправление касается дрейфа фиксированных атомов в L-ANC с фрагментированными
  гессианами и закрывает [issue #800](https://github.com/grimme-lab/xtb/issues/800).
  Это не доказывает наличие исправления в установленной conda-сборке или
  устранение всех возможных причин дрейфа.
- В issue #800 указан обходной путь: обычный ANC через `$opt engine=rf`.
  Текущий адаптер DVBFixer его явно не выбирает.

### Обязательные проверки

- [ ] Проверить актуальные релизы, conda-рецепты и backport #1429 по ревизии
  исходников/сборки, а не только строке версии.
- [ ] Сравнить 6.7.1 с текущим оптимизатором, 6.7.1 с ANC `engine=rf` и
  зафиксированную сборку с исправлением #1429.
- [ ] Проверить Linux x86_64 и macOS arm64 отдельно; сохранить версии,
  build strings, команды, параметры и полные логи.
- [ ] Измерять максимальный и RMS дрейф фиксированных атомов в сыром
  `xtbopt.xyz` до объединения. Различать общий перенос/поворот и внутреннюю
  деформацию; итоговую неподвижность проверять отдельно. Допуски определить
  до оценки результатов.
- [ ] Проверить нековалентный лиганд, ковалентный якорь/гликан и систему,
  реально активирующую путь с фрагментированным гессианом.
- [ ] Проверить число, порядок и элементы атомов, конечность координат,
  фактическую сходимость, обработку timeout и ошибок.
- [ ] Сравнить с OpenMM-only: связи/углы, хиральность, планарность, clashes,
  перемещение лиганда и контакты с карманом. Снижение энергии разных методов
  между собой не считать доказательством улучшения.
- [ ] Проверить ковалентный интерфейс после возврата исходных координат
  фиксированных атомов; не считать восстановление координат исправлением
  дрейфа внутри оптимизации.
- [ ] Проверить эффект отсутствующего кармана и совместного извлечения
  нескольких гетерогенов. Разделить внутреннюю оптимизацию молекулы и
  оптимизацию в окружении. Для последней отдельно спроектировать контекст,
  границы/кэппинг и фиксацию; не объявлять её уже реализованной.
- [ ] Оценить GFN2-xTB как отдельный явный режим: заряд, спин, растворитель,
  стоимость, поддерживаемые входы и пользу. Не подменять GFN-FF автоматически.
- [ ] Проверить переключение >5000 атомов и отсутствие гетерогенов;
  исключить незаметную смену смысла расчёта.
- [ ] Определить контракт явного запроса: ранняя ошибка либо явно отражённый
  пропуск при недоступности/сбое; сообщать фактический метод и статус.
- [ ] Dev-прогоны `minimize`/`zbs` выполнять с `--no-solvent`.

### Критерий решения

По результатам выбрать: сохранить текущую интеграцию, обновить/исправить её,
вынести в отдельный профиль либо предложить удаление. Нужны подтверждённая
польза относительно OpenMM-only, контроль фиксированных атомов, отсутствие
ухудшения химической геометрии и прозрачные ограничения по окружению кармана.
Успешный exit code сам по себе не доказывает корректность результата.

После проверки синхронизировать этот аудит, `docs/commands/minimize.md` и
инструкции установки. Подробные пояснения хранить в документации, а не
раздувать комментарии YAML.

## B.13. Принятые решения аудита и дополнительные проверки

Решения пользователя на 2026-10-27. Это согласованный объём будущей работы;
перечисленные проверки и миграция packaging ещё не выполнены.

- Параметризацию и все QM-интеграции пока сохранить. PSI4/psiresp — кандидат
  на будущее, не текущее удаление из-за сложности установки и платформенной
  поддержки. PySCF и Gaussian этим решением не отменяются.
- Сохранить все сценарии топологий: обычный `top`, ACPYPE, ParmEd, GRO,
  пользовательские force-field directories и совместимость с CHARMM-GUI.
- Сохранить MODELLER, диффузионные бэкенды, пользовательский PIR и legacy
  cc89 replay. Не подменять исторический replay профилями новых запросов.
- MDAnalysis, Plotly, `cluster` и текущие входные форматы сохранить, включая
  Plotly в текущей комплектации. Проверять зависимости до анализа, в том
  числе Plotly при `--plot`, совместимость пары topology/trajectory и выдавать
  понятные ошибки неподдерживаемого TPR.
- Сохранить RDKit, Open Babel, PROPKA, Reduce и opt-in `tleap-reduce`.
  Сохраняется и интеграция Probe до отдельного сравнительного исследования.

### Химический стек: проверки и единый источник Open Babel

- [ ] Выбрать один проверенный источник установки Open Babel на каждый
  поддерживаемый маршрут/профиль. Не накладывать conda `openbabel` и pip
  `openbabel-wheel` в одном окружении. Проверить транзитивные требования
  ACPYPE и package metadata, а не обходить конфликт без проверки через
  `--no-deps`. Если маршруты различаются, явно описать и проверить каждый.
- [ ] Синхронизировать environment specs, pip metadata, инструкции и doctor
  с выбранной политикой; проверить чистую установку и последующее обновление.
- [ ] Выполнить реальные import/feature smoke RDKit и Open Babel: bindings,
  доступные executables, CIF-плагины, извлечение SDF и гетероген-H пути.
- [ ] Проверить authoritative SMILES отдельно от автоматического восприятия
  химии и сохранить исключения для GLYCAM-шаблонного пути.
- [ ] Проверять executable Reduce отдельно от наличия AmberTools. Исправить
  утверждения, будто любой AmberTools обязательно предоставляет Reduce.
- [ ] Различать явное отключение PROPKA/Reduce пользователем, недоступную
  зависимость и ошибку расчёта. Показывать фактически выполненные этапы и
  fallback; определить раннюю проверку и политику ошибок без скрытого
  ухудшения результата.
- [ ] Для `tleap-reduce` до расчёта проверить executables и допустимость
  входной химии; legacy остаётся default и broad-input маршрутом.

### Probe: дополнительный компонент проверки, не обязательная зависимость

**Обновление 2026-10-28:** сравнение завершено, производственный адаптер исправлен
и проверен на ограниченном эталонном корпусе. См. ADR 0012. Это не общая
precision/recall-квалификация. Исторический checklist ниже не требует заново
открывать исследование перед реализацией профилей.

Текущий `diagnose` предпочитает executable MolProbity Probe, если он доступен,
и использует Python-проверку при его отсутствии/ошибке. Основной YAML Probe
не устанавливает. Это отдельный инструмент, не MolProbity Reduce.

- [ ] Сравнить Probe и Python engine на одинаковых входах с согласованной
  подготовкой водородов: чистые структуры, реальные clashes, H-bonds,
  ковалентные контакты, гликаны и лиганды.
- [ ] Зафиксировать версии, параметры, метрики и атомные пары находок;
  различия порогов, радиусов и классификации не выдавать за прямую разницу
  качества. Отдельно оценить дополнительные полезные находки, ложные
  срабатывания, пропуски, стоимость установки и время выполнения.
- [ ] Отражать фактический engine и fallback в отчёте; проверить ошибку,
  timeout, недоступность и несовместимый вывод Probe.
- [ ] Если значимой практической пользы относительно Python-проверки не
  обнаружено, предложить удаление интеграции отдельным решением. Сейчас
  Probe не удалять и не делать обязательным компонентом установки.

### Комментарии в YAML

- [ ] Очистить основной `environment.yml` и спецификации профилей от
  избыточных комментариев. Назначение пакетов, альтернативные установки,
  историю решений и подробные ограничения перенести в соответствующие
  разделы документации; оставить короткие необходимые пояснения к
  нетривиальным ограничениям версий и ссылки на подробности.

## B.14. Установочный helper и менеджеры окружений

Согласованное направление: единый helper для Standard/Desktop и дополнений,
с интерактивным режимом и неинтерактивными параметрами. Реализация одобрена
2026-10-28 после завершения исследований xTB/Probe, с PySCF в Standard/Desktop.
Продолжать по [handoff](installation-profiles-handoff.md); базовый helper уже
реализован, полная приёмка и отдельные runtimes ещё не завершены.

- Поддержать семейство conda/mamba: micromamba как рекомендуемый менеджер,
  conda и mamba как альтернативные маршруты. Проверять доступность, версию
  и пригодность менеджера; не требовать замены имеющейся conda на micromamba.
- Предусмотреть явный выбор менеджера. Автоопределение не должно незаметно
  менять выбранную пользователем установку менеджера или целевое окружение.
- uv и Poetry обсуждены и исключены из текущего объёма по решению
  пользователя. Не добавлять для них адаптеры, обязательные зависимости
  или отдельные lock-файлы. Обычный pip остаётся средством установки
  Python-пакетов внутри управляемого окружения там, где это требует спецификация,
  а не отдельным пользовательским профилем.
- Использовать одну модель комплектации, проверки и критерии готовности для
  всех адаптеров. Создавать новые изолированные окружения; существующие
  окружения и системный Python не изменять без явного согласия.
- Сохранять документированный ручной путь установки без helper.
- Первая версия — небольшой Python-helper в репозитории, не собственный
  менеджер пакетов и не отдельный GUI-установщик. Выбор Standard/Desktop,
  дополнений и каталога должен быть доступен интерактивно и для автоматизации.
- Desktop устанавливает Node/GUI; diffusion-дополнение подготавливает
  отдельные sampling/refinement runtimes, проверяет исходники/патчи и веса,
  регистрирует launcher. Не смешивать ML-стек с основным окружением.
- После установки запускать doctor и короткие функциональные smoke-проверки.
  Показывать «готово / требуется действие / ошибка», а не считать успешное
  разрешение зависимостей доказательством готовности продукта.
- Явно объяснять ручные этапы: лицензия MODELLER, предоставление checkpoints
  по условиям распространения, внешний Gaussian и системный NVIDIA-драйвер
  при соответствующем выборе. Проверять совместимость ОС/архитектуры и
  ускорителя; не скачивать веса молча и не обещать непроверенную acceptance.
- Обычная установка не должна требовать sudo или исполнения `curl | bash`.

### Лицензия MODELLER

- Предусмотреть скрытый интерактивный ввод ключа со ссылкой на регистрацию:
  https://salilab.org/modeller/registration.html.
- Настраивать ключ только в целевом окружении способом, поддерживаемым
  выбранной версией MODELLER. Не записывать его в репозиторий, логи, отчёты,
  команды или манифесты; минимизировать доступ к файлу с ключом.
- Поддержать пропуск с явным статусом «моделирование не готово — требуется
  лицензия» и последующее завершение настройки.
- Для автоматизации предусмотреть защищённый источник секрета, не аргумент
  командной строки. Проверять лицензию реальным импортом MODELLER, не
  раскрывая ключ в диагностике.

### Проверки менеджеров

- [ ] Clean-install Standard/Desktop через micromamba, mamba и conda на заявленных
  платформах; проверять реальные команды, а не только успешный solve/import.
- [ ] Проверить запуск без activation через выбранный менеджер, PATH для
  дочерних executables, пути с пробелами, частичные установки и ошибки.
- [ ] Проверить отсутствие дрейфа состава между адаптерами, единый источник
  Open Babel и работоспособность лицензии MODELLER.
- [ ] Зафиксировать поддерживаемые версии и различия conda/mamba/micromamba;
  не объявлять адаптер проверенным до clean-install и функциональных smoke.
