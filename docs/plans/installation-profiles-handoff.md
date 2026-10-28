# Handoff: реализация профилей установки Standard/Desktop

Дата: 2026-10-28.

**Статус: реализация начата 2026-10-28; приёмка ещё не завершена.** Пользователь
одобрил требования и порядок работ. Реализация и фактическая матрица проверок
описаны в [profile acceptance](../installation-profile-acceptance.md).
Повторное согласование уже одобренного состава не требуется.

## 1. Стартовая точка и приоритет решений

- Ветка: `feature/installation-profiles`.
- Исходный HEAD: `cc0716cb4bdce2a1efe3d05944e43ab749a458b5`.
- На момент подготовки HEAD, локальный `main` и `origin/main` совпадают.
- `test/probe-xtb` уже влит в `main` и отправлен на GitHub.
- Текущая версия проекта: **0.10.1**. Бамп следующей версии не входит в handoff.
- Пакет этих handoff-документов пока находится в рабочем дереве, без коммита.
  Перед началом проверить реальный Git status и сохранить эти изменения.

Прочитать [AGENTS.md](../../AGENTS.md), этот документ, затем актуальные B.4,
B.13 и B.14 [аудита зависимостей](environment-dependency-audit.md).
Настоящее одобрение, включая поправку про PySCF, имеет приоритет над прежним
конструктором core/full/extras и предложением вынести PySCF. Описания старых
xTB-проверок в аудите — исторические, не новый блокирующий этап.

**Критическая поправка пользователя:** PySCF уже включён в `environment.yml`
через pip и должен остаться в **Standard и Desktop**. RESP/PySCF не превращать
в отдельное дополнение. MDAnalysis и Plotly также сохраняются в Standard.

### Git и существующая работа

- Корневой checkout уже находится на нужной новой ветке. `main` занят отдельным
  worktree `.claude/worktrees/integration-main`; не пытаться освобождать его
  удалением worktree или force-переключениями.
- В корне есть посторонние untracked-файлы/каталоги: `.claude/`, `opencode.json`,
  `docs/plans/diagnose-gap-geometry-metrics.md`, `gap_visual_review*`. Не удалять,
  не добавлять их в коммиты и не считать результатами этой задачи.
- `.artifacts/` содержит локальные окружения и доказательства прошлых опытов.
  Не очищать их и не коммитить рекурсивно.
- Подготовка handoff не означает запрос на автоматический commit/push. Выполнять
  Git-публикацию только согласно дальнейшему явному указанию пользователя.

## 2. Одобренный пользовательский контракт

| Комплектация | Обязательные возможности |
|---|---|
| **Standard** | Обычная научная функциональность CLI: подготовка белков и поддерживаемых комплексов, MODELLER/Homology, протонирование, OpenMM-минимизация, лигандная параметризация, топологии/ACPYPE, антительная нумерация, анализ траекторий, MAFFT и **RESP/PySCF** |
| **Desktop** | **Standard + GUI**, Node.js/npm и сборка/запуск интерфейса; научный состав наследуется, а не дублируется вторым списком |

В Standard сохраняются OpenMM/PDBFixer, Biopython, Gemmi, SciPy,
PROPKA/Reduce, MODELLER, RDKit/Open Babel, OpenFF/openmmforcefields,
AmberTools/ParmEd/ACPYPE, ANARCI/HMMER, MDAnalysis/Plotly, MAFFT и PySCF.
Точную матрицу пакетов, каналов и данных зафиксировать по реальным consumers
и metadata, а не по одному поиску import.

**Отдельные дополнения:**

- платформенные diffusion runtimes для Protenix CUDA и Protpardelle CUDA/MPS;
- PSI4/psiresp в отдельном окружении;
- интеграция с внешним Gaussian, установленным оператором;
- альтернативные MSA: MUSCLE **5** и Clustal Omega.

**Probe** остаётся опциональным диагностическим компонентом; учитывать наличие,
поддерживаемую версию и фактический fallback. **xTB не входит ни в один профиль**
и не возвращается в production CLI, metadata или environment specs. Решение и
проверки уже зафиксированы в [ADR 0012](../adr/0012-remove-xtb-refinement-and-repair-probe.md).
Опциональный OpenBabel post-pass не является новой гарантией улучшения геометрии.

### Менеджеры и способ установки

- Один небольшой Python-helper в репозитории, интерактивный и неинтерактивный.
- Micromamba рекомендован; conda и mamba — поддерживаемые альтернативы после
  соответствующих проверок. Явный выбор менеджера и целевого prefix/каталога.
- Одна модель состава, действий и готовности для всех адаптеров. Не заменять
  существующий менеджер пользователя молча.
- Создавать новые изолированные окружения. Существующие окружения и системный
  Python не изменять без отдельного явного согласия.
- Pip применяется внутри управляемого окружения по спецификации; отдельный
  пользовательский pip-only профиль не согласован.
- Сохранить ручной способ установки без helper. Не вводить uv/Poetry, собственный
  resolver, GUI-установщик, обязательный sudo или `curl | bash`.

### Готовность и ручные этапы

- После установки показывать **«готово / требуется действие / ошибка»**;
  успешный solve/import не равен готовности всего продукта.
- MODELLER: скрытый ввод ключа, настройка только в целевом окружении, проверка
  реальным импортом. Возможны пропуск и последующее завершение настройки.
  Секрет не должен попадать в CLI-аргументы, логи, отчёты, Git или манифесты.
  Для автоматизации нужен защищённый источник секрета, не его значение в argv.
- Diffusion: отдельные sampling/refinement runtimes, исходники/патчи, совместимость
  ОС/архитектуры/ускорителя, проверка предоставленных checkpoints и регистрация
  launcher. Не скачивать веса молча и не смешивать PyTorch/CUDA/MPS с Standard.
- Новые Protpardelle-запросы используют **cc94**; cc89 остаётся явным историческим
  replay. Не переносить cc89-интервенции/acceptance на cc94. MPS acceptance cc94
  пока не завершён; работоспособность установщика не изменяет этот статус.
- Gaussian и системный NVIDIA-драйвер — внешние операторские компоненты, если
  выбран соответствующий маршрут.

## 3. Что уже есть и где смотреть

Таблица ниже — исходный снимок handoff до реализации. Актуальные компоненты:
`src/dvbfixer/install.py`, `installation_profiles.py`, `profile_readiness.py`,
`data/installation-profiles.toml`, `scripts/install_dvbfixer.py`, generated
`environments/` и [acceptance](../installation-profile-acceptance.md).

| Файл/область | Текущее состояние и роль в реализации |
|---|---|
| `environment.yml` | Единый прежний полный набор: сейчас включает Node, все три MSA и PySCF; xTB уже удалён |
| `pyproject.toml` | Только extra `dev`; обязательные Python dependencies включают MDAnalysis, Open Babel wheel, ParmEd, ACPYPE и Plotly; это ещё не новая модель профилей |
| `src/dvbfixer/doctor.py` | Capabilities пакетов/executables и diffusion preflight; нет общего установщика Standard/Desktop |
| `scripts/install_diffusion_backend.py` | **Регистрация launcher для уже подготовленного backend**, а не создание его окружения, установка движка или доставка checkpoints |
| `deploy/protenix-v1/`, `deploy/protpardelle-1c/` | Существующие runners, исходники/патчи, environment inputs и инструкции; переиспользовать и проверять, не подменять новой научной реализацией |
| `deploy/openmm-cuda122/` | Существующий отдельный refinement runtime |
| `docs/installation.md` | Текущий ручной путь и ограничения; PySCF уже документирован как включённый в окружение |
| `tests/test_environment_policy.py` | Политика Python/NumPy/BLAS и native linalg smoke; учитывать при изменении состава/specs |
| `tests/test_doctor.py`, `tests/test_install_diffusion_backend.py` | Существующие контракты readiness и регистрации |
| `.github/workflows/ci.yml` | Нынешние fast/full/GUI/deployment lanes; ссылки на старый `environment.yml` и install-команды требуют согласованной миграции |

Имена вроде `scripts/install_dvbfixer.py`, каталога environment overlays и нового
модуля профилей — **предлагаемые точки реализации**, не существующие команды и
не уже принятые API. Сначала определить минимальную композиционную модель; не
копировать научный список по всем менеджерам и файлам.

## 4. Последовательность реализации

### WP1 — зафиксировать состав и контракты

- [x] Составить dependency/feature matrix: пакет, ограничение, канал, consumers,
  executable/data/лицензия, принадлежность профилю, доказательство проверки.
- [x] Выбрать один лёгкий machine-readable источник научного состава и композиции
  Standard/Desktop/дополнений; отделить его от adapter-specific команд.
- [x] Зафиксировать readiness-модель и ранние actionable ошибки: отсутствие
  инструмента, неподдерживаемая платформа, несовместимый ABI, лицензия и данные.
- [x] Проверить ACPYPE metadata и транзитивный Open Babel. Определить один
  проверенный источник bindings на маршрут: не наслаивать conda `openbabel` и
  `openbabel-wheel`. `--no-deps` не является заменой этой проверки.

WP1 evidence: [матрица зависимостей](../installation-dependency-matrix.md),
[контракт и результаты](../installation-profile-acceptance.md), canonical
`src/dvbfixer/data/installation-profiles.toml`. Выбран wheel-only Open Babel
с обычным pip ACPYPE; чистый conda-путь ACPYPE отвергнут из-за несовпадения metadata.

### WP2 — согласовать спецификации и metadata

- [x] Подготовить Standard и Desktop-наследование, а также specs дополнений.
- [x] Согласовать `pyproject.toml` с реально поддерживаемым маршрутом, не обещая
  conda-only/licensed capabilities через неподтверждённый pip-only install.
- [x] Сохранить Python `>=3.11,<3.14`, NumPy `<2.5` и согласованные
  `libblas`/`libcblas`/`liblapack` OpenBLAS wrappers.
- [x] Проверить реальные NumPy/SciPy solves. Standalone `ldd` не доказывает,
  какие библиотеки уже загружены в Python-процесс.
- [x] Сохранить совместимый ручной путь `environment.yml` до доказанной замены;
  не убирать Node/альтернативные MSA раньше готовности Desktop/дополнений.
- [ ] Убрать из старого YAML исторические объяснения при финальной миграции, оставив необходимые ограничения
  и короткие ссылки на документацию.

### WP3 — helper, адаптеры и настройка компонентов

- [x] Реализовать выбор профиля, in-environment дополнений, менеджера и prefix в обоих режимах.
- [x] Реализовать conda/mamba/micromamba adapters к одной модели действий;
  проверять manager version/совместимость, выводить понятный целевой состав.
- [x] Обработать пути с пробелами, запуск без activation, PATH дочерних
  executables, частичные установки, повторный запуск и неуспешные процессы.
  Prefix с пробелами рано отклоняется: реальный MAFFT smoke доказал upstream
  несовместимость. Spaced source/spec/input argv сохраняются.
- [x] Реализовать защищённую настройку MODELLER и последующее завершение
  отложенных шагов. Не читать/копировать ключ из существующего окружения.
- [x] Для Desktop установить/собрать GUI из lockfile, проверить standalone host.
- [x] Для diffusion переиспользовать существующие runners/registration; отличать
  подготовку runtime от регистрации и от научной/hardware acceptance.

### WP4 — doctor и функциональные preflight

- [ ] Добавить проверки выбранной комплектации и дополнений с actionable
  результатом, не объявляя отсутствие невыбранного optional компонента ошибкой.
- [ ] Проверять Reduce независимо от AmberTools; различать явное отключение,
  недоступность зависимости и реальную ошибку расчёта PROPKA/Reduce.
- [ ] Проверять ANARCI, HMMER и профильные данные вместе. Ошибка инструмента не
  равна отсутствию антительного домена; не подменять явно выбранную схему SEQRES.
- [ ] MAFFT — основной MSA; MUSCLE 5/Clustal Omega — явные альтернативы.
- [ ] Сохранить все заявленные форматы, `cluster`/`align`, MDAnalysis/Plotly;
  рано проверять topology/trajectory compatibility и выдавать понятную ошибку TPR.
- [ ] Проверить реальный RESP/PySCF путь как часть Standard, а не только import.

### WP5 — чистые установки и приёмка

- [ ] Unit/contract tests composition/manager adapters/ошибок, затем чистые
  изолированные установки Standard/Desktop каждым заявленным менеджером.
- [ ] Прогоны preparation legacy и opt-in `tleap-reduce`, protonation,
  diagnose, representative `minimize --no-solvent`, MODELLER/Homology.
- [ ] Лигандная химия: RDKit/Open Babel bindings и executables, CIF/SDF,
  автоматическое восприятие и authoritative SMILES, GAFF/OpenFF и GLYCAM.
- [ ] Topology/ACPYPE/ParmEd/GRO/user force fields/CHARMM-GUI-compatible paths,
  антительная нумерация, MSA, анализ/plotting и PySCF RESP.
- [ ] Desktop typecheck/tests/build и запуск standalone; дополнения проверять
  отдельно, включая их отсутствие и отложенную ручную настройку.
- [ ] Хранить команды, версии/build strings, profile composition и результаты
  без секретов. Отмечать неисполненные платформенные проверки как not-tested,
  а не как успешные skipped tests. Linux доступен на текущей машине; macOS/MPS
  или другие заявленные платформы требуют соответствующих хостов.

### WP6 — публикация документации и CI, затем завершение миграции

- [ ] Обновить installation matrix, README/quick start, known issues, doctor,
  CI lanes и необходимые агентские контракты; сохранить ручной маршрут.
- [ ] После доказанной замены перенести GUI-only/альтернативные MSA компоненты.
  PySCF, MDAnalysis/Plotly, ACPYPE/ParmEd и химия остаются в Standard.
- [ ] Подтвердить, что нет xTB route, ML-пакетов в Standard или зависимости
  научного CLI от Desktop/Node.
- [ ] Записать фактическую matrix приёмки и обновить статусы, не объявляя
  выполненными пункты только потому, что подготовлены YAML или mocks.

Общий прогресс отмечать в [каноническом checklist](implementation-checklist.md#installation-profiles-standarddesktop),
а детальные пункты и evidence — здесь. При начале/завершении реализации обновлять
статусы handoff и заметку в AGENTS.md; документальное одобрение не означает
выполненную установку.

### Продолжение после первого implementation increment

WP1/WP2 завершены. Базовый WP3 helper, защищённая настройка и отложенное завершение
MODELLER, Desktop build и standalone проверены. Оба профиля прошли реальные
scientific smokes и лицензированные model/Homology на всех трёх Linux менеджерах.
Пользователь предоставил отдельный защищённый файл ключа; значение нигде не
публикуется. Не переносить ключ из старого окружения и не объявлять
лицензирование текущим blocker: оно уже проверено на новых prefix.

Первый increment закоммичен как `2ffc9fc`. Следующий increment реализовал
**WP3 подготовку diffusion runtimes и проверенную регистрацию**: cc94 CUDA и
Protenix прошли чистые установки, kernels и frozen preflight на трёх Linux
менеджерах, а новые micromamba runtimes — публичный checkpoint-backed `model`.
Нативный MPS остаётся not-tested. Источники/патчи/веса/frozen identities
проверяются без модификации старых окружений и исходников.

Следующий незавершённый этап — **полный WP4**, затем оставшаяся научная/
платформенная матрица WP5 и WP6 migration. Подробности и пути evidence находятся
в acceptance; реализованный helper не означает завершение всего плана.

Не менять научные алгоритмы/default prep backend, правила идентичности/нумерации,
SMILES, хиральности, `.dat` или проектную лицензию ради этой миграции.
Существующие notices и границы распространения сторонних компонентов сохраняются.

## 5. Локальный запуск и проверки

В текущем shell `python` может отсутствовать вне окружения. На этой машине
научное окружение называется `dvbfixer`, в нём Python 3.13.15. Его установленный
editable package ранее указывал на другой worktree; **явно выбирать исходники**:

```bash
export PYTHONPATH="$PWD/src"
micromamba run -n dvbfixer python -c "import dvbfixer; print(dvbfixer.__version__, dvbfixer.__file__)"
```

Проверить, что `__file__` относится к текущему checkout, а дочерние bare
`dvbfixer`/AmberTools команды видят `bin` выбранного окружения на PATH.
Наличие старого `dvbfixer` не является clean-install proof нового Standard.

Focused tests запускать первыми, затем обязательные проверки:

```bash
pytest -m 'not slow' -q
ruff check src/dvbfixer
mypy src/dvbfixer/cli.py src/dvbfixer/ffutils src/dvbfixer/pdbutils src/dvbfixer/align.py
python scripts/check_agent_docs.py
python scripts/check_user_docs.py
python scripts/check_versions.py
python scripts/gen_cli_reference.py --check
python scripts/gen_gui_spec.py --check
git diff --check
```

**CLI reference генерировать на Python 3.11, как в CI.** Argparse 3.13 меняет
выравнивание/показ aliases и даёт посторонний diff по множеству страниц. На этой
машине для генерации/проверки уже подготовлено временное tools-окружение:

```bash
PYTHONPATH="$PWD/src" micromamba run -p /tmp/opencode/dvbfixer-docs311 python scripts/gen_cli_reference.py --check
PYTHONPATH="$PWD/src" micromamba run -p /tmp/opencode/dvbfixer-docs311 python scripts/gen_gui_spec.py --check
```

Это локальное удобство, а не новый обязательный installation profile. Если
окружение не сохранилось, подготовить новый изолированный аналог baseline CI.
Справку `docs/reference/` и `gui/server/generated-dvbfixer-spec.ts` не править
вручную: после argparse changes запускать соответствующие генераторы.

Для GUI из `gui/`:

```bash
npm run typecheck
npm test -- --run
npm run build
```

Все dev-прогоны `minimize`/`zbs` — с `--no-solvent`. На Docker/VirtioFS создавать
окружения на native overlay, не на проблемном host bind mount; см. AGENTS.md.
Для GitHub здесь доступен настроенный `gh` **внутри** `micromamba run -n dvbfixer`.
Git-токены, MODELLER keys и checkpoints не включать в отчёты/коммиты.

## 6. Известная исходная проверка — не приёмка будущих профилей

На объединённом `main` `cc0716c` перед началом этой задачи:

- **1659 passed, 7 skipped, 31 deselected** в fast Python suite; шесть известных
  SWIG/closed-stream warnings. Skips: один biotite и шесть GROMACS.
- **403 GUI tests passed**, typecheck прошёл.
- Ruff, strict-target mypy, версии, docs/agent checks и генераторы прошли.
- Логи локально: `.artifacts/probe-xtb/merged-main-fast.log` и
  `.artifacts/probe-xtb/merged-main-fast.xml`; они игнорируются Git.

Этот раздел описывает исторический baseline до реализации. Текущие результаты:
WP1/WP2 реализованы; базовый helper, лицензирование и Desktop проверены на трёх
Linux менеджерах. WP3 отдельные runtimes, полноценный WP4 и остаток WP5/WP6
ещё не завершены; см. актуальный [acceptance](../installation-profile-acceptance.md).

## 7. Готовый стартовый запрос после обновления сессии

```text
Продолжай в feature/installation-profiles реализацию одобренного плана из
docs/plans/installation-profiles-handoff.md. Прочитай AGENTS.md и актуальные
решения в docs/plans/environment-dependency-audit.md. План уже одобрен:
Standard/Desktop, PySCF/RESP и MDAnalysis/Plotly остаются в обоих, Desktop
наследует Standard и добавляет GUI, xTB удалён, Probe опционален.
Начни с WP1 и затем проходи этапы по порядку с чистыми установками и smoke.
Не ограничивайся YAML/mocks; сохраняй научные возможности и существующие
окружения. Не трогай посторонние untracked-файлы. Сначала проверь Git status
и путь импортируемого dvbfixer. Подготовленные handoff-документы могут быть
ещё не закоммичены. Не запрашивай повторно одобрение уже согласованного плана.
```
