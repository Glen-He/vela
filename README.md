# Vela

Vela 用于研究裸 P15 环肽与人源 CK2α、CK2α′催化亚基的潜在直接结合。项目先在
不提供结合位置约束的条件下发现候选 site，再通过结合态复现、局部全原子精修、
多状态序列设计和多重复 MD 逐层验证。计算结果用于形成可实验检验的结构假设，
不能单独证明酶抑制、细胞效应或抗肿瘤疗效。

阶段一已经完成，并可把 10 个受体的原始/审计资产、P15 化学记录、4 个
apo/apo-like 主受体和 1 个资格控制受体直接交给阶段二。阶段二当前使用
CABS-dock 3.0.12 的官方
`--ca-rest-add` 表达二硫环的粗粒化拓扑，不再依赖第三方源码补丁。旧 SG 约束协议的
三个运行批次均只保留为历史诊断：它们暴露了科学 QC 异常传播、粗粒化拓扑命名、
Top-1000/Top-10 候选丢失和小样本 seed 投票不稳定等问题。当前候选从完整 TRAF 做
site-first、pose-second 去重。新 CA 协议的 `120655..120662` 开发批次已完整结束：
完整轨迹能够重复采样实验姿态，固定候选预算能够保留实验位点，但不保证保留同一低频
精确姿态；因此阶段二资格已明确拆成“完整池精确姿态采样”和“有限候选位点保留”，
精确姿态恢复由阶段三局部方法控制。该开发批次同时作为全原子拓扑校准来源；校准冻结后
必须使用 `120663..120670` 执行正式留出资格。阶段三已经实现配置驱动的局部恢复控制、
去配体结合态全表面复现、粗粒化 pose 全原子交接、正式候选 FlexPepDock 运行以及
QC/聚类报告，并补齐跨状态 site 比较、实验骨架 guided 起点、全酶环境映射和候选
事实卡；真实控制、阶段二正式 candidate 和校准阈值尚未产生，因此各放行门槛仍会
阻止正式采样或精修。阶段四的软件链路也已建立，但只接受阶段三 blind-supported
多模板证据；当前没有真实上游 candidate，因此尚未启动正式序列优化。

## 环境与入口

项目要求 Python 3.12，并统一使用 `uv`：

```bash
uv sync
uv run vela config check
```

项目依赖三个项目外的外部运行时，其可执行文件路径在当前机器的 `configs/` 下配置为本机
绝对路径，更换环境（机器）后需同步调整：

- CABS-dock：`configs/discovery.toml` 的 `[discovery.cabsdock]` 的 `executable` 与 `source_dir`；
- cg2all 与 Rosetta：`configs/validation.toml` 的 `executable`、`scripts_executable`、`package_metadata` 与 `checkpoint`；
- Pyright 的 `scripts` 类型解析：项目根 `.env` 的 `PYTHONPATH` 指向本机 CABS 源码根（模板见 `.env.example`，更换环境后更新）。

`uv run vela validation tool-check` 会核对 cg2all 与 Rosetta 的来源、版本和哈希，可用来确认路径配置正确。

阶段一建议使用阶段级入口：

```bash
uv run vela preparation status
uv run vela preparation run
```

`status` 只检查阶段一化学、下载、审计和制备产物，不检查阶段二方法参数；`run` 一次
执行化学记录、10 个受体下载和审计、4 个 apo/apo-like 主受体及 `3Q9X_A` 资格控制
受体制备，再按 manifest
路径、参数和 SHA-256 检查交付完整性。只有阶段一产物不完整或过期时才返回退出码 3。
需要单独排查时仍可运行底层命令：

```bash
uv run vela chemistry record
uv run vela receptors download
uv run vela receptors audit
uv run vela receptors prepare
```

`chemistry record` 会区分“配置结构有效”和“科学字段已经关闭”。当前裸 P15 主状态
记录为 `production_ready=true`；这只关闭化学身份，不代表全表面对接方法已经通过
资格验证。

受体下载通过 RCSB 的 mmCIF 下载接口和 Data API 获取原始结构及条目元数据。下载采用
有上限的指数退避重试、临时文件校验和原子安装；再次执行时仍会重新校验本地缓存，
不会只按文件是否存在跳过。接口依据见 [RCSB 文件下载服务](https://www.rcsb.org/docs/programmatic-access/file-download-services)
和 [RCSB Data API](https://data.rcsb.org/index.html)。

阶段二当前入口为：

```bash
uv run vela discovery qualification-plan \
  --target ck2_alpha \
  --run-id <qualification-run-id>
uv run vela discovery qualification-run \
  --run-dir outputs/discovery/qualifications/<qualification-run-id>
uv run vela discovery qualification-analyze \
  --run-dir outputs/discovery/qualifications/<qualification-run-id>

uv run vela discovery topology-calibration-plan \
  --source-run outputs/discovery/qualifications/<development-run-id> \
  --run-id <topology-calibration-run-id>
uv run vela discovery topology-calibration-run \
  --run-dir outputs/discovery/topology_calibrations/<topology-calibration-run-id>
uv run vela discovery topology-calibration-analyze \
  --run-dir outputs/discovery/topology_calibrations/<topology-calibration-run-id>

uv run vela discovery status --target ck2_alpha
uv run vela discovery plan --target ck2_alpha --run-id <run-id>
uv run vela discovery run --run-dir outputs/runs/<run-id>
uv run vela discovery analyze --run-dir outputs/runs/<run-id>
```

资格流程同样一次只处理一个 target。每个资格 seed 在同一 worker 内先把 Pc 全局对接
到公开环肽基准与 4IB5 配对的 unbound `3Q9X_A`，以 4IB5 实验复合物做事后 native
评价，再运行该 target 的 apo 技术
先导；控制的科学 QC 不合格只会产生
`unqualified` 证据，不会中断后续 pilot，只有技术完整性错误才终止任务。实验配体坐标
只用于控制输出的事后评分，不进入 CABS-dock 搜索命令。完整轨迹的采样资格要求精确
姿态进入 5.5 Å FlexPepDock 吸引域；固定候选集只要求在实验位点质心 4 Å 内并恢复至少
20% 的原生 Cα 接触。候选中的精确姿态召回仍会报告，但不作为阶段二位点发现门槛。
site-first 参数在正式资格 seed 运行前已经冻结。独立 `topology-calibration-*` 流程用开发 seed 的分层样本校准粗粒化
Cα 判据能否恢复为合理全原子二硫环，并给出可采用的最大阈值。全原子化时只由 cg2all
重建肽链，受体使用按 CABS 受体对齐的实验 receptor-only 结构，避免神经网络重建受体
产生异常 Rosetta 主链应变。它不读取实验配体坐标，
也不把 Rosetta 引入正式阶段二全表面采样。资格 seed 与正式 seed 完全分离。
旧 SG 约束协议得到的 7 Å 校准不授权当前 CA 约束协议。当前协议已用
`120655..120662` 的完整 TRAF 完成 native-free 分层校准：局部精修通过肽 C-alpha
平底坐标约束保持原 site，6/7/8/10 Å 四个候选层均为 8/8 通过，10 Å 已冻结为最大
连续合格包络；10–12 Å 层只作非门控外部对照。
第一个 target 通过后，后续 target 可在 `qualification-plan` 增加
`--control-run outputs/discovery/qualifications/<passed-run-id>`，复用已经通过且逐文件
校验 SHA-256 的同方法 Pc 控制，只运行新 target 自己的 apo 先导；重复的确定性
control/seed 任务不会再次计算。

`status` 和 `plan` 每次必须显式选择一个 target。一次完整 run 只包含该 target 的
多个受体构象；CK2α 和 CK2α′必须使用不同 run ID 独立运行。`plan` 只有在 P15 化学、
当前 target 的主受体、方法资格报告、独立 seeds 和聚类规则全部冻结后才会创建运行
目录。`run` 按冻结计划执行或复用哈希验证通过的已完成任务，使用
CABS-dock 官方 `--ca-rest-add` Cys C-alpha 距离约束、关闭全原子重建，并写出 `sampling_manifest.json` 和规范
`pose_evidence.tsv`；`analyze` 随后才生成 site 报告。正式 pose 由完整 10,000 帧 TRAF
先经过 `topology_feasible`/Cα 受体接触过滤，再按 site 和 site 内姿态分层去重。每个
seed 最多保留 64 个按人口优先的 site、每个 site 最多 4 个 pose 簇，每簇保留 medoid
与最低相互作用能代表，因此单任务上限为 512；CABS 原始 Top-1000/Top-10 单独保留为
非门控基线。所有阶段二 pose
都是粗粒化证据，不冒充全原子结构；
Rosetta/FlexPepDock 只在阶段三出现。

所有阶段的不可变计划同时记录发行版本和 `src/vela` 源码/包内资源 SHA-256。执行、恢复
和分析都会核对该摘要；修改源码后必须重新建立计划，不能用未变化的 `0.1.0` 版本号
继续消费旧计划。

当前每个 target 使用 8 个正式 seed（`120623..120630`）和 2 个 apo/apo-like
受体构象，因此一次 run 有 16 个单进程 CABS 任务。调度单元是 seed：
最多 8 个 seed worker 并发，每个 worker 串行完成所领取 seed 下当前 target 的所有受体
构象；单进程 CABS-dock 内部不会使用 8 个 worker。
资格任务使用另一组 seed，不能与正式支持度混算；已经用于协议诊断的 calibration seed
也不能在修改选择算法后重新宣布为正式资格证据。

阶段三已经开始重构，当前可先完成不依赖阶段二候选的结构和工具准备：

```bash
uv run vela validation status
uv run vela validation prepare
uv run vela validation tool-check
uv run vela validation control-plan --run-id <control-run-id>
uv run vela validation control-run --run-dir outputs/validation/controls/<control-run-id>
```

`prepare` 从原始实验复合物生成彼此分离的去配体受体、成对实验参考、4IB5/Pc
局部方法控制输入，以及 `4MD7`/`1JWH` 的全酶 assembly 参考；`tool-check` 核对 cg2all 与
Rosetta 的来源、版本、二进制或 checkpoint 哈希和所需选项。局部恢复控制的具体结构、
扰动、seed、排名项和通过标准全部来自 `[[validation.local_controls]]`；当前项目选用
4IB5/Pc，但代码不绑定该对象。

正式 FlexPepDock 对每个起始姿态运行 4 个 seed，每个 seed 传入
`-nstruct 128`，总计 512 个 decoy。4IB5/Pc 方法资格使用两个不重叠的
`4 × 128` 批次；每批先合并全部 seed 再统一排分和聚类，seed 只用于检查
近天然构象簇是否获得跨随机流支持。资格通过要求两批都完成采样回收、
排序回收且主要恢复姿态一致，不再使用逐 seed Top-10 投票。
`validation status` 分别报告独立工具/资产准备、结合态全表面复现和正式候选局部精修
三个门槛，避免把“文件已经准备好”误读为两种方法均已放行。

全表面方法资格通过后，所有登记为 `bound_state_blind_replication` 的去配体受体由
同一组通用命令执行：

```bash
uv run vela validation replication-plan \
  --target ck2_alpha \
  --run-id <replication-run-id>
uv run vela validation replication-run \
  --run-dir outputs/validation/replications/<replication-run-id>
uv run vela validation replication-analyze \
  --run-dir outputs/validation/replications/<replication-run-id>
uv run vela validation replication-compare \
  --discovery-run outputs/runs/<discovery-run-id> \
  --run-dir outputs/validation/replications/<replication-run-id>
```

计划从受体角色、`validation.bound_states` 派生资产和冻结的阶段二全表面参数生成，
Python 不列举具体 PDB、配体、受体数量或 seed。执行器只读取冻结清单中的显式任务，
不会重新替换成 apo 主集合；结果固定标记为 `bound_state_blind_replication`，也不能被
阶段二主候选交接入口接受。

不依赖阶段二 candidate 的 Pc 骨架 guided 起点可以先独立生成：

```bash
uv run vela validation guided-plan --run-id <guided-run-id>
uv run vela validation guided-run \
  --run-dir outputs/validation/guided/<guided-run-id>
```

它按 `[[validation.guided_templates]]` 的逐位映射将 P15 线程化到标准实验环肽骨架，
并显式恢复 His tautomer、端基和二硫键。FMP37 属于非标准环肽模拟物，只保留为
site reference，不会被自动线程化。guided 起点不能增加 blind 发现支持度。

阶段二完成并冻结 candidate 后，blind 与 guided 起点共用同一局部精修内核：

```bash
uv run vela validation handoff-plan \
  --discovery-run outputs/runs/<discovery-run-id> \
  --candidate-id <reviewed-candidate-id> \
  --run-id <handoff-run-id>
uv run vela validation handoff-run \
  --run-dir outputs/validation/handoffs/<handoff-run-id>

uv run vela validation refinement-plan \
  --source-run outputs/validation/handoffs/<handoff-run-id> \
  --run-id <refinement-run-id>
uv run vela validation refinement-run \
  --run-dir outputs/validation/refinements/<refinement-run-id>
uv run vela validation refinement-analyze \
  --run-dir outputs/validation/refinements/<refinement-run-id>
uv run vela validation environment-map \
  --run-dir outputs/validation/refinements/<refinement-run-id>
uv run vela validation candidate-review \
  --discovery-run outputs/runs/<discovery-run-id> \
  --replication-run outputs/validation/replications/<replication-run-id> \
  --run-dir outputs/validation/refinements/<refinement-run-id>
```

guided 精修只需把 `--source-run` 改为 `outputs/validation/guided/<guided-run-id>`；
计划、运行和分析清单会自动保持 `known_site_information_used=true`，不会冒充 blind
candidate。`environment-map` 只报告叠合 RMSD、CK2β/另一催化链的接触数和最小距离，
不使用未校准的竞争或碰撞阈值；`candidate-review` 也只汇总已完成证据，不给亲和力、
抑制作用、抗肿瘤效果或最终等级。

`handoff` 使用配置指定的 cg2all checkpoint 恢复全原子结构，再由 RosettaScripts 恢复
端基和二硫键；候选来自多模型 `trajectory_candidates.pdb` 时，计划与执行会共同核对
`model_index` 并只重建该帧。`refinement` 对每个 blind 起点乘以独立 seed 做 prepack
和局部精修；
`refinement-analyze` 才执行配置门槛下的碰撞、界面接触、受体位移、site 保持和构象
聚类。当前这些命令的数据合同和合成案例已经验证，但真实运行仍等待 Pc 局部恢复资格、
冻结 seed/QC/聚类规则以及阶段二 candidate，不能把代码完成等同于科学阶段放行。

阶段四按“低成本完整覆盖 → 有限组合 → 昂贵柔性复核 → 可选父序列迭代”分层执行：

```bash
uv run vela design status
uv run vela design tool-check
uv run vela design single-plan \
  --source-run outputs/validation/refinements/<refinement-run-id> \
  --target-cluster <cluster-id> \
  --run-id <single-run-id>
uv run vela design screen-run --run-dir outputs/design/screens/<single-run-id>
uv run vela design screen-analyze --run-dir outputs/design/screens/<single-run-id>

uv run vela design combination-plan \
  --source-run outputs/design/screens/<single-run-id> \
  --run-id <combination-run-id>
uv run vela design screen-run --run-dir outputs/design/screens/<combination-run-id>
uv run vela design screen-analyze \
  --run-dir outputs/design/screens/<combination-run-id>

uv run vela design finalist-plan \
  --source-run outputs/design/screens/<combination-run-id> \
  --run-id <finalist-run-id>
uv run vela design finalist-run --run-dir outputs/design/finalists/<finalist-run-id>
uv run vela design finalist-analyze \
  --run-dir outputs/design/finalists/<finalist-run-id>
```

`single-plan` 完整枚举配置允许的单点替换，不固定 legacy 中没有依据的 P5/G9，也不
禁止新 His；每条 WT/候选在同一模板和 seed 下成对计算。只有通过单点门槛的替换会
提出有限二/三点组合，组合必须按完整序列重新计算。`finalist-*` 再用已通过阶段三
资格验证的 FlexPepDock 做多 seed 柔性复核，并分别检查两类 score 和结构 QC，最后
写出含 WT、候选、亚型/模板身份和结构哈希的 `md_queue.tsv`。

阶段四采用 `single_supported_target`：一次运行只优化一个获得阶段三支持的亚型，
CK2α 和 CK2α′可分别运行，任一轨道的合格候选都能进入后续队列；不要求同一序列同时
覆盖两个亚型，也不把另一个亚型当作必须削弱的反向目标。

通过柔性复核并进入 `md_queue.tsv` 的多个同代候选可以作为下一代显式父序列。下一代
不是直接枚举四/五点组合，而是对每个父序列生成增加、撤销或替换一个位点的一步邻域：

```bash
uv run vela design iteration-plan \
  --source-run outputs/design/finalists/<finalist-run-id> \
  --parent-candidate-id <candidate-id> \
  --run-id <iteration-screen-run-id>
uv run vela design screen-run \
  --run-dir outputs/design/screens/<iteration-screen-run-id>
uv run vela design screen-analyze \
  --run-dir outputs/design/screens/<iteration-screen-run-id>
uv run vela design finalist-plan \
  --source-run outputs/design/screens/<iteration-screen-run-id> \
  --run-id <next-finalist-run-id>
```

`--parent-candidate-id` 可以重复，但父序列必须来自同一代，且已经通过柔性门槛并被选入
阶段五队列。`md_selected` 只表示资源选择，不表示 MD 或实验已经成功；production 中应在
查看相应动态和实验记录后再决定是否启动下一代。每个后代都会重新走 WT 配对初筛和
FlexPepDock，父代 score 不会累加或直接继承。

序列净电荷、疏水/芳香残基计数、连续疏水段和简单化学风险 motif 会作为透明事实
输出，不在没有实验校准时自动宣布“可开发”或淘汰候选。游离环肽预组织、溶液聚集和
构象代价改由阶段五的无受体多重复 MD 及实验测定处理，不沿用 legacy 中未经当前体系
验证的 `simple_cycpep_predict` 放行规则。

Rosetta 的终端许可提示说明商业使用可能需要单独许可；如果这里的“生产”包含商业
研发或部署，运行正式批次前必须由使用方确认 Rosetta 许可范围。

## 配置、数据和文档边界

- `src/vela/resources/` 下四份阶段 TOML：随包安装的工作流默认参数；
- `configs/common.toml`、`discovery.toml`、`validation.toml`、`design.toml`：P15/CK2 项目输入及实际阶段覆盖；
- `data/`：由命令生成的原始下载、审计表和可追溯制备产物，不保存配置；
- `outputs/`：后续不可变正式运行目录；
- `docs/research.md`：实验结构事实、适用范围和限制；
- `docs/plan.md`：项目方向、阶段状态、证据规则和决策记录。
- `docs/config.md`：两个配置来源的职责、参数索引和允许覆盖的边界。

`data/` 与 `outputs/` 是可再生成的本地运行产物，已从版本控制中排除。原始 mmCIF
只读保存在 `data/receptors/raw/`；受体制备不会覆盖它，而是在
`data/receptors/prepared/` 生成带来源哈希的派生文件。

当前数据布局如下：

```text
data/
├── chemistry/<ligand-id>/chemistry_record.json
├── receptors/
│   ├── raw/           # 原始 mmCIF、RCSB 元数据、下载 manifest
│   ├── audit/         # 结构、链、组分、缺失记录和序列差异审计
│   └── prepared/      # 四个主受体的基础 mmCIF、替代构象决定和制备 manifest
└── validation/
    ├── bound_states/  # 结合态 receptor-only、pair reference 和标准方法控制派生物
    └── environments/  # 4MD7/1JWH 全酶 assembly 参考
```

配置文件中的相对路径以该配置文件所在目录解析。运行策略可通过环境变量覆盖；当前支持
`VELA_DATA_DIR`、`VELA_OUTPUTS_DIR`、`VELA_DOWNLOAD_RETRIES` 和
`VELA_DOWNLOAD_TIMEOUT_SECONDS`。

包内默认值位于 `src/vela/resources/` 下四份阶段 TOML，项目特有输入和覆盖位于
`configs/` 下固定的四份 TOML；`src/vela/config/` 只是加载与校验代码。完整参数分类见
[`docs/config.md`](docs/config.md)。

## 代码结构

`src/` 是 Python 标准 source layout 的包根目录，不代表所有代码属于同一模块。Vela
在 `src/vela/` 内按稳定职责组织：

```text
src/vela/
├── cli.py                 # 仅负责进程入口和错误边界
├── commands/              # CLI 语法与工作流适配
├── core/                  # 无阶段依赖的哈希、原子写入和类型边界
├── config/                # 配置合并、分 section 解析和领域对象组装
├── preparation/           # 阶段一
│   ├── chemistry.py
│   └── receptors/
│       ├── acquisition.py # 单文件下载、重试、校验、原子安装
│       ├── download.py    # 登记表下载与 manifest 编排
│       ├── audit/         # mmCIF、实体、晶体接触、报告、工作流
│       └── cleaning/      # 替代构象、结构清理、报告、工作流
├── discovery/             # 阶段二
│   ├── sampling/          # CABS 适配、计划、执行和证据提取
│   └── analysis/          # pose 边界、site 聚类和报告
├── validation/            # 阶段三
│   ├── bound_states/      # 实验结合态、控制和去配体复现
│   ├── refinement/        # 全原子交接、局部精修和分析
│   └── assessment/        # 全酶环境与候选综合复核
├── design/                # 阶段四
│   ├── sequence/          # 候选序列、邻域和迭代
│   ├── screening/         # 成对界面初筛
│   └── finalists/         # 柔性复核和 MD 队列
└── resources/             # 包内唯一默认配置
```

后续只在开始真实实现阶段五时增加 `dynamics/`。
不创建空的阶段目录，也不使用 `stage1/` 之类编号名称，以免把可复用的领域能力与
执行顺序耦合。

## 当前受体范围

阶段二正式盲发现主集合为 CK2α 的 `3Q04_A`、`3QA0_A`，以及 CK2α′的
`5YF9_X`、`5Y9M_A`。`3Q04_A` 和 `5YF9_X` 仅额外承担技术先导角色，不能单独形成
正式 site 结论。

`4IB5_A`、`9FBM_A` 和 `9FBI_A` 的原始复合物继续只读保留；阶段三已经从原始文件
独立派生 receptor-only 和成对参考。4IB5/Pc 另有标准局部方法控制输入，FMP37 只作
非标准配体位点参考。项目另外登记 `4MD7_E` 与 `1JWH_A` 作为全酶环境链，并已生成
对应 assembly 参考；它们不参加全表面搜索。4IB5/Pc 的 P15 guided 全原子起点已经
通过真实 Rosetta 化学恢复 smoke，但正式 guided 精修仍等待局部方法资格和正式 seed。

因此，登记结构总数现在是 9，但当前最低 P15 全表面搜索计划仍是 **7 个受体构象**：阶段二 4 个
apo/apo-like 主发现受体，加阶段三 3 个去配体结合态复现受体。现在只显示 4 个
`prepared`，是因为该数字只统计阶段一方法无关的 apo/apo-like 基础受体；另外 3 个
现已作为阶段三 receptor-only 派生资产存在。它们仍须在阶段二全局方法通过资格后才可
执行 P15 全表面复现，且结果不能与阶段二主证据混算。

## 开发检查

```bash
uv run ruff format --check src tests
uv run ruff check src tests
uv run --env-file .env pyright
uv run pytest
```

完整研究方案和阶段一剩余条件见 [`docs/plan.md`](docs/plan.md)，结构清单及已核查事实
见 [`docs/research.md`](docs/research.md)。
