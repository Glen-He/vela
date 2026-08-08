# Vela config 说明

Vela 只有两个参数来源，职责不能混用：

| 优先级 | 位置 | 职责 |
|---:|---|---|
| 低 | `src/vela/resources/` 下四份阶段 TOML | 随包安装的通用工作流默认参数，是默认值的唯一来源 |
| 高 | `configs/common.toml`、`discovery.toml`、`validation.toml`、`design.toml` | 当前 P15/CK2 研究对象和实际阶段覆盖 |

解析代码位于 `src/vela/config/`，它不是第三个参数来源。环境变量只覆盖部署相关路径和下载策略；命令行参数只表达本次运行身份或显式 config 目录。相同事实不得同时在两个活动配置中完整复制。

运行 `uv run vela config check` 可以查看合并后的下载、晶体接触、altloc 和主受体集合关键参数，同时验证全部 section 和项目条目。

## 包内默认参数

以下参数按阶段保存在 `src/vela/resources/preparation.toml`、`discovery.toml`、`validation.toml` 和 `design.toml`；项目文件未覆盖时自动生效：

| TOML section | 参数 | 作用 |
|---|---|---|
| `download` | `coordinate_base_url`、`metadata_base_url` | RCSB 坐标和元数据端点 |
| `download` | `retries`、`timeout_seconds` | 下载重试次数和单次超时 |
| `download` | `backoff_initial_seconds`、`backoff_multiplier` | 重试退避序列 |
| `download` | `chunk_size_bytes`、`user_agent` | 网络分块和请求身份 |
| `audit.crystal_contacts` | `distance_A` | 晶体接触搜索距离 |
| `audit.crystal_contacts` | `min_occupancy`、`include_hydrogens` | 接触原子占有率和氢原子策略 |
| `preparation.altloc` | `preferred_label` | 平均占有率相同时优先选择的替代构象标签 |
| `discovery.ensemble` | `min_receptors_per_target` | 每个靶标所需的最少独立受体构象 |
| `discovery.ensemble` | `allowed_structure_states` | 正式主发现允许的受体状态 |
| `discovery.cabsdock` | `seed_workers` | 同时处理的 seed 数；每个 worker 串行执行当前 target 的全部受体构象 |
| `discovery.cabsdock` | Monte Carlo、replica、温度和受体约束参数 | CABS-dock 全表面粗粒化采样预算 |
| `discovery.cabsdock` | `filtering_count`、`clustering_medoids`、`clustering_iterations` | 低能帧过滤和仅用于对照的上游 Top-10 基线 |
| `discovery.cabsdock` | 接触、粗粒化二硫拓扑和 site-first 选择参数 | 完整 TRAF 的模型级拓扑/接触过滤、位点优先和位点内姿态去重 |
| `validation` | 二硫键 SG 范围、实验界面接触距离 | 阶段三结合态链配对和结构输入质控 |
| `validation.rosetta` | 并发任务数、每 seed decoy 数、评分函数和低分辨率预优化开关 | FlexPepDock 局部协议；正式使用受资格门槛约束 |
| `validation.cg2all` | 表示、进程、batch、断链距离和 CA 保真门槛 | 阶段二 CA/SC pose 的阶段三全原子重建 |
| `validation.handoff` | 每个受体 site 的独立起点数、化学恢复 seed | 阶段二 candidate 到阶段三输入的选择预算 |
| `validation.refinement` | prepack seed、局部扰动、排名 score 和三个增量 seed 批次的大小 | candidate 的 FlexPepDock 运行参数 |
| `validation.funnel` | Stage 3A/B/C 候选预算、source 起点预算、任务单元门槛和跨受体原子兼容阈值 | P15 探索性粗到细精修漏斗 |
| `validation.analysis` | 局部精修 QC、site 保持和构象聚类门槛 | 已由 4IB5/Pc matched control 的两个独立 `4×128` 批次冻结；它们是几何完整性与重复性 QC，不是 native pose 判别器 |
| `design.screen` | 固定骨架局部邻域、评分函数、排名项和任务并发 | 完整单点与有限组合共用的低成本 WT 配对初筛 |
| `design.combination` | 组合突变数、每位置替换数和候选总预算 | 只限制提案规模，组合仍必须重新计算 |
| `design.iteration` | 同代父序列数、相对 WT 的累计突变上限和每代候选预算 | 每代只生成增加、撤销或替换一个位点的一步邻域 |
| `design.finalists` | 柔性复核候选/MD 预算、独立 seed、两类 score 和结构通过门槛 | 初筛后运行，并复用已资格化的阶段三 FlexPepDock 几何合同 |

这些参数会写入相应 audit、preparation 或 discovery manifest。修改审计或制备参数后，旧 manifest 会先被阶段一完整性检查判为过期，阶段二也会因阶段一前置条件失效而拒绝运行；必须重新运行阶段一，而不是静默复用旧产物。

## 项目 config

`configs/` 固定包含四份项目参数文件，加载器按以下顺序读取；每个具体字段只能在其中一份文件声明，重复声明会直接报错，避免静默覆盖。

| 文件 | TOML section | 内容 |
|---|---|---|
| `common.toml` | `paths`、`chemistry.ligand`、`receptors` | 跨阶段共享的路径、当前配体化学事实和唯一受体登记表 |
| `discovery.toml` | `discovery` | 阶段二 CABS-dock、受体 ensemble、资格、seed 和 site 分析规则 |
| `validation.toml` | `validation` | 阶段三 FlexPepDock、cg2all、结合态、局部控制、guided 模板和全酶环境参考 |
| `design.toml` | `design` | 阶段四序列空间、初筛、迭代和柔性复核的项目特有覆盖 |

当前已经选择 CABS-dock 作为阶段二候选方法，真实采样参数均已进入上述配置层，代码中不再散落 Monte Carlo、replica、过滤或聚类数值。正式主发现固定 8 个 seed：
`120623..120630`；`120631..120654` 已用于旧 SG 约束协议的开发、验证和根因诊断；
`120655..120662` 已用于新 CA 约束协议的位点/姿态边界开发，并作为拓扑校准来源；`120663..120678` 已用于两轮旧合同留出或开发分析；跨受体候选合同的新资格固定使用未查看的 `120679..120686`。一次 run 只选一个 target；当前每个 target 有2个受体构象，因此展开为16个 CABS 任务。`seed_workers=8` 表示最多8个 worker 并发领取 seed，每个 worker 内按配置顺序串行执行该 target 的2个受体构象；它不是 CABS-dock 单任务的内部并行数。
`discovery.targets.<target-id>` 分别保存坐标参考、技术先导受体、资格报告及哈希和 site 分析规则，不能用一个 target 的先导替另一个 target 放行。当前采用 `min_seed_support=2` 建立单受体site；双受体对应候选归为 `ensemble_consensus`，单受体至少4/8 seed支持归为 `conformation_specific`，其余为 `insufficient_evidence`。前两类分别限制为Top-32和Top-8，只有预算内候选允许进入handoff。

`discovery.qualification` 当前使用未见 seed `120679..120686`，并登记正式 apo 控制受体集合 `control_receptor_ids=["3Q04_A","3QA0_A"]`、历史公开 benchmark `benchmark_receptor_id="3Q9X_A"`和4IB5/Pc native参考。一次资格只运行16个 Pc 位点控制，不夹带 P15 pilot，也不接受共享旧控制入口。每个 seed worker 顺序完成两个受体，候选冻结后才使用 native 评价。
实验控制分别报告 topology-feasible 且接触受体的 eligible sampling pool、单受体site诊断、跨受体候选交付和上游Top-10基线。正式硬门要求两个受体的实验区域各有至少2个seed支持，并由一个预算内`ensemble_consensus`候选同时交付；单受体Top-32只作诊断。L-RMSD 5.5 Å继续报告为capture-range指标，不参与site资格。
跨受体匹配和排序不是含糊的“综合分”。匹配固定使用同一靶标/公共坐标系内的确定性 complete-linkage，距离为接触 Jaccard 与质心距离的最大归一化值，门槛为1，每个候选最多包含一个同受体site；冲突按最大簇间距离和稳定site ID处理。候选在证据等级内按无权重字典序排列：最弱受体seed支持降序、总seed支持降序、跨构象最大归一化距离升序、最弱受体`selected_pose_fraction`降序、总`selected_pose_fraction`降序、受体内site median CABS score的midrank分位升序和稳定candidate ID。原始CABS score不跨受体直接合并；`selected_pose_fraction`只表示项目侧压缩后候选占比，不是完整轨迹人口或结合概率。`candidate_id`只标识聚类结果，不表达排名；正式次序只读取`rank_within_tier`。运行计划和analysis manifest保存同一完整合同，分析参数或实现不一致时直接拒绝；旧schema不再读取。
该控制验证位点级全局召回，不是亲和力验证。

`discovery exploration-plan` 是在正式资格尚未完成、但需要评估 P15 技术链路和资源规模时使用的独立入口。它仍要求阶段一产物、方法身份、8个 seed、双受体集合和完整分析合同全部有效，只不要求把未完成的 Pc 留出伪装成 `qualified`。调用时必须显式提供一个位于 `outputs/discovery/qualifications/` 下的历史开发依据；计划会冻结其 plan、sampling、report 和 pose 证据哈希，并写入 `scope="development_only"`、`production_qualified=false` 和证据限制。探索结果不得直接进入正式候选结论或覆盖 `discovery.targets.<target-id>.qualification_status`；正式 `discovery plan` 与结合态全表面复现仍只接受当前配置登记的有效资格报告。

`discovery.qualification.topology_calibration` 把开发 seed 的完整 TRAF 候选按 Cys 端点 Cα 距离分层，按 seed、site、CABS 相互作用能和 Cα 距离分位覆盖抽样；选中后才调用 CABS 自身运行时物化 CA/SC 帧，SC 距离只作描述性特征。候选包络为 6/7/8/10 Å，12 Å 层只作外部对照。每个样本先由 cg2all 重建肽链，再把实验 receptor-only 结构对齐到 CABS 受体并替换神经网络重建的受体，然后经过 RosettaScripts `ForceDisulfides`、FlexPepDock prepack 和一个无 native、带肽 C-alpha 平底坐标约束的局部 refine；约束平底宽度 2 Å、标准差 1 Å、权重 1，只防止闭环优化变成第二次无约束位姿搜索，不替代事后的 site 与接触保持判据。

`validation.bound_states` 对标准肽控制显式声明 `ligand_n_terminus` 和 `ligand_c_terminus`；当前 Pc 使用实验合成定义 `Acetyl`/`CONH2`，P15 主化学仍由 preparation chemistry 独立定义。拓扑校准与 candidate 交接调用同一恢复流水线，交接在恢复端基后会再次执行完整拓扑 QC，避免资格方法与生产方法漂移。开发性 `qualification-handoff-plan` 的 Top-B、每 site 起点数、姿态簇代表规则、工具身份和并发数都会写入不可变计划；native 报告只在候选冻结后核对来源身份，不参与起点选择。包内默认每 site 保留 2 个起点；当前 P15/CK2 项目在 `configs/validation.toml` 中覆盖为 4 个，以覆盖两个主要姿态簇的 medoid 和确定性最远点。seed 多样性只在几何条件相同后参与排序，不再作为排除同 seed 其他构象的硬门槛。
每层仍要求至少 6/8 模型且覆盖 6 个 seed 通过，不能为保留 10 Å 包络而降低旧门槛。
报告同时检查受体 CA 保真、肽内部形变、site 质心、接触保持、S–S 距离、严重界面穿插和肽内非局部重原子穿插，并选择连续通过分层判据的最大候选阈值。Rosetta `fa_rep`、`omega` 和 `rama_prepro` 完整记录但不门控；旧的“总复合物分数/全部残基数”会被 300 多个受体残基稀释，不能代表局部环肽质量。
`min_models_for_selection=10` 只表示模型少于 10 个时不启动聚类，不参与 topology 科学资格判定。正式阶段二采样仍然只运行 CABS-dock；这里的全原子步骤只用于离线校准粗粒化筛选语义。
旧 SG 约束协议曾得到 `≤6 Å` 6/8、`6–7 Å` 6/8、`7–8 Å` 5/8；该结果没有迁移到当前 CA 约束协议。当前 `/5` 报告中 6/7/8/10 Å 四层均为 8/8 模型和 8/8 seed 通过，10–12 Å 外部对照为 7/8；因此 `max_reconstructable_disulfide_ca_distance_A=10.0` 与 `topology_calibration_status=qualified` 已由独立报告冻结。
第二个及后续 target 可以在计划时引用一个已通过的资格运行；计划会冻结并校验其 plan、sampling、pose、原生回收表和报告哈希，仅展开新 target 的 8 个先导任务。

`chemistry.ligand.ligand_id` 是运行时配体身份和化学记录目录名，当前值为 `p15`。
配体序列、二硫键、端基、His 微状态和净电荷都从该 section 读取；更换多肽时不需改 Python，但必须同步重新声明序列设计位点、guided 映射和相关科学门槛。化学记录位于 `data/chemistry/<ligand-id>/chemistry_record.json`。

`source_revision`、CABS 可执行文件 SHA-256 和影响约束、采样、TRAF 解码及过滤的关键源码哈希共同构成方法身份。项目不修改或补丁化外部 CABS 源码；运行前还会拒绝关键 CABS 文件相对声明提交存在未提交改动的环境。新环境只需安装声明版本，并确认官方 `--ca-rest-add` 接口可用。
不可变计划中的 Vela 摘要代表采样生产者身份。任务执行和 checkpoint 恢复必须与该摘要完全一致；只读分析允许在 plan/sampling schema 未改变时由修复后的分析器消费历史证据，报告分别写入 `sampling_software` 和 `analysis_software`。生产放行和共享控制只接受当前报告 schema 且 `analysis_software` 与当前源码一致的报告，不识别旧报告字段或旧 schema。
当前机器的 CABS 根目录是 `/home/glen/apps/cabsflex`，因此 `configs/discovery.toml` 的可执行文件和源码路径，以及 `pyproject.toml` 中仅供 Pyright 解析适配脚本导入的 `scripts` execution environment，均使用该绝对路径。更换机器时必须同时改成新机器的真实安装位置；`typings/CABS/` 只声明适配脚本实际调用的上游接口，不包含 CABS 实现，也不会参与运行时导入。

阶段二固定使用 CABS 粗粒化 `CA/SC` 输出，并固定传入原生 `--ca-rest-add`、`-A N` 和无已知 site 约束的随机起点。当前 Cys 端点约束为 5.5 Å、权重 1.0；这是粗粒化拓扑先验，不冒充真实 S-S 键。关闭全原子重建是阶段边界，不是隐藏参数：N/C 端电荷、C 端酰胺和 His7 微状态在阶段二的 CABS 表示中不会被逐原子区分，必须在阶段三全原子局部精修时恢复并复核。
当前 MC、REMC、过滤和 k-medoids 数值使用 CABS-dock 3.0.12 默认协议：
`20 × 50`快照/replica、`mc_steps=50`、10 replicas、`replicas_dtemp=0.5`、
`temperature=2.0→1.0`、每 replica 保留 100 个低能帧并生成 Top-10 medoid。Top-10 继续保存以便和公开方法比较，但正式候选从完整 10,000 帧 TRAF 先按 Cα 二硫拓扑可行性和 10 Å Cα 受体接触过滤，再按宽松 site 边界分组；两层均使用清单显式冻结的确定性、簇直径受限 leader clustering。每个任务按压缩后候选数、最低相互作用能和稳定身份最多保留 64 个 site，每个 site 内按肽 Cα RMSD 聚类后用同一规则最多保留 4 个 pose 簇；每簇保存 medoid 和最低 CABS 相互作用能代表。阶段三按跨 seed 支持度、簇内候选数和能量对姿态簇排序，在前 `ceil(起点预算/2)` 个簇间先轮转选择 medoid，再轮转选择距已选代表最远的构象。参数来源和不能直接搬用的 site 阈值见 [`research.md`](research.md)。

阶段三把实验结合态的用途分开保存：`receptor_only` 用于结合态受体的无位置约束复现，`pair_reference` 用于实验 site 和接触参考，只有 4IB5/Pc 这样的标准二硫环肽才生成 FlexPepDock 方法控制输入。9FBM/9FBI 的 FMP37 是非标准环肽模拟物，登记为 `site_reference_only`，不会因配置存在就被强行回对接。

局部恢复控制没有写死为 4IB5/Pc。`[[validation.local_controls]]` 通过 `bound_state_id` 引用任意登记为 `standard_cyclic_peptide` 的结合态，并分别声明 prepack seed、`seed_batches`、初始扰动、排名 score、恢复 RMSD score、聚类距离、Top cluster 窗口、簇的最小 seed 支持和跨批姿态一致性门槛。更换或增加验证结构时修改 `validation.toml`；包内默认参数不保存具体 PDB、配体或控制身份。

`validation qualification-refinement-plan` 是 `cross_receptor_pose_robustness_diagnostic` 开发入口，不是正式资格入口。调用方必须用可重复的 `--start-id` 显式列出已经通过 handoff 的任务 ID，并用 `--receptor-backbone-mode fixed|local_constrained` 明确冻结受体主链协议。3Q9X 起点相对 4IB5 native 的 RMSD 只能表示 cross-docking 鲁棒性；失败不否定 4IB5/Pc matched local control，也不强制启动 `lowres_preoptimize`。

`validation exploration-handoff-plan` 专门消费 `exploratory_discovery`，不会把它降格复制成正式 `main_discovery`。命令要求分别提供按冻结盲排名排序的 `--blind-candidate-id` 和 `--functional-candidate-id`；计划写入两条永久分离的candidate arm、`production_qualified=false`、48起点预算和事前晋级合同。全原子交接QC只判断化学、拓扑、重建保真、site保持和严重碰撞是否通过，不允许根据连续QC值或Rosetta分数重排candidate。当前合同要求每个来源受体site至少2个起点通过；深度精修只选择盲分支前2个合格candidate和功能分支前1个合格candidate，失败时只能在同一分支按原顺序递补。blind精修支持簇必须同时覆盖至少2个CABS起点、2个独立CABS source seed和2个独立Rosetta seed；同一CABS轨迹内的多个起点不再被当作独立全局采样证据。

`validation confirmation-handoff-plan` 只用于复核已经完成开发精修的单个候选，不修改候选原排名或生产资格。它从现有探索轨迹中为每个来源受体site选择4个不同CABS source seed，各seed只保留其排名最高姿态簇中的seed内medoid，并把选择方法、输入哈希和`production_qualified=false`写入不可变计划。该命令不重跑CABS；后续仍复用普通`handoff-run`、`refinement-plan`、`refinement-run`和`refinement-analyze`，但证据类别保持为`exploratory_source_seed_confirmation_*`。

`validation refinement-plan` 直接读取来源manifest，不接受手工候选过滤参数。对于`exploratory_discovery_handoff`，它会先复核handoff计划、来源证据、`production_qualified=false`和完整晋级合同，再按每个受体site的pass/fail数量执行冻结的同分支选择；连续QC值不参与重排。计划显式记录入选candidate、每个candidate/site的有效起点数、起点数、seed数、任务数和总decoy预算，并把探索身份传播到运行和分析产物。当前项目由44个通过起点确定性选择`ALPHA_C014`、`ALPHA_C016`和`ALPHA_C027`的23个起点，展开92个任务和11,776个decoy。

`validation.refinement.receptor_backbone_contact_A` 和 `receptor_backbone_sequence_padding` 只定义 `local_constrained` 的 native-free 自由度选择。当前协议以每个起始复合物自身的肽重原子为中心选择 6 Å 受体接触壳层，再沿受体序列向两侧扩展 2 个 pose 残基；MoveMap 保留全部受体侧链，只开放所选受体主链、全部肽主链/侧链和唯一对接 jump。FlexPepDock 的 `min_receptor_bb` 同时给受体 Cα 施加均值 0 Å、标准差 1 Å 的内部谐和坐标约束；全受体无约束主链最小化不属于当前协议。

`[[validation.guided_templates]]` 只允许引用登记为 `standard_cyclic_peptide` 的实验结合态，并且必须逐位列出与 P15 等长的 `ligand_positions`。加载配置时会同时核对映射长度、方向、范围和二硫键端点。当前 `4IB5_Pc_loop_to_P15` 将 P15 的 11 个残基映射到 Pc 的 Cys2–Cys12 环骨架；Python 不包含 `4IB5` 或位置 `2..12` 的特例。FMP37 没有标准肽化学合同，因此不能仅凭空间形状伪装成可线程化模板。

`[[validation.environment_references]]` 当前登记 `4MD7` assembly 1 和 `1JWH` assembly 1。参考项显式区分用于对齐的 CK2α 链、CK2β 链和另一催化链，并记录构建体限制。`4MD7_E`、`1JWH_A` 在 `[[receptors]]` 中只有 `full_enzyme_environment` 角色，不会进入 `blind_discovery` 或 `bound_state_blind_replication`。CK2α′姿态映射到这些 CK2α 实验排布时，输出固定标记为同源布局模型，而不是 CK2α′实验全酶结构。

结合态全表面复现也不维护另一份受体或采样参数。`[[receptors]]` 中带 `bound_state_blind_replication` 角色且被 `[[validation.bound_states]]` 引用的条目，会自动展开为“结合态 `state_id` × `discovery.seeds`”任务；受体文件来自对应 `receptor_only.cif`，全表面方法和分析阈值直接复用已资格化的 `discovery` 配置。这样更换验证结构只需要修改受体与结合态登记，不需要修改 Python，也不会产生两套可能漂移的 CABS-dock 参数。运行清单是执行器的任务权威来源，主发现与结合态复现的证据类别和输出目录始终分开。

正式候选链路同样没有隐藏参数：`validation.handoff` 决定常规交付的每个受体site起点数；`validation.funnel`另行冻结探索漏斗的候选、source起点和逐层晋级预算。两者都复用`discovery.cabsdock.pose_clustering_rmsd_A`，在受体对齐后按有序肽CA RMSD建立姿态家族。Stage 3A-0要求筛选起点来自同一个4.0 Å家族中的不同CABS source seed，先排除既有高等级阴性，再按最弱受体的家族source支持、总支持、阶段二等级内排名和candidate ID确定性排序。`validation.refinement.seed_batch_sizes=[1,1,2]`只切分`validation.seeds`这一处声明的四条随机流，依次供screening、confirmation和deep confirmation增量使用；后层复用前层任务，不重新计算相同`start×seed`。`validation.analysis`声明界面、碰撞、受体和site保持门槛；Stage 3A的跨source命中单列为`cross_source_screening_hit`，不会因只有一条Rosetta随机流而冒充正式`supported`。

Stage 3B要求同一原子簇至少覆盖两个source起点、两个Rosetta seed以及四种`source×Rosetta seed`任务单元中的三种。Stage 3C再要求至少两个source起点分别获得至少两个Rosetta seed支持。Route B的跨受体兼容同时要求受体对齐后的P15骨架RMSD不高于3.0 Å、受体接触指纹Jaccard不低于0.4；两个条件均满足才可把两个受体上的簇视为同一原子假设。这些阈值在Stage 3A结果产生前冻结，不能看完结果后调整。
精修分析把采样软件身份和当前分析软件身份分别记录。单个decoy的二硫键、端基或几何科学QC失败会写入`qc_failures`并从聚类中排除，不会中断整批分析；文件缺失、摘要不一致或无法建立必要几何仍属于技术错误。界面接触使用Gemmi空间邻域索引计算，构象聚类使用结果等价的距离缓存complete-linkage实现，避免模型数量增加时重复计算全部原子对和簇对距离。
阶段二候选引用 `trajectory_candidates.pdb` 时，handoff 计划、任务结果和重建器都会冻结并复核 `model_index`；阶段三只提取该模型，不能把多模型文件静默当作第一帧。

blind handoff 与 guided handoff 共用同一个局部精修执行器，但来源 manifest 决定 `evidence_category` 和 `known_site_information_used`，命令行不能覆盖这两个身份。全酶映射复用 `validation.interface_contact_A` 报告原始接触数，同时保留最小重原子距离；当前没有经过校准的全酶“碰撞/竞争”判定阈值，因此代码不会据此自动宣布 CK2β 竞争。

`validation.seeds` 已冻结为 `120623..120626`；`validation.qualification_status="qualified"` 引用当前源码重析后的 4IB5/Pc matched-control 报告，`validation.analysis` 也已冻结。这只放行“已有候选起点时可使用当前局部精修协议”；正式 P15 仍需先通过阶段二 site-discovery 交付。代码把“起点 × seed”展开为独立单进程任务，`decoys_per_seed=128` 直接传给 Rosetta `-nstruct`，即每个起始姿态 `4 seeds × 128 = 512` 个 decoy。

4IB5/Pc 资格控制使用两个完整且随机流不重叠的 `4 × 128` 批次。
每批的 512 个 decoy 先合并，然后按 `reweighted_sc` 统一排分，在实验受体对齐后按肽主链 RMSD 统一聚类。近天然簇必须进入 Top-10 clusters 且由至少 2/4 个 seed 共同支持；两批都必须分别达到采样成功和排序成功。最低能量近天然成员只证明排序能力，跨批次姿态一致性必须比较各自主簇的几何 medoid，二者不能混用。当前控制的 medoid RMSD 为 0.505 Å并通过 2.0 Å门槛。

阶段四不存在隐藏的 P5/G9 固定规则或“不得新增 His”规则。当前 `design.sequence.mutable_positions=2..10`，Cys1/Cys11 由二硫键合同固定，允许集合排除额外 Cys。新 His 使用 `candidate_histidine_state=HIE`，运行命令会按每条实际序列动态固定对应 pose 编号，而不是依赖 WT His7 的硬编码位置。

`design.analysis` 只判断低成本界面初筛；`design.finalists` 独立控制昂贵柔性复核。两层必须各自冻结 seed 和经独立控制校准的门槛，不能在看到候选后用同一组数字反复调到通过。`max_candidates`、`max_md_candidates` 是资源截断，因预算未推进的候选会记录为 `resource_deferred`，不能写成科学失败。序列疏水性、芳香性、Met 和简单化学风险 motif 是输出事实，不是当前未校准的配置门槛。

`design.combination` 首轮只生成并实际计算双突变。三突变及以上只能从已完成柔性复核、被独立选为迭代父序列的双突变逐代产生。
初筛中的相同 WT 输入按 `wt_context_id` 去重；`design.finalists.max_flexibility_required_candidates` 为 Gly/Pro、电荷类别、小侧链到大侧链、二硫键邻位和相邻多突变等固定骨架风险保留独立柔性复核预算。finalist 先按肽骨架 complete-linkage 聚类，并结合骨架 RMSD 与受体接触相似度匹配 WT/突变体 pose；只有 `matched_pose` 可以产生配对 delta。
`design.iteration` 默认允许最多 4 个同代父序列、相对 WT 累计最多 5 个突变、每代最多 256 条初筛候选。每代仍只改变父序列的一个位置；没有单独的“每代突变数”配置，因为一步邻域是当前算法合同。完整邻域都会写入候选表，超出预算的条目只标记 `resource_deferred`。直接父序列、父代编辑、代次和祖先 ID 随计划冻结并在重读时核对。

阶段二 CABS、阶段三 Rosetta、阶段四 RosettaScripts/FlexPepDock 当前都使用 8 个 worker，在独立受体、seed 或起点任务层并行；每个外部进程独占结果目录并在进程内部串行完成单个任务。当前不使用 MPI rank 分发：真实试跑表明本机 OpenMPI/PRRTE 会在 Rosetta 已完成后错误返回非零状态，而同一 Rosetta 构建直接运行可正确完成并退出；任务层并行也避免多个 rank 争写同一目录。

`design.objective="single_supported_target"` 表示每个阶段四运行只接受一个具有阶段三支持的亚型及其多个正向模板，不设置另一个亚型为反向模板。CK2α 与 CK2α′可以分别运行，任一轨道通过即可进入后续候选池；这表示计算准入，不表示已经证明治疗有效。

## 递归覆盖示例

项目只需写出要改变的默认字段，其余字段继续来自包内默认值：

```toml
[audit.crystal_contacts]
distance_A = 5.0

[preparation.altloc]
preferred_label = "B"

[discovery.ensemble]
min_receptors_per_target = 3
```

相对路径始终相对于包含它的外部配置文件解析。当前环境变量覆盖包括 `VELA_DATA_DIR`、`VELA_OUTPUTS_DIR`、`VELA_DOWNLOAD_RETRIES` 和 `VELA_DOWNLOAD_TIMEOUT_SECONDS`。

## 有意保留在代码中的固定合同

以下内容不是可调实验参数，不应为了“零硬编码”移入配置：

- manifest、TSV 和目录协议的文件名及 schema 版本；
- CABS-dock `filtering_mode=each`、粗粒化 `FMS` 输出和关闭全原子重建；这些共同定义当前已经验证的适配器合同，而不是可在同一方法 ID 下任意改变的参数；
- SHA-256 格式、PDB ID 语法、数值合法范围和必需字段；
- `planned`、`completed`、`passed` 等稳定状态枚举；
- 同一 seed 不能重复增加支持度、两个亚型不能混聚、坐标系不一致不能直接比较；
- 位点距离的归一化门槛 `1.0`，它是输入距离阈值的定义而非额外参数；
- altloc 按最高平均占有率选择，这是当前唯一实现的制备算法，不伪装成可调策略；
- 原子写入分块、JSON 缩进和报告小数位等不影响研究含义的实现细节。

如果某个固定值以后确实需要因研究方案、外部工具、计算资源或受体集合而变化，应先明确其语义和验证范围，再加入相应 TOML section；不能仅为了减少代码字面量而制造无实际变化依据的配置项。
