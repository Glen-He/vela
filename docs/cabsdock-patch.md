# CABS-dock 源码补丁与环境复现

## 目的与边界

阶段二只使用 CABS-dock 做二硫环肽的粗粒化全表面盲对接，不使用 Rosetta 生成起始
构象或补救拓扑失败。CABS 3.0.12 的 `DockTask` 已经能够处理
`disulfide_bonds`、生成 `SG` 约束并把约束交给模拟核心，但 `dock_dict` 没有把已有的
`-F/--disulfide-bonds` 选项暴露给 `CABSdock` 命令。因此本项目采用一个最小补丁修复
CLI 入口，不修改采样算法、势能函数或过滤规则。

补丁只包含两项有行为意义的变化：

- 在 `CABS/io/config.json` 的 CABSdock restraint 选项中加入
  `disulfide-bonds`；
- 在 `tests/test_cli.py` 增加一条解析回归测试，确认
  `-F 1:PEP1 6:PEP1` 被保存为一个二硫键端点对。

补丁同时把 `CABS/io/config.json` 的文件末尾规范为换行结尾；这不改变运行行为，但使
当前源码 diff 与保存的补丁完全一致。

项目补丁文件是
[`patches/cabsdock-3.0.12-disulfide-cli.patch`](../patches/cabsdock-3.0.12-disulfide-cli.patch)。
它不包含本机 CABS 仓库内已有的 `cg2all` 改动，也不要求阶段二启用全原子重建。

## 已验证基线

| 项目 | 值 |
|---|---|
| 上游仓库 | `https://github.com/LCBio/CABSflex_standalone.git` |
| Git 提交 | `36ec10e681d0b6c5c991101bc0bdfc6e224c5e0b` |
| CABS 版本 | `3.0.12` |
| 当前源码目录 | `/home/glen/apps/cabsflex` |
| 当前可执行文件 | `/home/glen/apps/cabsflex/.venv-cabs/bin/CABSdock` |

正式运行会同时检查可执行文件实际导入的源码目录、Git 提交、项目补丁是否已经应用、
CLI 是否显示 `--disulfide-bonds`，并把补丁 SHA-256 写入运行清单。只看到
`CABS version 3.0.12` 不足以证明环境正确。

## 在新环境中应用

下面的路径变量只用于说明；在新机器上替换成实际绝对路径：

```bash
CABS_REPO=/path/to/cabsflex
VELA_REPO=/path/to/vela
```

先准备与本项目一致的源码基线。目标目录必须是专用于 CABS 的源码仓库；若已有未提交
修改，应先审计，不能通过 reset 或覆盖来套用本说明。

```bash
git clone https://github.com/LCBio/CABSflex_standalone.git "$CABS_REPO"
git -C "$CABS_REPO" checkout 36ec10e681d0b6c5c991101bc0bdfc6e224c5e0b
git -C "$CABS_REPO" rev-parse HEAD
```

应用项目补丁：

```bash
git -C "$CABS_REPO" apply \
  "$VELA_REPO/patches/cabsdock-3.0.12-disulfide-cli.patch"
```

若 CABS 尚未安装，使用 uv 创建主环境并采用 editable 安装，使源码补丁直接生效：

```bash
uv venv "$CABS_REPO/.venv-cabs" --python 3.10
uv pip install --python "$CABS_REPO/.venv-cabs/bin/python" -e "$CABS_REPO"
```

CABS 的 Fortran 编译器和其他运行依赖仍按上游安装要求准备；本补丁不替代完整的 CABS
安装步骤。阶段二固定传入 `-A N`，因此不依赖 MODELLER、cg2all 或 Rosetta 做全原子
重建。

最后修改 [`configs/discovery.toml`](../configs/discovery.toml) 中的以下项目覆盖：

```toml
[discovery.cabsdock]
executable = "/path/to/cabsflex/.venv-cabs/bin/CABSdock"
source_dir = "/path/to/cabsflex"
source_revision = "36ec10e681d0b6c5c991101bc0bdfc6e224c5e0b"
patch_file = "../patches/cabsdock-3.0.12-disulfide-cli.patch"
```

## 验证补丁

首先验证补丁确实已经应用。`--reverse --check` 只检查能否撤销，不会修改文件：

```bash
git -C "$CABS_REPO" apply --reverse --check \
  "$VELA_REPO/patches/cabsdock-3.0.12-disulfide-cli.patch"
"$CABS_REPO/.venv-cabs/bin/CABSdock" --help | rg -- '--disulfide-bonds'
```

运行新增的 CLI 回归测试：

```bash
uv run --no-project --isolated \
  --python "$CABS_REPO/.venv-cabs/bin/python" \
  --with pytest \
  pytest -c /dev/null -p no:cacheprovider -q \
  "$CABS_REPO/tests/test_cli.py::TestCLI::test_dock_disulfide_capture"
```

预期结果是 `1 passed`。然后在已经完成阶段一受体制备的 Vela 项目中运行最小原生
二硫键 smoke：

```bash
CABS_SMOKE_DIR="$(mktemp -d)"
"$CABS_REPO/.venv-cabs/bin/CABSdock" \
  -i "$VELA_REPO/data/receptors/prepared/3Q04_A.cif" \
  -p CWMSPRHLGTC:CCCCCCCCCCC \
  -F 1:PEP1 11:PEP1 \
  -a 1 -y 1 -s 1 -r 1 -D 0.1 -t 2.0 1.0 \
  -n 1 -k 1 --clustering-iterations 1 -z 104729 \
  -A N -o FM -C --json-output --restraints-output --no-progress-bar \
  -w "$CABS_SMOKE_DIR"
rg -n 'disulfide-bonds.*1:PEP1 11:PEP1' "$CABS_SMOKE_DIR/config.ini"
rg -n '1:PEP1 11:PEP1 2\.0000 0\.0000 1\.00 SG' \
  "$CABS_SMOKE_DIR/output_data/restraints.txt"
```

两次 `rg` 都必须命中。smoke 只验证接口和约束生成，不能代替 `docs/plan.md` 中的
1,000 帧拓扑、10 medoid、4IB5/Pc 召回和多 seed 资格门槛。

回到 Vela 项目后运行：

```bash
cd "$VELA_REPO"
uv run vela config check
uv run ruff format --check src tests
uv run ruff check src tests
uv run pyright
uv run pytest
```

## 更新上游版本时

若更换 CABS 提交或版本，不得只把 `source_revision` 改成新值继续运行。应依次执行：

1. 在新的干净源码基线上尝试应用补丁；
2. 若补丁不再适用，重新核查 CABSdock 是否已经原生暴露二硫键选项；
3. 需要新补丁时只修改对应入口并重新增加回归测试；
4. 更新项目补丁文件、Git 提交和相关哈希；
5. 重跑 CLI 测试、原生二硫键 smoke 和完整方法资格验证；
6. 使用新的 method ID，不能把改变后的实现冒充当前方法。

## 仅撤销本补丁

需要移除本项目补丁时，可在确认目标仓库仍匹配补丁内容后执行：

```bash
git -C "$CABS_REPO" apply --reverse \
  "$VELA_REPO/patches/cabsdock-3.0.12-disulfide-cli.patch"
```

该命令只撤销本补丁中的两处变化，不会触碰其他用户修改。撤销后 Vela 的阶段二工具
检查应失败，这是预期行为。
