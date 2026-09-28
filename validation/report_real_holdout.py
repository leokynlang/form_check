"""Render the saved evaluation as a Chinese report and static research plots."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "runs" / "real_engineering_holdout"
LABELS = {"power_law": "幂律", "log_poly_2": "对数二次式", "log_poly_3": "对数三次式",
          "log_poly_9": "对数九次式（压力对照）", "group_cv_poly": "分组 CV 选多项式",
          "discover_auto": "程序自动选型", "nist_form_refit": "NIST 形式／训练集重拟合"}
SPLITS = {"random_rows": "随机行留出", "heldout_groups": "整组工况留出",
          "high_distance": "更大距离", "middle_to_tails": "中间范围→两端",
          "late_time": "更长老化时间", "high_temperature": "更高温度"}
COLORS = {"power_law": "#576579", "log_poly_2": "#007f79", "group_cv_poly": "#007f79",
          "discover_auto": "#d06b32", "log_poly_9": "#9a4b73", "nist_form_refit": "#5a64a8"}


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join(["---"]*len(headers)) + " |",
                      *["| " + " | ".join(map(str, row)) + " |" for row in rows]])


def main():
    s = pd.read_csv(OUT / "summary.csv")
    m = pd.read_csv(OUT / "metrics.csv")
    c = pd.read_csv(OUT / "checks.csv")
    p = pd.read_csv(OUT / "predictions.csv")
    audit = json.loads((OUT / "audit.json").read_text())
    plt.rcParams.update({"font.family": "Microsoft YaHei", "font.size": 10,
                         "axes.unicode_minus": False, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "#f7f8fa",
                         "axes.facecolor": "white", "savefig.facecolor": "#f7f8fa"})

    fig, axs = plt.subplots(2, 2, figsize=(14, 10), layout="constrained")
    fig.suptitle("真实工程数据：能发现形式问题，也暴露选型与外推边界", fontsize=19, fontweight="bold")
    ax = axs[0, 0]
    raw = pd.read_csv(OUT / "Chwirut1_clean.csv")
    one = p[(p.dataset == "Chwirut1") & (p.split == "high_distance")]
    ax.axvspan(4.125, 6.1, color="#fbe7db", alpha=.8, label="留出区间")
    for partition, color, label in [("train", "#8794a4", "训练观测"), ("test", "#b7443d", "留出观测")]:
        points = one[(one.model == "power_law") & (one.partition == partition)].merge(raw, on="row_id")
        ax.scatter(points.distance, np.exp(points.observed_log), s=13, color=color, alpha=.6, label=label)
    for model in ["power_law", "group_cv_poly", "discover_auto"]:
        points = one[one.model == model].merge(raw[["row_id", "distance"]], on="row_id")
        curve = points.groupby("distance").predicted_log.mean().sort_index()
        ax.plot(curve.index, np.exp(curve), color=COLORS[model], lw=2, label=LABELS[model])
    ax.set(xlabel="金属距离（保留数据源刻度）", ylabel="超声响应（对数坐标）",
           title="A  超声标定：训练区间好看，外推仍可偏离", yscale="log", xlim=(.4, 6.1))
    ax.legend(fontsize=8, ncol=2, frameon=False)

    ax = axs[0, 1]
    names = ["power_law", "group_cv_poly", "discover_auto", "log_poly_9"]
    for offset, split, color in [(-.17, "random_rows", "#a3b4c9"), (.17, "heldout_groups", "#007f79")]:
        rows = [s[(s.dataset == "Chwirut1") & (s.split == split) & (s.model == name)].iloc[0] for name in names]
        vals = np.array([r.test_rmse_log_median for r in rows])
        low = vals-np.array([r.test_rmse_log_q25 for r in rows])
        high = np.array([r.test_rmse_log_q75 for r in rows])-vals
        ax.barh(np.arange(len(names))+offset, vals, height=.3, color=color,
                xerr=np.stack([low, high]), error_kw={"elinewidth": 1, "capsize": 3}, label=SPLITS[split])
    ax.set_yticks(np.arange(len(names)), [LABELS[n].replace("（压力对照）", "") for n in names])
    ax.invert_yaxis()
    ax.set(xlabel="测试 log-RMSE：中位数与四分位范围（越低越好）",
           title="B  超声标定：九次式对新距离更不稳定")
    ax.legend(fontsize=9, frameon=False)

    ax = axs[1, 0]
    rows = [m[(m.dataset == "Chwirut1") & (m.split == "high_distance") & (m.model == name)].iloc[0] for name in names]
    vals = [r.test_rmse_log for r in rows]
    ys = np.arange(len(names))
    ax.barh(ys, vals, color=[COLORS[n] for n in names], height=.58)
    ax.set_yticks(ys, [f"{LABELS[n].replace('（压力对照）', '')}\n训练 R²={r.train_r2_log:.3f}" for n, r in zip(names, rows)])
    ax.invert_yaxis()
    ax.set(xscale="log", xlim=(.1, 17), xlabel="范围外测试 log-RMSE（横轴为对数刻度）",
           title="C  更大距离留出：训练 R² 提升不保证泛化")
    for y, val in zip(ys, vals):
        ax.text(val*1.07, y, f"{val:.3f}", va="center", fontsize=10)

    ax = axs[1, 1]
    keys = ["fc_alarm_native", "rej_fc", "rej_fitted_hc3"]
    labels = ["原生 form_check\n|t|>3", "form_check 包装层\np<0.05", "RESET fitted HC3\np<0.05"]
    for offset, split, color in [(-.17, "random_rows", "#a3b4c9"), (.17, "heldout_groups", "#007f79")]:
        rows = c[(c.dataset == "Nelson") & (c.split == split)]
        vals = np.array([int(rows[k].sum()) for k in keys])
        ax.bar(np.arange(3)+offset, vals, width=.3, color=color, label=SPLITS[split])
        for x, val in zip(np.arange(3)+offset, vals):
            ax.text(x, val+.35, f"{val}/20", ha="center", fontsize=9)
    ax.set_xticks(np.arange(3), labels)
    ax.set(ylim=(0, 23), ylabel="训练集报警次数（同一数据集重复拆分）",
           title="D  绝缘老化：单方向曲率检查可能沉默")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.savefig(OUT / "overview.png", dpi=150)
    fig.savefig(OUT / "overview.svg")
    plt.close(fig)

    # Full descriptive tables, including failed comparisons.
    summaries = []
    for _, row in s.iterrows():
        val = "失败" if row.successes == 0 else f"{row.test_rmse_log_median:.4f} [{row.test_rmse_log_q25:.4f}, {row.test_rmse_log_q75:.4f}]"
        summaries.append([row.dataset, SPLITS[row.split], LABELS[row.model],
                          f"{int(row.successes)}/{int(row.runs)}", val])
    alarm_rows = []
    for (dataset, split), rows in c.groupby(["dataset", "split"], sort=False):
        cells = []
        for k, fail in [("fc_alarm_native", "fail_fc"), ("rej_fc", "fail_fc"),
                        ("rej_fitted_f", "fail_fitted_f"), ("rej_fitted_hc3", "fail_fitted_hc3"),
                        ("rej_exog_hc3", "fail_exog_hc3")]:
            nfail = int(rows[fail].notna().sum())
            cells.append(f"{int(rows[k].sum())}/{len(rows)-nfail}" + (f"（{nfail} 失败）" if nfail else ""))
        alarm_rows.append([dataset, SPLITS[split], *cells])
    domain = []
    for _, row in c[~c.split.isin(["random_rows", "heldout_groups"])].iterrows():
        domain.append([row.dataset, SPLITS[row.split], f"{row.fc_t_native:.3f}",
                       f"{row.p_fc:.3g}", f"{row.p_wild_sensitivity:.3g}", f"{row.p_fitted_hc3:.3g}"])
    auto = []
    for _, row in m[(m.model == "discover_auto") & ~m.split.isin(["random_rows", "heldout_groups"])].iterrows():
        auto.append([row.dataset, SPLITS[row.split], row.selected_family,
                     f"{row.inner_cv_rmse_log:.4f}", f"{row.test_rmse_log:.4f}",
                     f"{row.test_rmse_log/row.inner_cv_rmse_log:.2f}×"])
    primary = []
    for split in ["random_rows", "heldout_groups", "high_distance"]:
        records = []
        for model in ["power_law", "group_cv_poly", "discover_auto"]:
            row = s[(s.dataset == "Chwirut1") & (s.split == split) & (s.model == model)].iloc[0]
            records.append(f"{row.test_rmse_log_median:.4f}")
        primary.append([SPLITS[split], *records])
    paired = []
    for dataset in ["Chwirut1", "Nelson"]:
        for split in ["random_rows", "heldout_groups"]:
            rows = m[(m.dataset == dataset) & (m.split == split)]
            wide = rows.pivot(index="seed", columns="model", values="test_rmse_log")
            for model in ["discover_auto", "group_cv_poly", "nist_form_refit", "log_poly_9"]:
                if model in wide:
                    ratio = wide[model]/wide.power_law
                    paired.append([dataset, SPLITS[split], LABELS[model],
                                   f"{int((ratio<1).sum())}/{len(ratio)}", f"{ratio.median():.3f}"])

    report = f"""# form_check：真实工程数据的优势、过拟合风险与边界

日期：2026-09-05。使用两组官方实测数据，84 个训练／测试拆分、546 次模型拟合。
这是按任务目的选择的探索性案例验证，不是随机抽样的大规模优越性评测。
产品 `dimreflex/` 和 `common.py` 未修改，代码哈希核验通过。

## 主要判断

**它有价值，但价值目前主要是提示幂律可能失配，不能把自动选出的公式视为外推可靠。**
超声数据提供正面案例；老化数据提供单方向检查的边界；同一组超声数据又直接暴露了自动选型和复杂模型的泛化风险。
本次没有证明它优于 RESET：超声案例中各检验都报警，老化案例中 RESET 更稳定地报警。

![验证总览](overview.png)

## 数据与独立验证方式

- **Chwirut1：214 条超声标定实测数据**，输入为金属距离、输出为超声响应。
  来源：[NIST 数据说明](https://www.itl.nist.gov/div898/strd/nls/data/LINKS/i-chwirut1.shtml)、
  [原始数据](https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Chwirut1.dat)。
- **Nelson：128 条绝缘性能退化实测数据**，输入为老化时间、温度，输出为击穿强度。
  来源：[NIST 数据说明](https://www.itl.nist.gov/div898/strd/nls/data/LINKS/i-nelson.shtml)、
  [原始数据](https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Nelson.dat)。
  温度从摄氏度转换为开尔文后才进入幂律／对数模型。

两类随机拆分各重复 20 个固定种子：随机留出 25% 的行；按相同输入组合整组留出 25% 工况。
相同距离／相同时间温度组合的重复观测在分组拆分中不会跨训练和测试集。
这些重复拆分来自同一数据集，不是 20 份独立数据，频率不应当解释成检验功效或误报率。

范围外验证预先固定为：超声距离 **0.5–4.0 训练、4.25–6.0 测试**；另做中间范围 1.25–4.25 训练、两端测试。
老化时间 **1–32 周训练、48–64 周测试**；另做 **180–250°C 训练、275°C 测试**。
拆分仅依赖输入，报警、拟合、预处理、内层选型均只使用训练集。

主指标为 **log-RMSE，越低越好**。它是自然对数尺度上的均方根误差，不是百分比。
同时保存 log-R²、相对误差中位数与逐条预测。所有模型都拟合 log 响应，保持损失口径一致。
NIST 给出的函数形式只作参照，其参数全部在训练集重拟合；官方认证参数仅在计算完成后的解析审计中使用。
NIST 形式不是已知真值，本次拟合的目标也不是复现 Chwirut 的原尺度最小二乘认证解。

## 1. 能体现优势的地方：高拟合度没有挡住报警

超声数据中，随机行拆分的幂律训练 log-R² 中位数为 **0.9530**；整组留出时为 **0.9510**。
两种拆分中原生 `form_check` 均 **20/20 报警**。在更大距离留出中，训练 log-R² 仍为 **0.9537**，
原生 t={float(c[(c.dataset=='Chwirut1') & (c.split=='high_distance')].iloc[0].fc_t_native):.3f}，包装层 p=0.001。
这说明它能在看似不错的拟合后继续提示形式问题。

下表均为测试 log-RMSE；随机拆分取 20 次中位数，范围外案例只做一个固定拆分。

{table(['超声验证方式', '幂律', '分组 CV 选多项式', '原程序 discover_auto'], primary)}

随机行测试中，自动选型的误差中位数比幂律低约 **21%**；整组工况测试中仅低约 **6%**。
这是“报警后尝试其他形式可能有帮助”的证据，但并不意味着每个拆分都改善，也不是检验本身自动提高预测精度。
分组 CV 多项式是本次新增实验对照，不能记成产品现有能力。

## 2. 最值得注意的产品风险：自动选型的内部 CV 不能代表范围外误差

更大距离案例中，原程序选中了 `rational_log`，训练 log-R² 提高到 **0.9709**，
但测试 log-RMSE 从幂律的 **0.3341** 增到 **0.4543**，约恶化 **36%**。
它的内部 CV 误差只有 **0.1158**。这次升级在训练范围内看起来合理，在范围外却失败。

{table(['数据', '留出方式', '原程序所选族', '内部 CV log-RMSE', '实际留出 log-RMSE', '比值'], auto)}

这些比值比较的是不同预测任务的误差，说明适用范围差异，不是严格的无偏过拟合估计。
随机行交叉验证能评估相近输入范围内的预测，却不能验证没有见过的距离或温度。
分组 CV 也不能自动保证外推：例如其更高温度测试仍很差。

## 3. 明确的复杂度压力案例：训练更好，新工况更差

预先加入九次对数多项式作为刻意灵活的压力对照。它不是产品内置模型，也不用于宣传产品成绩。
在超声整组留出中，它的训练 log-R² 中位数达到 **0.9695**，优于幂律的 **0.9510**；
但测试 log-RMSE 中位数为 **0.2394**，差于幂律的 **0.1808**，且四分位范围扩大到 **[0.1886, 1.2139]**。
在更大距离测试中误差达到 **9.0178**。中间范围训练、两端测试时更出现原尺度溢出，已保留并标记。
它说明小幅训练拟合改善可以伴随非常不稳定的泛化；不能据此声称 form_check 会检查任意复杂模型的过拟合。
当前 form_check 检查的是幂律形式，复杂模型的风险仍需独立验证。

## 4. 不利结果：绝缘老化数据暴露检查盲区

Nelson 随机行拆分中，原生 3σ 报警为 **0/20**，包装成 p<0.05 后为 **12/20**；
整组拆分中分别为 **4/20** 和 **14/20**。RESET fitted F 与 fitted HC3 在两种拆分中均为 **20/20**。
两种 form_check 阈值不同，必须分列，不能把包装层成绩说成原生报警能力。

尤其是 **1–32 周训练、48–64 周测试**：原生不报警，包装层 p=**0.085**，
但幂律测试 log-RMSE 达到 **0.8212**，三次对数模型为 **0.3028**。
这是一例“没有报警也可能存在很大预测风险”。
当前实现只检查多特征平方项平均后的一个方向；对方向、交互和输入覆盖的敏感性是可能解释，
本次没有单独隔离其因果贡献。幂律本身在此例的训练拟合也不高，因此它主要检验边界，不是高 R² 漏检案例。

{table(['数据', '方式', '原生 3σ', '包装层 .05', 'RESET F .05', 'RESET fitted HC3 .05', 'RESET exog HC3 .05'], alarm_rows)}

“次数/可运行次数”是稳定性描述。高温留出的 exog HC3 发生秩亏，表中按失败报告，不能算成未报警。
真实数据没有已知条件均值，不能由这些频次计算真实功效或误报率。

## 5. 残差零分布敏感性

新增的实验性 HC3 残差 wild bootstrap 只对同一个曲率方向做敏感性检查，结果不属于产品功能。
它保留了超声范围外案例的报警，也没有改变老化长时间案例的沉默。
这能检查结论对一种替代重采样方案是否敏感，不能证明噪声独立或全面稳健。

{table(['数据', '方式', '原生 t', '包装层 p', '实验 wild p', 'RESET fitted HC3 p'], domain)}

## 6. 当前证据支持的改进方向

1. 在自动选型接口中增加用户可指定的工况分组和范围留出验证；分别报告区间内与范围外成绩。
2. 将“未发现幂律失配”与“模型可用于外推”明确分开，不把低内部 CV 误差当作范围外通行证。
3. 多变量检查应增加独立方向、交互方向或与 RESET 组合；后续必须在新数据上验证，避免围绕本次结果调参。
4. 超出训练包络时使用已有监视能力提示，但需要另行验证监视输出与预测入口是否形成完整流程。

本次只验证形式检查和模型选型，未验证完整量纲推导、生产监视或真实工程控制。
两套数据经过目的性选择且已知有非线性结构，不能据此给出整体准确率、误报率或普遍优越性结论。
重复观测的实验批次／试件标识不足，分组只能按已提供的输入组合进行，不能保证独立实验室泛化。
随机行和整组拆分的测试难度及部分训练样本数不同，二者差距不能单独归因于数据泄漏。

## 7. 全部模型成绩与配对稳定性

表内为 log-RMSE 中位数 [Q25,Q75]；单次范围留出的三个数相同，不是置信区间。

{table(['数据', '方式', '模型', '成功/尝试', '测试 log-RMSE'], summaries)}

下面使用相同拆分进行配对比较。比值小于 1 表示优于幂律；这是重复拆分的描述性统计，不做独立样本显著性推断。

{table(['数据', '方式', '对照模型', '优于幂律次数', '配对误差比中位数'], paired)}

## 8. 审计与复现

- 原始文件、来源 URL 与 SHA256：`raw/`、`manifest.json`。
- 输入清洗表：`Chwirut1_clean.csv`、`Nelson_clean.csv`；本次未删除观测。
- 原始产品文件 SHA256 前后完全一致；原生 form_check／autodiscovery 相关测试 **17 passed**。
- 独立审计复算了 **{audit['prediction_partitions_recomputed']}** 个预测分区的 log-RMSE，
  检查全部 **84** 个拆分无行重叠、分组拆分无工况重叠，并验证实验幂律与产品 OLS 预测一致。
- 数据解析经 NIST 认证参数对应 RSS 交叉核对：Chwirut1 **{audit['source_rss_chwirut1']:.10f}**，
  Nelson **{audit['source_rss_nelson']:.12f}**。这些参数不参与训练和选型。
- 一次拟合失败：Nelson 高温留出的三次多项式只有三个训练温度层，设计矩阵秩为 9/10。
  一个压力对照出现原尺度相对误差溢出；log 预测与 log-RMSE 仍完整保存，未裁剪或掩盖。
- 首次计算发现新增参照拟合器的初始点超出设定边界；修正参数边界后全部重算。
  未据结果改变数据、拆分、模型集合或检验阈值。

在仓库根目录运行（原始数据已缓存）：

```powershell
python validation/real_engineering_holdout.py
python validation/audit_real_holdout.py
python validation/report_real_holdout.py
python -m pytest tests/test_form_check.py tests/test_autodiscovery.py -q
```

实验计划：`validation/real_data_protocol.md`。
逐拆分指标：`metrics.csv`；检验明细：`checks.csv`；逐条预测：`predictions.csv`；
样本拆分：`splits.json`；模型参数和内层 CV：`models.json`；审计：`audit.json`。
"""
    (OUT / "report.md").write_text(report, encoding="utf-8")
    print("Created report.md, overview.png, overview.svg")


if __name__ == "__main__":
    main()
