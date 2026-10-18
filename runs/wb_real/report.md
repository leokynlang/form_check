# form_check vs RESET — 世界银行真实宏观数据

- 数据: World Bank Open Data API, 跨国横截面, 年份 ['2019', '2021', '2022', '2023'], 剔除聚合体, 每格 5 个排列种子(form_check 经验 p 的稳定性)
- 统一 α=0.05; 方法 A=fc(排列经验 p), B=RESET fitted-F, C=RESET fitted-HC3, D=RESET exog-HC3

| 年份 | 关系 | n | fc p (中位) | fitted-F p | fitted-HC3 p | exog-HC3 p | 一致结论 |
|---|---|---|---|---|---|---|---|
| 2019 | co2pc__gdppc | 189 | 0.0050 | 0.0060 | 0.0004 | 0.0004 | REJ / REJ / REJ / REJ |
| 2019 | gdp__pop | 199 | 0.7050 | 0.3334 | 0.2990 | 0.2990 | no-rej / no-rej / no-rej / no-rej |
| 2019 | infmort__gdppc | 188 | 0.0030 | 0.0018 | 0.0039 | 0.0039 | REJ / REJ / REJ / REJ |
| 2019 | lifeexp__gdppc | 199 | 0.0350 | 0.0520 | 0.3408 | 0.3409 | REJ / no-rej / no-rej / no-rej |
| 2019 | lifeexp__gdppc_health | 186 | 0.2940 | 0.0806 | 0.6249 | 0.8244 | no-rej / no-rej / no-rej / no-rej |
| 2019 | netusers__gdppc | 182 | 0.0010 | 0.0000 | 0.0000 | 0.0000 | REJ / REJ / REJ / REJ |
| 2021 | co2pc__gdppc | 189 | 0.0040 | 0.0070 | 0.0006 | 0.0006 | REJ / REJ / REJ / REJ |
| 2021 | gdp__pop | 199 | 0.8040 | 0.3254 | 0.2623 | 0.2623 | no-rej / no-rej / no-rej / no-rej |
| 2021 | infmort__gdppc | 188 | 0.0090 | 0.0037 | 0.0019 | 0.0019 | REJ / REJ / REJ / REJ |
| 2021 | lifeexp__gdppc | 199 | 0.9280 | 0.5718 | 0.7016 | 0.7014 | no-rej / no-rej / no-rej / no-rej |
| 2021 | lifeexp__gdppc_health | 186 | 0.5410 | 0.4088 | 0.4770 | 0.5452 | no-rej / no-rej / no-rej / no-rej |
| 2021 | netusers__gdppc | 176 | 0.0010 | 0.0000 | 0.0000 | 0.0000 | REJ / REJ / REJ / REJ |
| 2022 | co2pc__gdppc | 189 | 0.0040 | 0.0070 | 0.0003 | 0.0003 | REJ / REJ / REJ / REJ |
| 2022 | gdp__pop | 199 | 0.8610 | 0.2774 | 0.2195 | 0.2195 | no-rej / no-rej / no-rej / no-rej |
| 2022 | infmort__gdppc | 188 | 0.0430 | 0.0569 | 0.1594 | 0.1594 | REJ / no-rej / no-rej / no-rej |
| 2022 | lifeexp__gdppc | 199 | 0.0080 | 0.0002 | 0.5034 | 0.5038 | REJ / REJ / no-rej / no-rej |
| 2022 | lifeexp__gdppc_health | 185 | 0.1080 | 0.0002 | 0.5020 | 0.8230 | no-rej / REJ / no-rej / no-rej |
| 2022 | netusers__gdppc | 178 | 0.0010 | 0.0000 | 0.0000 | 0.0000 | REJ / REJ / REJ / REJ |
| 2023 | co2pc__gdppc | 187 | 0.0010 | 0.0003 | 0.0003 | 0.0003 | REJ / REJ / REJ / REJ |
| 2023 | gdp__pop | 197 | 0.9120 | 0.2548 | 0.1891 | 0.1891 | no-rej / no-rej / no-rej / no-rej |
| 2023 | infmort__gdppc | 187 | 0.0030 | 0.0025 | 0.0002 | 0.0002 | REJ / REJ / REJ / REJ |
| 2023 | lifeexp__gdppc | 197 | 0.9400 | 0.3584 | 0.2089 | 0.2084 | no-rej / no-rej / no-rej / no-rej |
| 2023 | lifeexp__gdppc_health | 184 | 0.9080 | 0.6931 | 0.6278 | 0.6822 | no-rej / no-rej / no-rej / no-rej |
| 2023 | netusers__gdppc | 177 | 0.0010 | 0.0000 | 0.0000 | 0.0000 | REJ / REJ / REJ / REJ |

## 备注
- 真实数据无已知真值, 本表只呈现判定差异, 不含功效/校准结论。
- p_fc 为排列经验 p(5 种子中位数), RESET 三个 p 值为确定性解析值。
