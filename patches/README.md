# patches/0001-form-check-battery — 形式合格性电池(未入主线,补丁形态)

## 这是什么
`discover.py` 的 core 升级:**form_check 报警时回答"那什么才合适"**。
四族电池(power/exponential/logarithmic/rational),每族在各自的
线性化坐标 (u(x), v(y)) 下跑同构的二次曲率检验,wild bootstrap
(Rademacher) 零分布,viable = 未被拒绝的族(候选列表语义)。

**真实 core 零改动**(2026-09-04 用户权限指示:core 变更跳过)。
本补丁已在临时 worktree 全量验证,采用即 `git apply` 一次。

## 验证状态(全部实测,无推测)
- `git apply --check` PASS(对 f9156e1 基线,应用前再 check 一次,
  若 families_analyze.py 等并行改动先落,不影响本补丁——只碰 discover.py)
- worktree 副本上 **153 测试全绿**(144 旧 + 9 新 test_form_check_battery.py)
- PK 对数律数据上:报警文本升级为"幂律被拒,可行族=[...]"/"无可行族→升级"
- 确定性:固定种子 viable 逐位可复现

## 应用
```
cd C:\Users\16910\Desktop\dimreflex
git apply experiments/form_check_vs_reset/patches/0001-form-check-battery.patch
# 测试文件单独拷贝(不在 diff 里):
cp experiments/form_check_vs_reset/patches/test_form_check_battery.py tests/
"$LOCALAPPDATA/Programs/Python/Python313/python.exe" -m pytest tests/ -q
```
回滚:`git checkout discover.py && rm tests/test_form_check_battery.py`

## 实现要点(写代码时的三个决定)
1. **wild 不是排列**:v(y)=y / v(y)=1/y 坐标下乘性噪声异方差,排列
   可交换性破裂(rational 假升级 30%→wild 4.5-9%)。零分布合法性是
   坐标依赖性质——battery_probe vs battery_probe_wild 对比实验在
   runs/battery_probe*/。
2. **幂律门不自我背书**:form_check 报警时 power 从 viable 移除
   (同一统计量同源);不报警时电池照常跑,留痕不报警(不破坏干净输出)。
3. **每族各自 α=0.05,不做族间多重校正**:viable 是"各自不拒绝"的
   可行集,非拒绝语义下多重性不膨胀第一类错误;假排除/假包含分开报告。

## 已知边界(docstring 与论文节均已声明)
- 小样本尺寸膨胀:n≲60 时 wild size≈6.5%,n=120 收敛 5.0%(200 种子实测)
- 假包含已知来源:H0_MULTIVAR 上 rational 63% 假包含(均值平方列只捕
  对角曲率);OPPOSING/INTERACT 携带 fc 同源盲区
- 单特征版;多特征(逐列 u、均值平方列)代码已支持但未正式验证
- exponential 盲区(时间常数≫采样量程)是域内信息极限,电池不解决

## 测试文件
test_form_check_battery.py(9 测试)在 patches/ 内,**没有**打进
patch 文件(git diff 只覆盖 discover.py)——应用时手动拷贝。
其中 2 个测试用多种子中位断言而非单种子(单固定种子会抽进 wild
尺寸尾巴,n=60 实测 6.5%,那是校准行为不是缺陷)。
