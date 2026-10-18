# 门卫对照: form_check vs 二次族模型比较(CV/AICc)

- 网格: n∈[50, 200], rel=0.03, 200 种子; 竞争族=log 全二次展开; CV=5 折; 校准=H0 池逐 bench×n 调阈值至 FPR≤5%

## srsd n=50
- H0(n=2600): fc 0.051 | cv 朴素δ0 0.167 → 校准 0.050 | aicc 朴素m0 0.086 → 校准 0.050
- 结构H1 校准TPR(均值/任务): fc 0.348 | cv 0.462 | aicc 0.525 | OR(fc,cv) 0.547  (任务数 19)
- 逐任务: feynman-bonus.15[srsd-H1]: fc0.09/cv0.04/aicc0.04, feynman-bonus.20[srsd-H1]: fc0.83/cv0.99/aicc1.00, feynman-bonus.5[srsd-H1]: fc0.79/cv1.00/aicc1.00, feynman-i.10.7[srsd-H1]: fc0.26/cv0.15/aicc0.21, feynman-i.26.2[srsd-H1]: fc0.35/cv0.55/aicc1.00, feynman-i.27.6[srsd-H1]: fc0.81/cv1.00/aicc1.00, feynman-i.30.5[srsd-H1]: fc0.84/cv0.68/aicc1.00, feynman-i.37.4[srsd-H1]: fc0.60/cv0.41/aicc0.81, feynman-i.39.11[srsd-H1]: fc0.06/cv0.99/aicc1.00, feynman-i.41.16[srsd-H1]: fc0.07/cv0.10/aicc0.09, feynman-i.43.43[srsd-H1]: fc0.06/cv0.99/aicc1.00, feynman-i.48.2[srsd-H1]: fc0.20/cv0.20/aicc0.22, feynman-ii.11.27[srsd-H1]: fc0.08/cv0.03/aicc0.02, feynman-ii.11.28[srsd-H1]: fc0.06/cv0.09/aicc0.04, feynman-ii.13.23[srsd-H1]: fc0.20/cv0.20/aicc0.22, feynman-ii.13.34[srsd-H1]: fc0.20/cv0.20/aicc0.22, feynman-ii.38.14[srsd-H1]: fc1.00/cv1.00/aicc1.00, feynman-iii.4.32[srsd-H1]: fc0.05/cv0.10/aicc0.05, feynman-iii.4.33[srsd-H1]: fc0.05/cv0.10/aicc0.05

## srsd n=200
- H0(n=2600): fc 0.038 | cv 朴素δ0 0.152 → 校准 0.050 | aicc 朴素m0 0.123 → 校准 0.050
- 结构H1 校准TPR(均值/任务): fc 0.493 | cv 0.628 | aicc 0.669 | OR(fc,cv) 0.675  (任务数 19)
- 逐任务: feynman-bonus.15[srsd-H1]: fc0.06/cv0.06/aicc0.07, feynman-bonus.20[srsd-H1]: fc0.99/cv1.00/aicc1.00, feynman-bonus.5[srsd-H1]: fc1.00/cv1.00/aicc1.00, feynman-i.10.7[srsd-H1]: fc0.74/cv0.61/aicc0.73, feynman-i.26.2[srsd-H1]: fc0.43/cv0.86/aicc1.00, feynman-i.27.6[srsd-H1]: fc1.00/cv1.00/aicc1.00, feynman-i.30.5[srsd-H1]: fc1.00/cv0.86/aicc1.00, feynman-i.37.4[srsd-H1]: fc0.94/cv0.90/aicc1.00, feynman-i.39.11[srsd-H1]: fc0.04/cv1.00/aicc1.00, feynman-i.41.16[srsd-H1]: fc0.03/cv0.08/aicc0.10, feynman-i.43.43[srsd-H1]: fc0.04/cv1.00/aicc1.00, feynman-i.48.2[srsd-H1]: fc0.64/cv0.80/aicc0.86, feynman-ii.11.27[srsd-H1]: fc0.05/cv0.02/aicc0.03, feynman-ii.11.28[srsd-H1]: fc0.04/cv0.06/aicc0.07, feynman-ii.13.23[srsd-H1]: fc0.64/cv0.80/aicc0.86, feynman-ii.13.34[srsd-H1]: fc0.64/cv0.80/aicc0.86, feynman-ii.38.14[srsd-H1]: fc1.00/cv1.00/aicc1.00, feynman-iii.4.32[srsd-H1]: fc0.04/cv0.06/aicc0.07, feynman-iii.4.33[srsd-H1]: fc0.04/cv0.06/aicc0.07

## synth n=50
- H0(n=600): fc 0.080 | cv 朴素δ0 0.212 → 校准 0.050 | aicc 朴素m0 0.132 → 校准 0.050
- 噪声层H1(应安静): fc 0.072 | cv 0.063 | aicc 0.061
- 结构H1 校准TPR(均值/任务): fc 0.774 | cv 0.929 | aicc 0.930 | OR(fc,cv) 0.937  (任务数 11)
- 逐任务: COLLINEAR_CURV[collinear-quad]: fc1.00/cv1.00/aicc1.00, CUBIC_LIN[cubic]: fc0.32/cv0.28/aicc0.31, DILUTION[dilution]: fc1.00/cv1.00/aicc1.00, EXP_MODEL[exponential]: fc1.00/cv1.00/aicc1.00, INTERACT[interaction]: fc0.24/cv1.00/aicc1.00, LOG_MODEL[logarithmic]: fc1.00/cv1.00/aicc1.00, OMITTED_M[omitted-var]: fc0.93/cv0.94/aicc0.92, OPPOSING[opposing]: fc0.03/cv1.00/aicc1.00, QUAD_LIN[quadratic]: fc1.00/cv1.00/aicc1.00, RATIONAL[rational]: fc1.00/cv1.00/aicc1.00, THRESHOLD[piecewise]: fc1.00/cv1.00/aicc1.00

## synth n=200
- H0(n=600): fc 0.040 | cv 朴素δ0 0.183 → 校准 0.050 | aicc 朴素m0 0.135 → 校准 0.050
- 噪声层H1(应安静): fc 0.046 | cv 0.041 | aicc 0.044
- 结构H1 校准TPR(均值/任务): fc 0.782 | cv 0.934 | aicc 0.940 | OR(fc,cv) 0.941  (任务数 11)
- 逐任务: COLLINEAR_CURV[collinear-quad]: fc1.00/cv1.00/aicc1.00, CUBIC_LIN[cubic]: fc0.33/cv0.28/aicc0.34, DILUTION[dilution]: fc1.00/cv1.00/aicc1.00, EXP_MODEL[exponential]: fc1.00/cv1.00/aicc1.00, INTERACT[interaction]: fc0.26/cv1.00/aicc1.00, LOG_MODEL[logarithmic]: fc1.00/cv1.00/aicc1.00, OMITTED_M[omitted-var]: fc1.00/cv1.00/aicc1.00, OPPOSING[opposing]: fc0.02/cv1.00/aicc1.00, QUAD_LIN[quadratic]: fc1.00/cv1.00/aicc1.00, RATIONAL[rational]: fc1.00/cv1.00/aicc1.00, THRESHOLD[piecewise]: fc1.00/cv1.00/aicc1.00

## 成本(拟合次数): fc=2 全量拟合 + 1 pinv(999 排列向量化); cv_quad=12 次拟合(2 模型×6); aicc_quad=2 次拟合
