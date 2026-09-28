"""
dimreflex
=========

A dimensional-aware engineering model discovery framework.
量纲感知的工程模型发现框架。

核心理念:
- 数据进入框架前必须携带单位、量纲、变量角色等元信息
- 函数发现在无量纲空间进行,而非原始变量空间
- 数据不足时框架显式降级,绝不硬输出
- 仅服务工程认知,不接管生产端控制
"""

from .anomaly import (
    check_data_health,
    detect_collinearity,
    detect_leverage,
    detect_residual_outliers,
    fit_with_guard,
    ridge_fit,
    robust_irls,
)
from .autodiscovery import assess_family, discover_auto, predict_candidate
from .core import System, Variable
from .discover import (
    compute_groups,
    discover,
    discover_auto_pi,
    discover_multi_pi,
    fit_power_law,
    form_check,
)
from .families import (
    FAMILIES,
    CandidateModel,
    FunctionFamily,
    PiecewisePowerLawFamily,
    PowerLawFamily,
    RationalLogFamily,
    SymbolicRegressionFamily,
    discover_multi,
    register_family,
)
from .monitoring import (
    EngineeringMonitor,
    population_stability_index,
)
from .scenarios import fit_scenarios, predict_scenarios
from .templates import Group, fit_with_templates
from .uncertainty import propagate_uncertainty

__version__ = "0.5.0"
__all__ = [
    'FAMILIES',
    'CandidateModel',
    'EngineeringMonitor',
    'FunctionFamily',
    'Group',
    'PiecewisePowerLawFamily',
    'PowerLawFamily',
    'RationalLogFamily',
    'SymbolicRegressionFamily',
    'System',
    'Variable',
    'assess_family',
    'check_data_health',
    'compute_groups',
    'detect_collinearity',
    'detect_leverage',
    'detect_residual_outliers',
    'discover',
    'discover_auto',
    'discover_auto_pi',
    'discover_multi',
    'discover_multi_pi',
    'fit_power_law',
    'fit_scenarios',
    'fit_with_guard',
    'fit_with_templates',
    'form_check',
    'population_stability_index',
    'predict_candidate',
    'predict_scenarios',
    'propagate_uncertainty',
    'register_family',
    'ridge_fit',
    'robust_irls',
]
