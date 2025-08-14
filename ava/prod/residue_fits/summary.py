# %%
import numpy as np
from uncertainties import ufloat
# %%

# values of first 3 foldons
first_vals = [0.00351, 0.00326, 0.00386]
first_errors = [0.00129, 0.00095, 0.00099]

ufs = [ufloat(v, e) for v, e in zip(first_vals, first_errors)]
np.mean(ufs)
# 0.0035433333333333337+/-
# 0.0006277561451533372

# %%
second_vals = [0.00204, 0.00184]
second_errors = [0.00033, 0.00040]

ufs = [ufloat(v, e) for v, e in zip(second_vals, second_errors)]
np.mean(ufs)
# 0.00194+/-
# 0.0002592778432492834
# %%
