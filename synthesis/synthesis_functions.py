import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass

from scipy import stats
from scipy.stats import genextreme as gev
from scipy.stats import genextreme
from scipy.optimize import minimize_scalar


#defining the functions used for synthesis plotting, showcased in synthesis_example.ipynb

#calculating linear regression between GMST and the extreme event, used both in models and obs
def linreg(x,y):
    a, b = np.polyfit(x, y, 1)
    y_pred = a * x + b
    return a, b, y_pred


# based on the attribution_serbia.ipynb code piece, used for callculating the attribution estimate and the lower and upper bound of the confidence intervals
# modified to provide output in the form usable by other functions

def get_confidence_levels(extreme_data, mu, mu0, dist_name='genextreme', ci=95 ):
    """
    Getting the shifted observational data for the mu value of choice and finding the lower and upper confidence intervals of the observed event of choice
    via bootstraping.

    Parameters:
    - extreme_data: observed extreme event values
    - mu0: the observed mean of the distribution
    - mu: the shifted value of mu
    - dist_name: distribution name (default 'genextreme')
    - ci: confidence interval percentage (default 95)
    """
    extreme_data = extreme_data + mu - mu0
    bootstrapped_levels = np.array([])

    for _ in range(1000):
        resampled_data = np.random.choice(extreme_data, size = len(extreme_data), replace=True)

         #modified to return the 5 and 95% confidence intervals only based on the data, not the conversion to RP
        bootstrapped_levels = np.append(bootstrapped_levels,resampled_data)
    
    lower_bound = np.percentile(bootstrapped_levels, (100 - ci) / 2, axis=0)
    upper_bound = np.percentile(bootstrapped_levels, 100 - (100 - ci) / 2, axis=0)
    extreme_data1 = extreme_data
    return lower_bound, upper_bound

## NOT YET IMPLEMENTED
def get_confidence_levels_abs(event, extreme_data, mu1, mu2, mu0, dist_name='genextreme', ci=95 ):
    """
    Getting the shifted observational data for the mu value of choice and finding the lower and upper confidence intervals of the observed event of choice
    via bootstraping, expressed in absolute values.

    Parameters:
    - event: observed value of the studied extreme event
    - extreme_data: observed extreme event values
    - mu0: the observed mean of the distribution
    - mu1,2: the shifted value of mu, with indices 1 and 2 corresponding to different climates in the observed event
    - dist_name: distribution name (default 'genextreme')
    - ci: confidence interval percentage (default 95)
    NOT YET IMPLEMENTED
    """

    return "NOT YET IMPLEMENTED


## NOT YET IMPLEMENTED
def get_confidence_levels_rc(event, extreme_data, mu1, mu2, mu0, dist_name='genextreme', ci=95 ):
    """
    Getting the shifted observational data for the mu value of choice and finding the lower and upper confidence intervals of the observed event of choice
    via bootstraping, expressed in relative change.

    Parameters:
    - event: observed value of the studied extreme event
    - extreme_data: observed extreme event values
    - mu0: the observed mean of the distribution
    - mu1,2: the shifted value of mu, with indices 1 and 2 corresponding to different climates in the observed event
    - dist_name: distribution name (default 'genextreme')
    - ci: confidence interval percentage (default 95)
    NOT YET IMPLEMENTED
    """

    return "NOT YET IMPLEMENTED"


def get_confidence_levels_pr(event, extreme_data, mu1, mu2, mu0, dist_name='genextreme', ci=95 ):
    """
    Getting the shifted observational data for the mu value of choice and finding the lower and upper confidence intervals of the observed event of choice
    via bootstraping, expressed as probability ratios.

    Parameters:
    - event: observed value of the studied extreme event
    - extreme_data: observed extreme event values
    - mu0: the observed mean of the distribution
    - mu1,2: the shifted value of mu, with indices 1 and 2 corresponding to different climates in the observed event
    - dist_name: distribution name (default 'genextreme')
    - ci: confidence interval percentage (default 95)
    Uses the function below to calculate the return periods of both climates and then returns their probability ratio. 
    The function is taken from the original attribution code, where it was used similarily.
    """
    def confidence_level(event, extreme_data, mu, mu0, dist_name='genextreme', ci=95 ):
        extreme_data = extreme_data + mu - mu0
        dist = getattr(stats, dist_name)
        params = dist.fit(extreme_data)
    
        #shifting the observational data to the desired year
        
        x_range_max = dist.ppf(0.9999, *params)
        x_range = np.linspace(min(extreme_data), x_range_max, 100) #max(extreme_data)+5, 100)
        #Return periods for all observed data
        sorted_data = (np.sort(extreme_data))
        n = len(sorted_data)
        emp_cdf_vals = np.arange(1, n + 1) / (n + 1)
        emp_return_periods = 1 / (1 - emp_cdf_vals)
    
        #re-fitting of the shifted data to the selected distribution
        cdf_vals = dist.cdf(x_range, *params)
        return_periods = 1 / (1 - cdf_vals)
    
        #finding the cdf of the event in the shifted distribution
        cdf_event = dist.cdf(event, *params)
        
        #Bootstrap confidence intervals
        bootstrapped_levels = np.array([])
        boot_cdf = []
            
        for _ in range(1000):
            resampled_data = np.random.choice(extreme_data, size = len(extreme_data), replace=True)
            boot_params = dist.fit(resampled_data)
            boot_cdf_vals = dist.cdf(x_range, *boot_params)
            boot_cdf_vals = np.clip(boot_cdf_vals, 0, 0.99999)
            boot_return_periods = 1 / (1 - boot_cdf_vals)
            
             #modified to return the 5 and 95% confidence intervals only based on the data, not the conversion to RP        
            if not np.inf in boot_return_periods:
                bootstrapped_levels = np.append(bootstrapped_levels, boot_return_periods)
            boot_cdf.append(boot_cdf_vals)
    
        # Selecting the bootstrapped values (cdf & return periods) corrseponding to the selected confidence levels
        lower_boot_cdf =  np.percentile(boot_cdf, (100 - ci) / 2, axis=0)
        upper_boot_cdf = np.percentile(boot_cdf, 100 - (100 - ci) / 2, axis=0)
        
        index = np.where(cdf_vals <= cdf_event)[0]
        lower = lower_boot_cdf[index[-1]]
        upper = upper_boot_cdf[index[-1]]
    
        index = np.where(cdf_vals >= cdf_event)[0]
        lower1 = lower_boot_cdf[index[0]]
        upper1 = upper_boot_cdf[index[0]]
        return cdf_event, (lower+lower1)/2, (upper+upper1)/2
    now, lower_now, upper_now = confidence_level(event, extreme_data, mu1, mu0 )
    past, lower_old, upper_old = confidence_level(event, extreme_data, mu2, mu0 )

    return (1-now)/(1-past), (1-lower_now)/(1-lower_old), (1-upper_now)/(1-upper_old)

# the core of the synthesis process, based on the R code by WWA (https://github.com/WorldWeatherAttribution/rwwa/blob/main/R/synthesis.R)

def getsynmean(data: pd.DataFrame, sig_mod: float = 0.0) -> pd.Series:
    """
    Weighted mean of climate-model results with optional model representation error (sig_mod).
    data must have columns: 'est', 'lower', 'upper' with a single value, (?next part)corresponding to the mean of the variable 
    Returns: pd.Series(index=['est','lower','upper'])
    Should be used after sig_mod is evaluated with the help of chi^2 (this can be skipped)
    """
    # inverse-variance weights from 95% CI width + model representation error
    # (upper-lower) = 2 * 1.96 * sigma  => sigma = (upper-lower) / (2*1.96)
    sigma = (data['upper'] - data['lower']) / (2 * 1.96)
    w = 1.0 / (sigma*sigma + sig_mod*sig_mod) #respective weights for each model
    w1 = w.sum() # the cumulative weight for all models

    # weighted mean of point estimates
    s1 = (w * data['est']).sum() / w1

    # “weighted interval” by adding variances in quadrature on each side
    sig_lower = np.sqrt((w * (((data['est'] - data['lower']) / 1.96) ** 2 + sig_mod*sig_mod)).sum() / w1)
    sig_upper = np.sqrt((w * (((data['est'] - data['upper']) / 1.96) ** 2 + sig_mod*sig_mod)).sum() / w1)

    out = pd.Series(
        [s1, s1 - 1.96 * sig_lower, s1 + 1.96 * sig_upper],
        index=['est', 'lower', 'upper']
    )
    return out


def getsynchi2(data: pd.DataFrame, sig_mod: float = 0.0) -> float:
    """
    Chi^2 used to estimate sig_mod (find sig_mod so that chi^2 / dof ~ 1).
    """
    s1 = getsynmean(data, sig_mod)['est']

    # for each model: use lower- or upper-side sigma depending on the sign of (est - s1)
    # lower-side sigma = (est - lower)/1.96 ; upper-side sigma = (est - upper)/1.96 (note: upper < est)
    over = data['est'] > s1 # True when in data['est] > s1, else False
    sig_side = np.where(
        over, #condition
        ((data['est'] - data['lower']) / 1.96) ** 2 + sig_mod**2, # if the condition applies
        ((data['est'] - data['upper']) / 1.96) ** 2 + sig_mod**2 # else (if the condition doesn't apply)
    )
    chi2 = np.sum(((data['est'] - s1) ** 2) / sig_side)
    return float(chi2)

# ------------------------------------------------------
# Main synthesis (observations + models)
# ------------------------------------------------------
@dataclass
class SynthesisResult:
    synth_type: str
    sig_obs: float | None
    chi2_over_dof: float
    sig_mod: float
    df: pd.DataFrame
    uw_mean: float
#ABS AND REL NOT YET IMPLEMENTED
def synthesis(obs_in: pd.DataFrame | None,
              models_in: pd.DataFrame,
              synth_type: str = "abs") -> SynthesisResult:
    """
    Combine observational and model attribution estimates.
    - obs_in/models_in: DataFrames with columns ['est','lower','upper'] and index = dataset/model names.
      If obs_in is None or empty, a dummy obs is used and the final table includes only model lines.
    - synth_type: 'abs' (absolute), 'rel' (percent changes), 'PR' (probability ratios)

    Returns SynthesisResult with:
      - sig_obs: observational representation error (None if no obs provided)
      - sig_mod: model representation error (optimized so chi2/dof ~ 1 when needed)
      - chi2_over_dof: initial ratio with sig_mod=0
      - df: long-form table with individual items and synthesized rows
      - uw_mean: unweighted mean of obs & models center (after inverse-transform)
    """

    no_obs = False
    if obs_in is None or len(obs_in) == 0:
        # dummy single-row obs to preserve code path; dropped later
        no_obs = True
        obs_in = pd.DataFrame({'est': [0.0], 'lower': [0.0], 'upper': [0.0]}, index=['dummy'])

    # standardize column names
    obs_in = obs_in[['est', 'lower', 'upper']].copy()
    models_in = models_in[['est', 'lower', 'upper']].copy()

    # ensure we have a 'model' column (like R code uses rownames)
    if 'model' not in obs_in.columns:
        obs_in = obs_in.copy()
        obs_in['model'] = obs_in.index.astype(str)
    if 'model' not in models_in.columns:
        models_in = models_in.copy()
        models_in['model'] = models_in.index.astype(str)

    # transform to log-space if PR, or log(1+x/100) for relative (%)
    def _forward(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        if synth_type == "PR":
            df[['est', 'lower', 'upper']] = np.log(df[['est', 'lower', 'upper']])
        elif synth_type == "rel":
            df[['est', 'lower', 'upper']] = np.log(1.0 + df[['est', 'lower', 'upper']] / 100.0)
        return df

    def _backward(values: pd.DataFrame | pd.Series | float) -> pd.DataFrame | pd.Series | float:
        if synth_type == "PR":
            return np.exp(values)
        elif synth_type == "rel":
            return 100.0 * (np.exp(values) - 1.0)
        else:
            return values

    obs_in_t = _forward(obs_in)
    models_in_t = _forward(models_in)

    # ---------- Observational representation error (scatter across datasets)
    nobs = len(obs_in_t)
    obs_mean = obs_in_t[['est', 'lower', 'upper']].mean(axis=0).to_numpy()
    if nobs == 1:  ##if we have a single observational dataset, which holds for our case
        sig_obs = 0.0
    else:
        s2 = np.sum((obs_in_t['est'].to_numpy() - obs_mean[0]) ** 2)
        sig_obs = np.sqrt(s2 / (nobs - 1))

    # widen each obs interval by adding (1.96*sig_obs)^2 in quadrature
    obs_in_t['l_wb'] = obs_in_t['est'] - np.sqrt((obs_in_t['est'] - obs_in_t['lower']) ** 2 + (1.96 * sig_obs) ** 2)
    obs_in_t['u_wb'] = obs_in_t['est'] + np.sqrt((obs_in_t['est'] - obs_in_t['upper']) ** 2 + (1.96 * sig_obs) ** 2)

    # also extend the *pooled* obs interval
    obs_c = obs_mean.copy()  # [est, lower, upper] (on transformed scale)
    obs_c[1] = obs_c[0] - np.sqrt((obs_c[0] - obs_c[1]) ** 2 + (1.96 * sig_obs) ** 2)
    obs_c[2] = obs_c[0] + np.sqrt((obs_c[0] - obs_c[2]) ** 2 + (1.96 * sig_obs) ** 2)

    # ---------- Model representation error (optimize sig_mod so chi2/dof ~ 1, if needed)
    chi2_0 = getsynchi2(models_in_t, sig_mod=0.0)
    mdof = max(len(models_in_t) - 1, 1)
    if chi2_0 / mdof > 1.0:
        # bounded Brent search in [0, 5], as in R code https://en.wikipedia.org/wiki/Brent%27s_method
        def _obj(x):
            return (getsynchi2(models_in_t, sig_mod=x) - (len(models_in_t) - 1)) ** 2

        res = minimize_scalar(_obj, bounds=(0.0, 5.0), method='bounded') #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalar.html 
        #python implementation of the Brent search algorithm. since method is 'bounded' though, it is not used here
        # _obj is the single variable scalar function to be minimized, passed to the minimize_scalar function as defined in the documentation 
        sig_mod = float(res.x)
    else:
        sig_mod = 0.0

    models_c = getsynmean(models_in_t[['est', 'lower', 'upper']], sig_mod=sig_mod)

    # add widened bands to individual models
    models_in_t['l_wb'] = models_in_t['est'] - np.sqrt((models_in_t['est'] - models_in_t['lower']) ** 2 + (1.96 * sig_mod) ** 2)
    models_in_t['u_wb'] = models_in_t['est'] + np.sqrt((models_in_t['est'] - models_in_t['upper']) ** 2 + (1.96 * sig_mod) ** 2)

    # ---------- Combine obs + models (weighted by inverse of interval width)
    # weights from interval width (upper-lower) on transformed scale
    w_obs = 1.0 / ((obs_c[2] - obs_c[1]) ** 2)
    w_mod = 1.0 / ((models_c['upper'] - models_c['lower']) ** 2)

    wmean = (w_obs * obs_c[0] + w_mod * models_c['est']) / (w_obs + w_mod)

    sig_lower = np.sqrt((w_obs * ((obs_c[0] - obs_c[1]) / 1.96) ** 2 + w_mod * ((models_c['est'] - models_c['lower']) / 1.96) ** 2) / (w_obs + w_mod))
    sig_upper = np.sqrt((w_obs * ((obs_c[0] - obs_c[2]) / 1.96) ** 2 + w_mod * ((models_c['est'] - models_c['upper']) / 1.96) ** 2) / (w_obs + w_mod))
    synth_c = pd.Series(
        [wmean, wmean - 1.96 * sig_lower, wmean + 1.96 * sig_upper],
        index=['est', 'lower', 'upper']
    )

    # Unweighted center + “white band” (average of side variances)
    umean_c = 0.5 * (obs_c[0] + models_c['est'])
    l_wb_synth = umean_c - np.sqrt(((obs_c[0] - obs_c[1]) ** 2 + (models_c['est'] - models_c['lower']) ** 2) / 2.0)
    u_wb_synth = umean_c + np.sqrt(((obs_c[0] - obs_c[2]) ** 2 + (models_c['est'] - models_c['upper']) ** 2) / 2.0)

    # ---------- Assemble output table (on transformed scale)
    obs_tab = obs_in_t.copy()
    obs_tab['group'] = 'obs'

    obs_row = pd.DataFrame({
        'model': ['Observations'],
        'group': ['obs_synth'],
        'est': [obs_c[0]],
        'lower': [obs_c[1]],
        'upper': [obs_c[2]],
        'l_wb': [obs_c[0] - (obs_c[0] - obs_c[1])],  # not used visually (kept for structure)
        'u_wb': [obs_c[0] + (obs_c[2] - obs_c[0])],
    })

    mod_tab = models_in_t.copy()
    mod_tab['group'] = 'models'

    mod_row = pd.DataFrame({
        'model': ['Models'],
        'group': ['model_synth'],
        'est': [models_c['est']],
        'lower': [models_c['lower']],
        'upper': [models_c['upper']],
        'l_wb': [models_c['est'] - (models_c['est'] - models_c['lower'])],  # structure
        'u_wb': [models_c['est'] + (models_c['upper'] - models_c['est'])],
    })

    synth_row = pd.DataFrame({
        'model': ['Synthesis'],
        'group': ['synth'],
        'est': [synth_c['est']],
        'lower': [synth_c['lower']],
        'upper': [synth_c['upper']],
        'l_wb': [l_wb_synth],
        'u_wb': [u_wb_synth],
    })

    # make sure all needed columns exist before concat
    cols = ['group', 'model', 'est', 'lower', 'upper', 'l_wb', 'u_wb']
    for df_ in (obs_tab, obs_row, mod_tab, mod_row, synth_row):
        for c in cols:
            if c not in df_.columns:
                df_[c] = np.nan
        df_[:]  # no-op, just being explicit

    res = pd.concat(
        [obs_tab[cols], obs_row[cols], mod_tab[cols], mod_row[cols], synth_row[cols]],
        axis=0,
        ignore_index=True
    )

    # drop dummy obs if we started with no observations
    if no_obs:
        res = res[res['group'].str.contains('model', na=False)].reset_index(drop=True)
        sig_obs_out = None
    else:
        sig_obs_out = sig_obs

    # ---------- Back-transform to user scale
    for col in ['est', 'lower', 'upper', 'l_wb', 'u_wb']:
        res[col] = _backward(res[col])

    # also back-transform sigmas and unweighted mean (to mirror R code)
    def _maybe_back(x):
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return x
        return float(_backward(x))

    sig_obs_bt = _maybe_back(sig_obs_out)
    sig_mod_bt = float(_backward(sig_mod))
    uw_mean_bt = float(_backward(umean_c))

    return SynthesisResult(
        synth_type=synth_type,
        sig_obs=sig_obs_bt,
        chi2_over_dof=float(chi2_0 / mdof),
        sig_mod=sig_mod_bt,
        df=res.reset_index(drop=True),
        uw_mean=uw_mean_bt
    )

    
def plot_synthesis(synth: SynthesisResult | pd.DataFrame,
                   xlim: tuple[float, float] | None = None,
                   lwd: float = 10.0,
                   xlab: str = "Probability ratio",
                   season: str = "",
                   add_space: bool = True,
                   log: bool | None = None,
                   hide_labels: bool = False):
    """
    Recreates the WWA bar visualization with outer white band, colored CI, and point estimate.

    - If synth is SynthesisResult, uses synth.df; otherwise expects a dataframe with columns:
      ['group','model','est','lower','upper','l_wb','u_wb'].
    - If log is None, auto: PR => log x-axis, else linear.
    """
    main = "Synthesis plot "+season #change compared to the original, added season as a variable, main moved down
    if isinstance(synth, SynthesisResult):
        df = synth.df.copy()
        if log is None:
            logaxs = 'x' if synth.synth_type == 'PR' else ''
        else:
            logaxs = 'x' if log else ''
    else:
        df = synth.copy()
        logaxs = '' if (log is None) else ('x' if log else '')

    # color map (semi-transparent “groups”)
    gcols_map = {
        'obs': (0, 0, 1, 0.5),         # blue alpha
        'obs_synth': (0, 0, 1, 1.0),   # blue
        'models': (1, 0, 0, 0.5),      # red alpha
        'model_synth': (1, 0, 0, 1.0), # red
        'synth': (1, 0, 1, 1.0)        # magenta
    }
    df = df.copy()
    df['group'] = df['group'].astype(str)

    # x-limits
    if xlim is None:
        vals = pd.to_numeric(df[['lower', 'upper', 'l_wb', 'u_wb']].to_numpy().ravel(), errors='coerce')
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            xlim = (0.0, 1.0)
        else:
            lo, hi = np.min(vals), np.max(vals)
            pad = 0.05 * (hi - lo if hi > lo else (abs(hi) + 1))
            xlim = (lo - pad, hi + pad)

    # y positions (replicates the “spacer line” logic)
    nobs = int((df['group'] == 'obs').sum())
    nmod = int((df['group'] == 'models').sum())
    if add_space and nobs > 0:
        # Build like the R code: first obs, then models, then synth, with gaps
        yy = np.arange(len(df), 0, -1)
    else:
        yy = np.arange(len(df), 0, -1)

    # plot
    fig, ax = plt.subplots()
    if logaxs == 'x':
        ax.set_xscale('log')
        vline = 1.0
    else:
        vline = 0.0

    ax.set_xlim(*xlim)
    ax.set_ylim(0.5, len(df) + 0.5)
    ax.set_xlabel(xlab)
    ax.set_title(main)
    ax.grid(axis='x', color=(0, 0, 0, 0.1))
    ax.axvline(vline, linestyle='--', color='k', linewidth=1)

    # draw bands and points
    # map per-row color from group
    colors = [gcols_map.get(g, (0, 0, 0, 1)) for g in df['group']]

    # outer white band (black border effect in R done by two segments; here: white band atop background)
    # We'll emulate by drawing a thick black, then slightly thinner white, then the colored band.
    for i, row in df.reset_index(drop=True).iterrows():
        y = yy[i]
        lwb, uwb = row['l_wb'], row['u_wb']
        low, up = row['lower'], row['upper']
        c = colors[i]

        # safety: skip rows with NaNs
        if not np.all(np.isfinite([lwb, uwb, low, up, row['est']])):
            continue

        # black base band
        ax.plot([lwb, uwb], [y, y], solid_capstyle='butt', linewidth=lwd, color='k')
        # white inner band
        ax.plot([lwb, uwb], [y, y], solid_capstyle='butt', linewidth=max(lwd - 2, 1), color='w')
        # colored CI
        ax.plot([low, up], [y, y], solid_capstyle='butt', linewidth=lwd, color=c)
        # point estimate
        ax.plot(row['est'], y, marker='o', markersize=max(lwd, 8) / 2.5, markeredgewidth=2,
                markerfacecolor=c, markeredgecolor='k')

    # y labels
    if not hide_labels:
        ax.set_yticks(yy)
        ax.set_yticklabels(df['model'])
    else:
        ax.set_yticks([])
    return fig, ax

