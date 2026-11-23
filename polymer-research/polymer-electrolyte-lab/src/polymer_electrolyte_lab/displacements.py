"""Einstein slopes from unwrapped 3D PEO-frame ion displacements with formal transported charges."""
import numpy as np

from .model import BOLTZMANN, ELEMENTARY_CHARGE
from .shell import aggregate, integer, object_record, plot, positive, result, table, text_field


def ion_moments(positions_nm, polymer_com_nm, dt_ps, n_li, lag_steps):
    """Time-origin averages; order is all Li ions then all whole TFSI anion COMs."""
    r = np.asarray(positions_nm, dtype=float)
    com = np.asarray(polymer_com_nm, dtype=float)
    n = integer(n_li, "n_li", 1, 10**6)
    dt = positive(dt_ps, "dt_ps")
    if r.ndim != 3 or r.shape[1:] != (2*n, 3) or com.shape != (len(r),3):
        raise ValueError("Use (frame, 2*n_li, 3) unwrapped positions and (frame, 3) unwrapped polymer COM.")
    if not np.isfinite(r).all() or not np.isfinite(com).all() or len(r)<3:
        raise ValueError("Provide finite, unwrapped ion positions and polymer COM for each frame.")
    r = r - com[:,None,:]
    out = []
    for lag in lag_steps:
        j = integer(lag, "lag steps", 1, len(r)-2)
        dr = r[j:]-r[:-j]
        plus, minus = dr[:,:n], dr[:,n:]
        dq = plus.sum(axis=1)-minus.sum(axis=1)
        def mean_sq(v):
            return float(np.mean(np.sum(v*v,axis=-1)))
        out.append([j*dt, mean_sq(plus), mean_sq(minus), mean_sq(dq)])
    return out


def window_fit(t, y):
    if np.any(y <= 0):
        return None, None, None
    slope = float(np.polyfit(t, y, 1)[0])
    exponent = float(np.polyfit(np.log(t), np.log(y), 1)[0])
    halves = np.array_split(np.arange(len(t)), 2)
    slopes = [float(np.polyfit(t[idx], y[idx], 1)[0]) for idx in halves]
    if slope <= 0:
        return exponent, None, None
    difference = abs(slopes[0]-slopes[1])/slope
    accepted = (0.8 <= exponent <= 1.2 and difference <= 0.3 and min(slopes)>0)
    return exponent, difference, slope if accepted else None


def displacement_result(record):
    rec = object_record(record, "Displacement record")
    if rec.get("frame") != "polymer_com" or rec.get("transport_charge") != "formal_1e":
        raise ValueError("Use unwrapped neutral polymer-COM frame and formal transported charges ±e.")
    state = text_field(rec,"state_id")
    temp = positive(rec.get("temperature_K"),"temperature_K")
    volume = positive(rec.get("volume_nm3"),"volume_nm3")*1e-27
    n = integer(rec.get("n_li"),"n_li",1,10**6)
    window = rec.get("fit_ps")
    if not isinstance(window, list) or len(window) != 2:
        raise ValueError("fit_ps must give two positive lag limits.")
    low, high = [positive(v, "fit lag") for v in window]
    if low >= high:
        raise ValueError("The fit window must have increasing lag limits.")
    runs = rec.get("runs")
    if not isinstance(runs,list) or not 1 <= len(runs) <= 64:
        raise ValueError("Supply 1–64 independent replicas.")
    rows, diagnostics, groups, seen, series, provenance = [], [], {}, set(), [], {}
    for run in runs:
        if not isinstance(run, dict) or run.get("state_id") != state:
            raise ValueError("The displacement input must match its material-state identifier.")
        model, replica = text_field(run,"model"), text_field(run,"replica")
        src = text_field(run,"force_field_source")
        text_field(run,"trajectory_source")
        if provenance.setdefault(model,src) != src or (model,replica) in seen:
            raise ValueError("Use consistent force-field provenance and unique replicate IDs per model.")
        seen.add((model,replica))
        try:
            data = np.asarray(run.get("moments"), dtype=float)
        except (ValueError, TypeError) as exc:
            raise ValueError("Moment rows need lag, Li self, TFSI self and collective moment.") from exc
        if data.ndim != 2 or data.shape[1] != 4 or not 6 <= len(data) <= 10000 or not np.isfinite(data).all():
            raise ValueError("Supply 6–10000 finite four-column moment rows.")
        if np.any(data[:,0]<=0) or np.any(np.diff(data[:,0])<=0) or np.any(data[:,1:]<0):
            raise ValueError("Lags must be positive and strictly increasing; moments must be nonnegative.")
        # Cauchy-Schwarz for |sum signed dr|² <= (2n) sum |dr|², for shared time origins.
        if np.any(data[:,3] > 2*n*n*(data[:,1]+data[:,2])*(1+1e-8)+1e-14):
            raise ValueError("Collective displacement exceeds the neutral ion population's moment bound.")
        subset = data[(data[:,0]>=low)&(data[:,0]<=high)]
        if len(subset)<6:
            raise ValueError("Include at least six lags in the requested fit window.")
        display = data[np.linspace(0, len(data)-1, min(1200, len(data)), dtype=int)]
        coefficients = []
        for j,label in enumerate(("Li self","TFSI self","Collective charge"),1):
            exponent, diff, slope = window_fit(subset[:,0],subset[:,j])
            diagnostics.append([model,replica,label,exponent,diff,"accepted" if slope is not None else "withheld"])
            coefficients.append(slope)
            series.append(dict(name=f"{model} / {replica} / {label}", x=display[:,0].tolist(), y=display[:,j].tolist()))
        # nm²/ps → m²/s = 10^-6. S/m → mS/cm = 10.
        dp, dm, collective = [s*1e-6/6 if s is not None else None for s in coefficients]
        prefactor = 10*ELEMENTARY_CHARGE**2/(volume*BOLTZMANN*temp)
        sigma = prefactor*collective if collective is not None else None
        sigma_ne = prefactor*n*(dp+dm) if dp is not None and dm is not None else None
        rows.append([model,replica,dp,dm,sigma,sigma_ne,sigma/sigma_ne if sigma is not None and sigma_ne else None])
        groups.setdefault(model,[]).append((sigma,sigma_ne))
    summary=[]
    for model, values in groups.items():
        for j, observable in enumerate(("Collective","Nernst–Einstein")):
            accepted=[v[j] for v in values if v[j] is not None]
            mean, sd = aggregate(accepted) if accepted else (None,None)
            summary.append([model,observable,len(accepted),len(values),mean,sd])
    return result([
        table("3D conductivity per model",["Force field","Observable","Accepted runs","Submitted runs","Mean (mS/cm)","Run SD (mS/cm)"],summary),
        table("Transport estimates per run",["Force field","Replica","D Li (m²/s)","D TFSI (m²/s)","Collective (mS/cm)","NE (mS/cm)","Ratio"],rows),
        table("Fixed-window diagnostics",["Force field","Replica","Moment","Log slope","Half-slope difference / full slope","Accepted"],diagnostics)],
        [plot("Displacement moments, including the extensive charge sum", "Lag (ps)", "Moment (nm²)", series)],
        ["The isotropic Einstein charge moment is the squared displacement of the signed sum over all ions. The Li and TFSI moments are per-ion means. The charge basis for the current is physical ±e, independently of electrostatic charge scaling in the dynamics.",
         "A positive 3D diffusive slope is required. The log slope must lie between 0.8 and 1.2, with relative first-half/second-half slope difference ≤0.30. Failure suppresses that estimate. These numerical diagnostics do not establish equilibration, finite-size convergence or a physical plateau.",
         "Use homogeneous constant-volume records from one state at the specified T. Repeat with different fit windows and independently equilibrated trajectories. Report the accepted-run count together with the SD. Collective fluctuations usually require longer sampling."])
