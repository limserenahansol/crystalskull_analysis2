function D = extract_cohort_delta_from_session_pair(bl_path, inj_path, params)
%EXTRACT_COHORT_DELTA_FROM_SESSION_PAIR  Δ metrics for one mouse (inj − baseline).
%   Same preprocessing as extract_analysis2_core via extract_session_metrics.
%
%   params.match_inj_last_to_baseline_dur — if true, injection window is the last
%   m_bl.dur_s seconds (after baseline crops), same rule as extract_analysis2_core.

m_bl = extract_session_metrics(bl_path, params, 'baseline');
p_inj = params;
if cohort_getp(p_inj, 'match_inj_last_to_baseline_dur', false)
    p_inj.use_last_sec_inj = m_bl.dur_s;
    p_inj.use_first_sec_inj = 0;
end
m_inj = extract_session_metrics(inj_path, p_inj, 'injection');
D = compute_delta_fields(m_inj, m_bl);

end

function D = compute_delta_fields(m_inj, m_bl)
fn = {'mean_mean_dff', 'median_mean_dff', 'mean_std_dff', 'median_std_dff', ...
    'mean_skew_dff', 'median_skew_dff', 'mean_event_rate', 'median_event_rate', ...
    'mean_pop_frac', 'mean_pop_rate_z', 'median_pw_corr', 'pc1_var_pct', 'pc1_to_pc3_var_pct'};
D = struct();
for k = 1:numel(fn)
    f = fn{k};
    D.(['d_' f]) = m_inj.(f) - m_bl.(f);
end
D.d_n_neurons = m_inj.n_neurons - m_bl.n_neurons;
D.bl_n_neurons = m_bl.n_neurons;
D.inj_n_neurons = m_inj.n_neurons;
end

function v = cohort_getp(s, name, default)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
