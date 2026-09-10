function S = ccf_dorsal_compose(S_a, S2)
%CCF_DORSAL_COMPOSE  Combine FOV->1p and 1p->CCF into one FOV->CCF transform.

T1 = S_a.tform.T;
T2 = S2.tform.T;
S = struct();
S.name = "fov_to_ccf_via_1p";
S.moving_points = S_a.moving_points;
S.fixed_points = S2.fixed_points;
S.transform_type = "composed";
if exist("rigidtform2d", "file") == 2 || exist("affinetform2d", "file") == 2
    try
        S.tform = affinetform2d(T1 * T2);
    catch
        S.tform = affine2d(T1 * T2);
    end
else
    S.tform = affine2d(T1 * T2);
end
S.step_fov_to_1p = S_a;
S.step_1p_to_ccf = S2;
end
