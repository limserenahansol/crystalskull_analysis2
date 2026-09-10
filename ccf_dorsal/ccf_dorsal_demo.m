function ccf_dorsal_demo()
%CCF_DORSAL_DEMO  How to register a crystal-skull FOV to dorsal CCF.
%
% Interactive: opens cpselect. Close the landmark-guide figure, click
% matching points, then File > Export Points / close cpselect to continue.

thisDir = fileparts(mfilename("fullpath"));
addpath(thisDir);
addpath(fileparts(thisDir));

% --- 1) One-step: mesoscope summary -> dorsal atlas ---
% Use FlipLR=true for mirror2 (same as lab maps: P left, A right).
%
% S = ccf_dorsal_register(extFile, 'FlipLR', true, ...
%     'OutDir', fullfile(fileparts(thisDir), 'ccf_dorsal_output', 'example'));

% --- 2) Better two-step if you have a 1p skull-off / large-FOV photo ---
% S = ccf_dorsal_register(extFile, 'FlipLR', true, ...
%     'Photo1p', 'C:\path\to\skull_off_1p.png', ...
%     'OutDir', fullfile(fileparts(thisDir), 'ccf_dorsal_output', 'example'));

% --- 3) Real Allen top-down PNG (AllenAtlasTopDown screenshot, etc.) ---
% S = ccf_dorsal_register(extFile, 'FlipLR', true, ...
%     'AtlasPng', 'C:\path\to\allen_dorsal_Aup.png');

% Preview the built-in atlas (no clicks):
Ref = ccf_dorsal_make_reference("Orient", "Pleft_Aright");
figure("Color", "w");
imshow(Ref.image);
title("Built-in dorsal sketch (P left, A right). Replace with Allen PNG for area names.");
end
