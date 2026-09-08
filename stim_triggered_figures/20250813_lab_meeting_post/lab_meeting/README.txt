Lab-meeting figure pack
========================

per_recording/<mouse_session_power_mirror>/
  fig1 population PSTH for that movie
  fig2 per-neuron STA heatmap
  fig3 excited / inhibited counts
  fig4 poke response on FOV: red = excited, blue = inhibited
       (left anatomy, middle all cells, right p<0.05 only)
  fig5 per-trial raster + ITI

cohort_average/
  Same 5 figure numbers, averaged (or pooled) across all recordings.
  fig4 is a gallery of per-movie response maps (FOVs not aligned).

Spatial footprints were reconstructed from Omer's spatial_weights
(ndSparse CSC inside each *_ext.mat). No re-export was needed.
