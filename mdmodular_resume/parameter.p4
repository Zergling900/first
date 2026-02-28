# ===== Continue / Resume Control =====
# 0 = start a new run from First_Data_file and reinitialize velocities
# 1 = continue from Data_file (+ Et0 file for s/ps/target T)
continue_run = 0

# 0 = use the last frame/row in files
# 1 = try to resume from a specified target time (fs)
resume_use_target_time = 0
resume_time_fs = -1.0

# When resume_use_target_time = 1, Data and Et0 must match target time within 0.01*dt.
# (This tolerates floating-point roundoff while remaining effectively exact for MD output grids.)
# resume_require_exact_time is kept for future stricter modes; current implementation already enforces the 0.01*dt tolerance.
resume_require_exact_time = 0

# 1 = when continue_run=1, output to next index (recommended, avoids overwrite)
# 0 = reuse current Data_file / Et_file paths (will overwrite)
continue_write_next_index = 1

# Force model dispatch (selected at runtime in Force_Current)
# Current supported: BeW_ABOP_LCL
force_model = BeW_ABOP_LCL
