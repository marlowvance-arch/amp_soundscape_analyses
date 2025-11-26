library(reticulate)

# Point to your project root
setwd("C:/Users/rockn/OneDrive/Desktop/amp_soundscape_analyses")

# Source the Python script
source_python("C:/Users/rockn/OneDrive/Desktop/amp_soundscape_analyses/scripts/indices_py_test.r")

# Run indices with parallel processing
run_acoustic_indices(parallel = TRUE, n_workers = 4)
