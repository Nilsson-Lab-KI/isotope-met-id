#
# Setup directories
#
# This is needed since git does not track empty directories
#

source("common.R")


create_dir_if_not_exists(preprocessed_data_path)
create_dir_if_not_exists(mi_data_path)
create_dir_if_not_exists(mid_distance_path)
