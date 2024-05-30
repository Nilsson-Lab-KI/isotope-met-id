#
# Setup directories
#
# This is needed since git does not track empty directories
#

source("common.R")


create_dir_if_not_exists <- function(dir_path)
{
    if(!file.exists(dir_path))
        dir.create(file.path(mainDir, subDir))
}

create_dir_if_not_exists(preprocessed_data_path)
create_dir_if_not_exists(mi_data_path)
