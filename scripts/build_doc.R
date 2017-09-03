library(ezknitr)
ezspin(file = "scripts/model_fn.R", out_dir = "doc", keep_md = FALSE,)
ezspin(file = "scripts/model_run.R", out_dir = "doc", keep_md = FALSE, chunk_opts = list(eval=FALSE))

