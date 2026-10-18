make_tab <- function(model, directory, title, pred_labels, dv_labels) {
  # Written by Danae Papadopetraki, Jan 2024
  
  library('webshot')
  library('sjPlot')
  
  # Set default value for pred_labels if none provided
  if (missing(pred_labels) || length(pred_labels) == 0) {
    pred_labels <- NULL
  }
  
  # Set default value for dv_labels if none provided
  if (missing(dv_labels) || length(dv_labels) == 0) {
    dv_labels <- NULL 
  }
  
  save_dir <- file.path(directory, 'figures')
  
  # Convert title to lowercase and replace spaces with underscores
  formatted_title <- gsub(" ", "_", tolower(title))
  file_name <- file.path(save_dir, sprintf("%s", formatted_title))
  
  # `model` may be a single model or a list of models. A list produces a
  # side-by-side table (one column per model), e.g. deviance and RT together.
  models <- if (inherits(model, "list")) model else list(model)

  # Build the tab_model() call referencing each model by an index expression
  # (models[[i]]) rather than splicing the objects in via do.call. tab_model
  # deparses its unnamed model arguments to auto-name the columns; a spliced-in
  # model object deparses to hundreds of lines, triggering
  # "'names' attribute [n] must be the same length as the vector [1]".
  # An index expression deparses to a single token, so column naming stays sane.
  model_args <- lapply(seq_along(models), function(i) bquote(models[[.(i)]]))
  # For linear mixed models use Satterthwaite degrees of freedom (matching the emmeans
  # marginal contrasts reported in the text) instead of the Wald normal approximation,
  # so the table p-values agree with the text (e.g. the drug effect p = .043, not .032).
  # generalized (binomial) models use asymptotic Wald z, which is already what glmer uses.
  all_lmer <- all(vapply(models, function(m) inherits(m, "lmerMod") && !inherits(m, "glmerMod"), logical(1)))
  tab_call <- as.call(c(
    quote(sjPlot::tab_model),
    model_args,
    list(
      pred.labels = pred_labels,
      title = title,
      show.aic = T,
      dv.labels = dv_labels,
      p.val = if (all_lmer) "satterthwaite" else "wald",
      p.style = "numeric_stars",
      show.ci = FALSE,
      show.se = T,
      show.stat = T,
      use.viewer = FALSE,  # Disable viewer for better compatibility with scrolling
      CSS = list(
        css.depvarhead = 'font-weight: bold;',
       css.fontsize = 'font-size: 690px;',
        css.table = 'border-collapse: collapse; width: 100%; font-family: Cambria;',
        css.tr = 'margin-bottom: 5px;',
        css.td = 'text-align: center;'  # Center-align the text within table cells
      )
    )
  ))
  table_result <- eval(tab_call)

  # sjPlot writes the interaction separator as a raw Unicode "×" (U+00D7). The
  # saved .html has no charset declaration, so webshot reads the UTF-8 bytes as
  # Latin-1 and renders "Ã—". Swap it for a charset-safe HTML entity so the
  # interaction rows show "setSize × type" correctly. (σ/τ/ρ/R² are already
  # written as entities by sjPlot, hence they render fine.)
  table_result$knitr <- gsub("\u00d7", "&#215;", table_result$knitr, useBytes = FALSE)

  # Rename selected predictor labels wherever they appear, including inside
  # interaction terms (so "setSize \u00d7 type" -> "effort level \u00d7 task-type").
  # Done as word-boundary text replacement on the rendered labels rather than via
  # sjPlot's positional pred.labels, which would need a per-model vector and break
  # whenever a model's term list changes. Names kept lower-case first letter to
  # match the surrounding label convention (setSize, type, age z, ses order, ...).
  # Consistent predictor labels across ALL model tables: capitalize proper terms
  # (Age, UPDRS) and DROP the "z" suffix that scale()/_z pre-scaling leaves behind,
  # so z-scored and non-z versions of a term read identically. Order matters: put
  # each z-suffixed / longer key BEFORE its shorter substring so the substring rule
  # doesn't leave a dangling "z" (perl \b word-boundary matching).
  # Predictor labels: every entry starts with a capital letter so ALL tables read
  # consistently (Pupil, Age, Sex, Effort level, Task-type, ...). Order matters: put
  # each z-suffixed / longer key BEFORE its shorter substring. Values are never
  # re-matched by later keys (checked), so capitalisation is stable.
  label_map <- c(
    # deviation-coded design factors (±0.5) -> marginal-effect coefficient rows.
    # Listed first so they map before the shorter "drug"/"type" keys below.
    "effort dev"     = "Effort level",       # effort_dev: high vs low load
    "type dev"       = "Task-type",          # type_dev: Update vs Ignore
    "drug dev"       = "Drug",               # drug_dev: propranolol vs placebo
    "task dev"       = "Effort task",        # task_dev: no-effort-vs-effort vs low-vs-high
    "group dev"      = "Tremor group",       # group_dev: tremor vs non-tremor
    "pupil"          = "Pupil",              # common renamed pupil column (by-window TD table)
    "wz whole"       = "Pupil",              # scale(wz_whole): within-session z pupil
    "wz enc"         = "Pupil",              # per-phase pupil tables (encoding)
    "wz del"         = "Pupil",              # per-phase pupil tables (delay)
    "wz int"         = "Pupil",              # per-phase pupil tables (interference)
    "raw whole"      = "Pupil residuals",    # scale(raw_whole): gaze-residual pupil
    "effort task"    = "Effort task",        # merged model: nr vs sz decision contrast
    "mean pup int z" = "Mean pupil",        # interference-phase pupil covariate
    "mean pup"       = "Mean pupil",         # scale(mean_pup)
    "updrs total z"  = "UPDRS total",        # z-version first so the "z" is stripped
    "updrs total"    = "UPDRS total",
    "updrs nt z"     = "UPDRS (non-tremor)",
    "updrs nt"       = "UPDRS (non-tremor)",
    "updrs tot"      = "UPDRS total",        # scale(updrs_tot) in the pupil TD models
    "age z"          = "Age",                # scale(Age)/age z -> "Age"
    "Gender"         = "Gender",             # data column Gender -> "Gender"
    "tremor group"   = "Tremor group",       # tremor_group -> "Tremor group"
    "ses order"      = "Session order",      # ses_order -> "Session order"
    "median dev"     = "Performance error",  # group choice covariate / subject-level DV
    "median rt"      = "Reaction time",      # subject-level RT DV
    "setSize"        = "Effort level",
    "type"           = "Task-type",
    "easyOffer"      = "Easy offer",         # also fixes "easyOffer z" -> "Easy offer z"
    "wm dev"         = "WM error",            # scale(wm_dev) -> WM performance-error covariate
    "wm rt"          = "WM RT",               # scale(wm_rt) -> "wm rt" -> "WM RT"
    "drug"           = "Drug"                # scale/factor drug -> "Drug [PROPRANOLOL]"
  )
  for (nm in names(label_map)) {
    table_result$knitr <- gsub(paste0("\\b", nm, "\\b"),
                               label_map[[nm]], table_result$knitr, perl = TRUE)
  }
  # "type dev" -> "Task-type" (my key), which the later "type" -> "Task-type" key then
  # re-hits (the output "Task-type" contains a bare "type"), yielding "Task-Task-type".
  # Collapse that one self-inflicted double-replacement back to the intended label.
  table_result$knitr <- gsub("Task-Task-type", "Task-type", table_result$knitr, fixed = TRUE)

  # The label replacements above are applied to the whole rendered table, so they also
  # hit the title/caption (e.g. "drug" -> "Drug" inside "Pupil drug effects"). Restore
  # the original title verbatim: recompute what the mapping did to the title, then swap
  # that mangled form back for the caller-supplied title.
  if (!is.null(title)) {
    mangled_title <- title
    for (nm in names(label_map))
      mangled_title <- gsub(paste0("\\b", nm, "\\b"), label_map[[nm]], mangled_title, perl = TRUE)
    if (!identical(mangled_title, title))
      table_result$knitr <- sub(mangled_title, title, table_result$knitr, fixed = TRUE)
  }

  # if table is long make scrollable file
    if (length(pred_labels) > 10){
  # Wrap the HTML content in a <div> container with scrolling styles
  scrollable_table <- paste0(
    '<div style="overflow-x:auto; max-height:400px; overflow-y:scroll;">',
    table_result$knitr,
    '</div>')
    cat(scrollable_table, file = paste0(file_name,'_scroll.html'))
    
  }
  
  # save
  cat(table_result$knitr, file = paste0(file_name,'.html'))
  
  # Capture a screenshot of the HTML file
  webshot::webshot(
    paste0(file_name,'.html'),
    file.path(save_dir, sprintf("%s.png", formatted_title))
  )
}
