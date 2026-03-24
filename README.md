# Statistical Methods TP02


## 🚀 Getting Started (first time using this repo)
To ensure all file paths and environments load correctly, always open the **`Statistical_Methods_A1.Rproj`** file first.

1.  **Activate Environment:** Opening the project should automatically activate `renv`. If it doesn't, run `renv::activate()` in your R console.
2.  **Install Dependencies:** Run the following command in your R console to install the required packages: `renv::restore()`
3.  **Run Analysis:** Open and knit **`main.Rmd`** to reproduce the full analysis and generate the final report.

---

## 🔄 Daily Workflow

### After Pulling Updates
Whenever you `git pull` new changes from the repository, the package versions might have changed. To stay in sync, run: `renv::restore()`.

### Before Committing Changes
If you've installed or updated any packages, run the following to save the current state of your environment in `renv.loc` before committing: `renv::snapshot()`.
