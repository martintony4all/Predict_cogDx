# =============================================
# MULTIPLE IMPUTATION & EXPLORATORY ANALYSIS
# =============================================

# ----------------------------
# 1. LOAD LIBRARIES
# ----------------------------

required_packages <- c("mice", "MASS", "VIM", "lattice", "dplyr", 
                       "ggplot2", "rlang")

missing <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if (length(missing) > 0) install.packages(missing)

suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
})

# ----------------------------
# 2. LOAD & CLEAN DATA
# ----------------------------

df <- read.csv("/Users/martinnwadiugwu/Downloads/RAC_data_sheet1.csv", 
               stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN", ".", " "))

df <- df %>%
  mutate(
    Age = as.numeric(gsub("\\+", "", Age)),
    Gender = factor(Gender),
    APOE_Status = factor(APOE_Status),
    Pathology = factor(Pathology),
    Diagnosis = factor(Diagnosis),
    Highest_Education_Level = factor(Highest_Education_Level),
    Braak = factor(Braak, levels = 0:6, ordered = TRUE),
    Years_Education = factor(Years_Education, levels = 9:25, ordered = TRUE),
    last_MMSE = factor(last_MMSE, levels = 0:30, ordered = TRUE)
  )



# ----------------------------
# 3. MISSINGNESS DIAGNOSTICS
# ----------------------------

cat("Missing data pattern:\n")
md.pattern(df)

cat("Missing data summary:\n")
aggr(df, numbers = TRUE, prop = FALSE)



# ----------------------------
# 4. IMPUTATION SETUP
# ----------------------------

meth <- make.method(df)
pred <- make.predictorMatrix(df)

if ("ID" %in% names(df)) {
  meth["ID"] <- ""
  pred[, "ID"] <- 0
  pred["ID", ] <- 0
}

# Assign methods
meth["Age"] <- "pmm"
meth["PMI.h."] <- "pmm"
meth["KCNH8"] <- "pmm"
meth["ZEB2"] <- "pmm"
meth["ST18"] <- "pmm"
meth["ZNF536"] <- "pmm"
meth["Gender"] <- "logreg"
meth["Diagnosis"] <- "polyreg"
meth["Pathology"] <- "polyreg"
meth["Highest_Education_Level"] <- "polyreg"
meth["APOE_Status"] <- "polyreg"
meth["structure"] <- "polyreg"
meth["Braak"] <- "polr"
meth["Years_Education"] <- "polr"
meth["last_MMSE"] <- "polr"


# ----------------------------
# 5. RUN MULTIPLE IMPUTATION
# ----------------------------

set.seed(123)
imp <- mice(df, m = 5, method = meth, predictorMatrix = pred)

#colSums(is.na(df)) == nrow(df)

# Just the numerical summary, no plot
#aggr(df, plot = FALSE, numbers = TRUE, prop = FALSE)

#table(df$APOE_Status, useNA = "ifany")
#table(df$Diagnosis, useNA = "ifany")

#df <- df[, colSums(is.na(df)) < nrow(df)]

#meth <- make.method(df)
#pred <- make.predictorMatrix(df)

#imp <- mice(df, m = 5, method = meth, predictorMatrix = pred)


# ----------------------------
# 6. VISUALIZATION
# ----------------------------

# Density plots for numeric variables
numeric_vars <- c("Age", "PMI.h.")
for (v in numeric_vars) {
  for (i in 1:imp$m) {
    imp_data <- complete(imp, i)
    vals <- imp_data[[v]]
    if (sum(!is.na(vals)) >= 2) {
      print(
        ggplot(imp_data, aes(x = !!sym(v))) +
          geom_density(fill = "skyblue", alpha = 0.5) +
          labs(title = paste("Density plot of", v, "- Imputation", i),
               x = v, y = "Density") +
          theme_minimal()
      )
    } else {
      message(paste("Skipping", v, "in imputation", i, "- not enough data"))
    }
  }
}


# Bar plots for categorical variables
categorical_vars <- c("APOE_Status", "Braak", "Diagnosis", "Gender", "Pathology", "last_MMSE", "Highest_Education_Level", "Years_Education", "structure")
for (v in categorical_vars) {
  for (i in 1:imp$m) {
    imp_data <- complete(imp, i)
    vals <- imp_data[[v]]
    if (length(unique(na.omit(vals))) >= 2) {
      print(
        ggplot(imp_data, aes(x = !!sym(v))) +
          geom_bar(fill = "lightgreen") +
          labs(title = paste("Count plot of", v, "- Imputation", i),
               x = v, y = "Count") +
          theme_minimal()
      )
    } else {
      message(paste("Skipping", v, "in imputation", i, "- not enough levels"))
    }
  }
}

# ----------------------------
# 7. SAVE IMPUTED DATA
# ----------------------------

completed_data <- complete(imp, 1)
write.csv(completed_data, "/Users/martinnwadiugwu/Downloads/RAC_data_imputed_new_new11.csv", row.names = FALSE)

# ----------------------------
# 8. OPTIONAL: MODEL FITTING
# ----------------------------

# Check Braak levels before modeling
if (nlevels(df$Braak) >= 2) {
  fit <- with(imp, lm(PMI.h. ~ Age + Gender + Diagnosis + APOE_Status + Braak))
} else {
  warning("⚠️ Braak has fewer than 2 levels. Excluding from model.")
  fit <- with(imp, lm(PMI.h. ~ Age + Gender + Diagnosis + APOE_Status))
}

pooled <- pool(fit)
summary(pooled)

