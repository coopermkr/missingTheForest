#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' Missing the forest for the trees
#' Case study code
#' 
#' @author Cooper Kimball-Rhines
#' @date 2025-01-06
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

#' Data used in this case study is from two Forest Inventory and Analysis database
#' spreadsheets. The raw data is in NC_TREE and plot metadata is in NC_COND. I will
#' assume these data are saved into a subdirectory called "data".
#' This code walks through:
#' 1. Data manipulation and metadata joining
#' 2. Calculating diversity with rarefaction curves
#' 3. Parallel trends analysis
#' 4. Difference-in-Difference modelling
#' 5. Placebo test
#' 6. Data visualization

#### 1. Data manipulation ####
# These libraries are active throughout all parts of this code
library(tidyverse)
library(ggpubr)

# Start with the plot-level tree data
trees <- read.csv("data/NC_TREE.csv")

# Keep useful columns
filtrees <- trees |>
  select(CN, INVYR, UNITCD, COUNTYCD, PLOT, SUBP,
         TREE, STATUSCD, SPCD, SPGRPCD, RECONCILECD) |>
  filter(STATUSCD == 1)

# Fix year mistake (replace 9999 with 1999)
filtrees$INVYR <- gsub(9999, 1999, filtrees$INVYR)

# Save the intermediate for future use
#write.csv(filtrees, file = "data/NC_filtrees.csv")

# Load in filtered data if picking up back here
#filtrees <- read.csv("data/NC_filtrees.csv") |>
#  mutate(bind = paste(sep = "_", INVYR, COUNTYCD, PLOT))

# Now filter and join the metadata
condition <- read.csv("data/NC_COND.csv")

filtCond <- condition |>
  select(INVYR, UNITCD, COUNTYCD, PLOT,
         OWNCD, RESERVCD, OWNGRPCD) |>
  mutate(bind = paste(sep = "_", INVYR, COUNTYCD, PLOT))

ownTree <- merge(filtrees, filtCond)

# Create one code for all national land
unique(ownTree$OWNCD)# Everything less than 31 is federal/national
ownTree$OWNCD <- replace(ownTree$OWNCD, ownTree$OWNCD < 31, 1)
# 1 = federal, 31 = state, 46 = private


#### 2. Rarefaction ####
library(iNEXT)

# This program requires a matrix with species abundance by site
treeSummary <- ownTree |>
  filter(INVYR %in% c(1974, 1984, 1990, 2002, 2007, 2012, 2017, 2022),
         #INVYR %in% c(2003, 2004, 2005, 2006, 2009, 2010, 2011, 2013),
         #INVYR %in% c(2014, 2015, 2016, 2018, 2019, 2020, 2021),
         #NOTE: iNEXT can only calculate for so many plots, so we split the years up and run the code multiple times
         OWNCD %in% c(1, 31, 46)) |>
  group_by(INVYR, COUNTYCD, SPCD, OWNCD) |>
  summarize(count = n())

# plot total number of trees per plot to pick sample size
treeSummary |>
  ungroup(SPCD) |>
  group_by(COUNTYCD, INVYR) |>
  summarize(total = sum(count)) |>
  filter(total < 1000) |>
  ggplot(mapping = aes(x = total)) +
  geom_histogram()

# Looks like 200-250 is a good range for the median, we need to specify a sample size
# for iNEXT to set a rarefaction cutoff- so we will use 200.
m = 200

# Do q = 1 for Shannon (Hill number) with 100 bootstraps
treeShan <- iNEXT(as.data.frame(abundance[2:ncol(abundance)]), q = 1, datatype = "abundance", 
                  size = m, se = TRUE, conf = 0.95, nboot = 100)

# Pull out the rarefaction and separate the plot-ID columns
head(treeShan$iNextEst$coverage_based)

treeRare <- treeShan$iNextEst$coverage_based |>
  separate(col = Assemblage, into = c("year", "county", "owner"), sep = "-")

# Write the rarefaction numbers to a csv to save progress
write.csv(treeRare, file = "data/rarefiedShannon.csv")


#### 3. Parallel Trends ####
library(did)
library(ggpubr)

# Load and mutate data to indicate treatments
shanRare <- read_csv("data/rarefiedShannon.csv") |>
  mutate(treated = owner == '31' &
           year > 1997) |>
  # Place policy implementation at time zero
  mutate(timeFrom = year - 1997,
         owner = as.factor(owner),
         county = as.factor(county)) |>
  arrange(year)


shanRare |>
  filter(m == 200) |>
  #Method %in% c("Observed", "Estimated")) |>
  group_by(owner, year) |>
  summarize(mean = mean(qD),
            se = sd(qD)/n()) |>
  ggplot(mapping = aes(x = year, y = mean, color = as.factor(owner))) +
  geom_point()
# Looks mostly good- much more variance after the policy

# Now a real test:
# Estimate group-group time average treatment effects
# Note, this function has some weird object requirements so we have to mess with
# our data formats a bit to get it to work
ptShan <- shanRare |>
  # Select a sample size
  filter(m == 200,
         year %in% c('1974', '1984', '1990', '2002', '2007', '2012', '2017', '2022')) |>
  # Create a variable that tells when a county was first treated- 5 if it was, 0 if it was not
  mutate(firstTreat = as.numeric(owner == '31')*5,
         # Combine county and owner to create a unique ID for each entry,
         # then coerce it into a numeric so att_gt doesn't throw a fit
         idFac = as.factor(paste(county, owner, sep = "_")),
         idNum = as.numeric(idFac)) |>
  dplyr::select(qD, timeFrom, idNum, firstTreat, treated)

# Run the test
set.seed(1997)

trendAssess <- att_gt(yname = "qD",
                      tname = "timeFrom",
                      idname = "idNum",
                      gname = "firstTreat",
                      data = ptShan,
                      panel = TRUE,
                      bstrap = TRUE,
                      cband = TRUE)

summary(trendAssess) # This looks amazing

# Make a plot of the confidence intervals- Figure 2 in the paper
ptPlot <- tidy(trendAssess) |>
  ggplot(mapping = aes(x = time, y = estimate, color = time > 0)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high)) +
  labs(title = NULL,
       subtitle = NULL,
       tag = NULL,
       x = "Years From Implementation",
       y = "Effect Size",
       color = NULL) +
  theme_pubr(base_size = 20) +
  scale_color_manual(labels = c("Pre-amendment", "Post-amendment"), 
                     values = c('#648FFF', '#DC267F')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.5)


#### 4. Difference-in-Difference ####
library(fixest)
library(modelsummary)
library(performance)
library(broom)

# Build DID model with year and owner variables- we already took care of county by rarefying
treeDID <- shanRare |>
  filter(m == 200,
         # Note: we tested this with multiple sets of years and got similar models
         year %in% c('1974', '1984', '1990', '2002', '2007', '2012', '2017', '2022')) |>
  rename("Shannon" = qD,
         "Amendment " = treated,
         "Years From Implementation" = timeFrom,
         "Owner" = owner) |>
  feols(fml = Shannon ~ `Amendment ` | `Years From Implementation` + Owner)

check_model(treeDID) # Residuals looks good
r2(treeDID) # Fit is ok

# Assess coefficients
tidy(treeDID, conf.int = TRUE)
msummary(treeDID, stars = c('*' = .1, '**' = .05, '***' = .01))

# Effect of treatment is an increase of 3.054 on the Shannon Index at p < 0.01


#### 5. Placebo Test ####

# Placebo Test- only looking at the pretreatment data, set a placebo of the treatment
# occuring 10 years before it actually did
pl <- shanRare |>
  filter(m == 200,
         timeFrom < 0) |>
  mutate(treated = owner == '31' &
           timeFrom > -10) |>
  rename("Shannon" = qD,
         "Amendment " = treated,
         "Years From Implementation" = timeFrom,
         "Owner" = owner)

# Create a model with the placebo treatment
plDID <- feols(fml = Shannon ~ `Amendment ` | `Years From Implementation` + Owner,
               data = pl)

# Assess the placebo model like we did above
msummary(plDID, stars = c('*' = .1, '**' = .05, '***' = .01))
# There is no effect of treatment- placebo test passed.

# Print a summary table- Table 3 in the paper
msummary(list("Amendment" = treeDID, 
              "Placebo" = plDID), 
         title = "DID Model Summary",
         stars = c('*' = .1, '**' = .05, '***' = .01),
         statistic = "conf.int", conf_level = 0.95,
         gof_omit = "R2 Within|R2 Within Adj.|Std.Errors|FE: `Years From Implementation`|FE: Owner")

#### 6. Data Visualization ####
library(gridExtra)
library(grid)
library(ggpubr)

# Recalculate grand means for before and after policy implementation
didSeg <- shanRare |>
  filter(m == 200,
         year %in% c('1974', '1984', '1990', '2002', '2007', '2012', '2017', '2022')) |>
  # Calculate standard error from confidence intervals
  mutate(se = (qD.UCL - qD.LCL)/1.96) |>
  
  # Calculate mean and sum of squares for the standard error weighted by sample size
  group_by(owner == 31, timeFrom > 0) |>
  summarize(m = mean(qD),
            n = n(),
            v = sum(se^2+qD^2)/n() - m^2) |>
  mutate(stdev = sqrt(v),
         gse = stdev/sqrt(n))

## DID Plot
# Create a nice base plot showing the points for each plot in the background
didPlot <- shanRare |>
  mutate(`Land Owners` = owner == 31) |>
  filter(m == 200,
         year %in% c('1974', '1984', '1990', '2002', '2007', '2012', '2017', '2022')) |>
  
  ggplot(mapping = aes(x = timeFrom, y = qD,
                       color = `Land Owners`, shape = `Land Owners`)) +
  geom_point(position = position_jitter(width = .5, seed = 10), size = 1, alpha = 0.9) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', linewidth = 0.5) +
  scale_color_manual(labels = c("Non-State", "State-Owned"), 
                     values = c('#648FFF', '#DC267F')) +
  scale_shape_manual(labels = c("Non-State", "State-Owned"), 
                     values = c(0, 2)) +
  theme_pubr() +
  labs(x = "Years From Treatment",
       y = "Shannon Diversity Index") +
  theme(plot.title = element_text(hjust = 0.5))


# And lines showing the overall averages before and after treatment and error bars
didPlotSegs <- didPlot +
  # Unaffected Pretreatment
  geom_rect(xmin = -23, xmax = 0,
            ymin = didSeg$m[1] - didSeg$gse[1], 
            ymax = didSeg$m[1] + didSeg$gse[1],
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = min(timeFrom), xend = 0,
                   y = didSeg$m[1], yend = didSeg$m[1]),
               linewidth = 1, color = "#648FFF", linetype = "dashed") +
  
  #Unaffected Posttreatment
  geom_rect(xmin = 0, xmax = 25,
            ymin = didSeg$m[2] - didSeg$gse[2], 
            ymax = didSeg$m[2] + didSeg$gse[2],
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = 0, xend = max(timeFrom),
                   y = didSeg$m[2], yend = didSeg$m[2]),
               linewidth = 1, color = "#648FFF", linetype = "dashed") +
  
  # Affected Pretreatment
  geom_rect(xmin = -23, xmax = 0,
            ymin = didSeg$m[3] - didSeg$gse[3], 
            ymax = didSeg$m[3] + didSeg$gse[3],
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = min(timeFrom), xend = 0,
                   y = didSeg$m[3], yend = didSeg$m[3]),
               linewidth = 1, color = "#DC267F", linetype = "dashed") +
  
  # Affected Posttreatment
  geom_rect(xmin = 0, xmax = 25,
            ymin = didSeg$m[4] - didSeg$gse[4], 
            ymax = didSeg$m[4] + didSeg$gse[4],
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = 0, xend = max(timeFrom),
                   y = didSeg$m[4], yend = didSeg$m[4]),
               linewidth = 1, color = "#DC267F", linetype = "dashed") +
  
  # Add labels
  geom_text(x = 0, y = 15, label = "Non-State", color = "#648FFF") +
  geom_text(x = 0, y = 8.85, label = "State-Owned", color = "#DC267F") +
  labs(title = "Unadjusted Means", tag = "A")

didPlotSegs

# Now norm based on changes in the control to actually show Diff-In-Diff
# Subtract the change in the control from the post-treatment Impact group
# To calculate the combined standard error take the square root of the sum of the gse squares

didAdjusted <- didPlot +
  # Unaffected Pretreatment
  geom_rect(xmin = -23, xmax = 0,
            ymin = didSeg$m[1] - didSeg$gse[1], 
            ymax = didSeg$m[1] + didSeg$gse[1],
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = min(timeFrom), xend = 0,
                   y = didSeg$m[1], yend = didSeg$m[1]),
               linewidth = 1, color = "#648FFF", linetype = "dashed") +
  
  #Unaffected Posttreatment
  geom_rect(xmin = 0, xmax = 25,
            ymin = didSeg$m[2]-(didSeg$m[2]-didSeg$m[1]) - sqrt(didSeg$gse[1]^2 + 2*didSeg$gse[2]^2), 
            ymax = didSeg$m[2]-(didSeg$m[2]-didSeg$m[1]) + sqrt(didSeg$gse[1]^2 + 2*didSeg$gse[2]^2),
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = 0, xend = max(timeFrom),
                   y = (didSeg$m[2]-(didSeg$m[2]-didSeg$m[1])), 
                   yend = didSeg$m[2]-(didSeg$m[2]-didSeg$m[1])),
               linewidth = 1, color = "#648FFF", linetype = "dashed") +
  
  # Affected Pretreatment
  geom_rect(xmin = -23, xmax = 0,
            ymin = didSeg$m[3] - didSeg$gse[3], 
            ymax = didSeg$m[3] + didSeg$gse[3],
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = min(timeFrom), xend = 0,
                   y = didSeg$m[3], yend = didSeg$m[3]),
               linewidth = 1, color = "#DC267F", linetype = "dashed") +
  
  # Affected Posttreatment
  geom_rect(xmin = 0, xmax = 25,
            ymin = didSeg$m[4]-(didSeg$m[2]-didSeg$m[1]) - sqrt(didSeg$gse[4]^2 + didSeg$gse[1]^2 + didSeg$gse[2]^2), 
            ymax = didSeg$m[4]-(didSeg$m[2]-didSeg$m[1]) + sqrt(didSeg$gse[4]^2 + didSeg$gse[1]^2 + didSeg$gse[2]^2),
            fill = "lightgrey", color = "lightgrey", alpha = 0.1) +
  
  geom_segment(aes(x = 0, xend = max(timeFrom),
                   y = (didSeg$m[4]-(didSeg$m[2]-didSeg$m[1])), 
                   yend = (didSeg$m[4]-(didSeg$m[2]-didSeg$m[1]))),
               linewidth = 1, color = "#DC267F", linetype = "dashed") +
  
  # Add Labels
  geom_text(x = 0, y = 15.3, label = "Non-State", color = "#648FFF") +
  geom_text(x = 0, y = 12.2, label = "State-Owned", color = "#DC267F") +
  labs(title = "Adjusted Means", tag = "B")

didAdjusted

# Create panel plot- Figure 3 in the paper
grid.arrange(didPlotSegs, didAdjusted, ncol = 2, 
             top=textGrob("State Forest Diversity is Stable in a Backdrop of Decline",
                          gp=gpar(fontsize = 20),
                          just = "centre"))

