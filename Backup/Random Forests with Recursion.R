# -------------------------------------------------------------------------
# Implement Random Forests algorithm using recursion
# Train a classification tree to predict heart disease using out-of-bag error
# Data: South African Heart Disease <https://rdrr.io/cran/ElemStatLearn/man/SAheart.html>
# 
# Ben Claassen
# -------------------------------------------------------------------------


# NOTES -------------------------------------------------------------------
## All character classes in data set must be binary
## Outcome variable must be binary and 0/1 in format


# Libraries ---------------------------------------------------------------
library(tidyverse)


# -------------------------------------------------------------------------
# Functions ---------------------------------------------------------------

# -------------------------------------------------------------------------
## Function -> Get predicted outcome of input branch (most common unique value)
getmode <- function(v)
{
  uniqv <- unique(v) # List unique values in terminal branch
  uniqv[which.max(tabulate(match(v, uniqv)))] # Return most common of unique values in terminal branch
}


# -------------------------------------------------------------------------
## Function -> Create tree
branchTree <- function(inputList, usedVariables = NULL, maxObsPerBranch = 10)
{
  
  if(length(inputList) > maxObsPerBranch) # Confirm that there are "enough" obs to continue branching
  {
    # All tree entries are either indicator entries, or lists of unparsed obs
    #   So if not indicator entry, then length > 1, so split into new branches
    tmpList <- inputList
    
    if(length(usedVariables) > 0) # Don't re-check variables that have already been used to split the tree
    {
      variablesList <- names(heartDat)[-usedVariables]
    } else {
      variablesList <- names(heartDat)
    }
    
    # Test Residual Sum Of Squares (RSS) for median cut on each variable
    rssStorage <- matrix(, length(variablesList), 2) # Rows: Number of variables to check; Cols: [colNumber],[entropy]
    varCount <- 0 # Initialize a variable counter
    
    # Iterate over variables
    for(varNames in variablesList)
    {
      varCount <- varCount + 1
      # Set column for current varName
      v = which(names(heartDat) == varNames)
      
      ## Calculate individual category entropies
      if(class(heartDat[,v]) == "character") # String variables
      {
        
        regionRSS <- 0
        for(class in unique(heartDat[tmpList,v])) # Iterate over subsets in string variable
        {
          tmpRows <- tmpList[which(heartDat[tmpList,v] == class)] # Extract rows with current subset
          regionAvg <- mean(heartDisease[tmpRows]) # Calculate average incidence of heart disease across that subset
          regionRSS <- regionRSS + sum((heartDisease[tmpRows] - regionAvg)^2) # Calculate the RSS
        }
        # Store entropy
        rssStorage[varCount,1] <- v
        rssStorage[varCount,2] <- regionRSS
        
      } else { # Numeric variable
        
        # Split variable at the median value
        tmpRows <- tmpList[which(heartDat[tmpList,v] >= median(heartDat[tmpList,v]))] 
         # Store the average incidence of heart disease for observations above and below the median
        regionAvg_1 <- mean(heartDisease[tmpRows])
        regionAvg_2 <- mean(heartDisease[-tmpRows])
        # Calculate the RSS
        regionRSS <- sum((heartDisease[tmpRows] - regionAvg_1)^2)
        regionRSS <- regionRSS + sum((heartDisease[-tmpRows] - regionAvg_2)^2)
        
        # Store entropy
        rssStorage[varCount,1] <- v # Data column number
        rssStorage[varCount,2] <- regionRSS # Data column RSS
        
      }
    } # End variable check list {v, varNames}
    
    
    # Make cut (use first if matching entropies)
    tmpCutVarCol <- rssStorage[which(rssStorage[,2] == min(rssStorage[,2]))[1], 1] # Find minimum RSS/entropy to cut on. NOTE minimizing entropy maximizes similarity in classification in the resulting sub-tree
    usedVariables <- c(usedVariables, tmpCutVarCol) # Concatenate variable used for split to previous var list
    
    # print(rssStorage)
    # print(paste0("~ ", tmpCutVarCol, " ~"))
    
    ## Split list
    if(class(heartDat[,tmpCutVarCol]) == "character")
    {
      #Set variable cut to "C" for categories
      tmpCutValue <- "C" # Set storage value
      firstClass <- unique(heartDat[tmpList,tmpCutVarCol])[1] # Set first class in list for splitting
      
      # Make split
      tmpSplit_A <- tmpList[which(heartDat[tmpList,tmpCutVarCol] == firstClass)]
      tmpSplit_B <- tmpList[which(heartDat[tmpList,tmpCutVarCol] != firstClass)]
      
      
    } else { #numeric variable
      
      tmpCutValue <- median(heartDat[tmpList,tmpCutVarCol]) # Set storage value
      
      # Make split
      tmpSplit_A <- tmpList[which(heartDat[tmpList,tmpCutVarCol] >= tmpCutValue)]
      tmpSplit_B <- tmpList[which(heartDat[tmpList,tmpCutVarCol] < tmpCutValue)]
      
    }
    
    
    # Store split information and sub-trees
    inputList <- list(tmpCutVarCol, tmpCutValue, tmpSplit_A, tmpSplit_B)
     
    # Recursively run branching algorithm on each sub-tree
    inputList[[3]] <- branchTree(inputList[[3]], usedVariables)
    inputList[[4]] <- branchTree(inputList[[4]], usedVariables)
    
   
    return(inputList)
    
  } else { # Else, if there insufficient additional obs on this branch, return the input list
    
    return(inputList)
  }
}


# -------------------------------------------------------------------------
## Function -> Predict outcome given a tree and an out-of-bag observation
#  - Follow branching tree using data from given obs to find expected outcome
chooseTreeBranch <- function(tmpPoint, tmpTree)
{
  
  if(is.list(tmpTree)) # Check if [tmpTree] is a list, i.e. has multiple branches
  {
    if(tmpTree[[2]] == "C") # Check if current split variable is a class
    {
      if(tmpPoint[tmpTree[[1]]] == unique(heartDat[, tmpTree[[1]]])[1]) # If value in obs matches first of unique chars, split left, else split right
      {
        nextTree <- 3
      } else {
        nextTree <- 4
      }
      
    } else { # Else, current split variable is numeric
      if(tmpPoint[,tmpTree[[1]]] >= tmpTree[[2]]) # If tmpPoint has value above median for current top variable, split left, else split right
      {
        nextTree <- 3
      } else {
        nextTree <- 4
      }
    }
  } else { # Else, if [tmpTree] is not a list, then you've hit end of a branch, so output most likely outcome
    return(getmode(heartDisease[tmpTree]))
  }
  
  # If sub-tree branch ([nextTree]) is itself branching, recursively call [chooseTreeBranch] until a solution is reached
  if(length(tmpTree[[nextTree]]) > 0)
  {
    chooseTreeBranch(tmpPoint, tmpTree[[nextTree]])
  } else {
    # Else, the selected tree branch is empty. Then get outcome values for current branch and all sub-branches without splitting further
    output <- getFinalNodeVals(tmpTree, NULL)
    
    # print("WARNING: DEGENERATE TREE BRANCH - Length 0")
    
    return(getmode(heartDisease[output])) # Return the computed outcome
  }
  
  
}


# -------------------------------------------------------------------------
## Function -> If final branch selected is empty, approximte outcome by collating all answers on next-best branch, including all sub-branches
getFinalNodeVals <- function(tmpBadTree, outputVals)
{
  if( is.list(tmpBadTree[[3]]) ) # Find which branch is not degenerate
  {
    outputVals <- c(outputVals, getFinalNodeVals(tmpBadTree[[3]], outputVals)) # Recursively iterate over said branch to get all observation numbers
  } else {
    outputVals <- c(outputVals, tmpBadTree[[3]])
  }
  
  if( is.list(tmpBadTree[[4]]) ) # Find which branch is not degenerate
  {
    outputVals <- c(outputVals, getFinalNodeVals(tmpBadTree[[4]], outputVals)) # Recursively iterate over said branch to get all observation numbers
  } else {
    outputVals <- c(outputVals, tmpBadTree[[4]])
  }
  
  # Return list of next-best branch answers
  return(outputVals)
}



# -------------------------------------------------------------------------
# Read in data ------------------------------------------------------------
# setwd("~/Documents/Reno/Reno_2019-Spring/MachineLearning/Homework9/Data")
heartDat <- read.csv("SouthAfricanHeartDisease.csv", stringsAsFactors = F)
heartDisease <- heartDat %>% select(chd) %>% unlist() # Set the training key values [chd]
heartDat <- heartDat %>% select(!chd) # Drop the key


# Set up algorithm --------------------------------------------------------
numTrees <- 50 # Set number of trees for random forest implementation
sampleDat <- vector(mode = "list", length = numTrees) # Initialize storage vector for sampled data


# Fill sample list
# set.seed(4321)
for(d in 1:numTrees)
{
  sampleDat[[d]] <- sample(1:462, replace = TRUE) # Sample with replacement observation numbers for each tree
}

# Initialize storage array for confusion matrices; 2x2 grid for each tree
confusionMatrixArray <- array(0, dim = c(2,2,numTrees))



# Run algorithm -----------------------------------------------------------

for(t in 1:numTrees)
{
  tmpSampleList <- sampleDat[[t]] # Assign sample list
  # print(length(unique(tmpSampleList)) / 462) # Print how many observations are included in sample
  
  # Create decision tree
  tree <- branchTree( tmpSampleList ) 
  
  # Set list of 'out of bag' obs
  outOfBagObs <- c(1:462)[-unique(tmpSampleList)]
  sumOfErrors <- 0
  
  # Iterate over out of bag obs and calculate predicted outcome
  for(b in outOfBagObs)
  {
    # print(b)
    # errors <- chooseTreeBranch(heartDat[b,], tree)
    # print(chooseTreeBranch(heartDat[b,], tree))
    
    # Predict outcome
    predictedValue <- chooseTreeBranch(heartDat[b,], tree)
    
    # Check if outcome is correct
    if(predictedValue == 1)
    {
      if(heartDisease[b] == predictedValue)
      {
        confusionMatrixArray[2,2,t] <- confusionMatrixArray[2,2,t] + 1 # Increase True Positive by 1
      } else {
        confusionMatrixArray[2,1,t] <- confusionMatrixArray[2,1,t] + 1 # Increase False Positive by 1
      }
    } else if(predictedValue == 0)
    {
      if(heartDisease[b] == predictedValue)
      {
        confusionMatrixArray[1,1,t] <- confusionMatrixArray[1,1,t] + 1 # Increase True Negative by 1
      } else {
        confusionMatrixArray[1,2,t] <- confusionMatrixArray[1,2,t] + 1 # Increase False Negative by 1
      }
    }
    
  }
  
  # Print confusion matrix
  # print(confusionMatrixArray[,,t])
  
  # Print F1 score
  # (2 * true positives) / ( (2 * true positives) + (false positives) + (false negatives) )
  print(
    round(
      (2*confusionMatrixArray[2,2,t]) / ( 2*confusionMatrixArray[2,2,t] + confusionMatrixArray[2,1,t] + confusionMatrixArray[1,2,t] )
    , 3)
  )
  
}

# Print average F1 score (1 = perfect prediction; 0 = precision or sensitivity are 0%)
mean(
  (2*confusionMatrixArray[2,2,]) / ( 2*confusionMatrixArray[2,2,] + confusionMatrixArray[2,1,] + confusionMatrixArray[1,2,] )
)


