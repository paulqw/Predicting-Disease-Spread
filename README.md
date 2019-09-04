# **DengAI - Predicting Disease Spread**  
  
Date: 03/09/2019  
  
## Project summary
**Project description:** Participating the data-science competitions 'DengAI: Predicting Disease Spread' on www.drivendata.org in a one-days hackathon.  
Dengue fever, occurring in tropical and sub-tropical regions around the globe, in its mild course of disease is similar to a flu - fever, rash, and muscle and joint pain - but can, however, be extremely life-threatening in its sever peculiarity. Therefore, the competition seeks to accurately predict the number of weekly dengue fewer cases in two highly affected places - San Juan(sj), Puerto Rico and Iquitos(iq), Peru.
As dengue is carried by mosquitoes its transmission dynamics are highly related to the climate conditions in the particular region; especially, temperature, precipitation, and vegetation. The environmental data available within the competition has been collected by various U.S. Federal Government agencies (e.g. Centers for Disease Control and Prevention, the National Oceanic and Atmospheric Administration)  
The competition allows for any platform and ML-technique but evaluates the final result based on the mean absolute error (MAE). In this analysis the statistical computing programming language is used applying prophet algorithm for time series and various machine learning techniques such as support vector machine, tree-based techniques, and neural networks. 
  
## Technical approach of analysis  
  
**1. Exploration and preparation of data**
- Set time-zone
- Identify and tread/eliminate NAs
- Create time series
- Normalize data
- Conduct correlation analysis and identify relevant features
- Conduct cross-correlation analysis and identify relevant lags of features
    
**2. Feature-Engineering**
- Create lagged data of relevant features
  
**3. Data splitting**
- Split ratio: 75/25
  
**4. Train models, evaluate, forecast**
- Select a bunch of algorithms
- Train models using caret::caretList and automatic parameter tuning
- Create and test model ensemble
- Decide on model to be used for final decision based on performance metrics
- Tune hyperparameters for the final model
  
**5. Forecast and submit results**
  