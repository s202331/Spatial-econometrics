---
title: "Short-Term Forecasting and Simulation of Bitcoin Price Based on Market, Macro and Sentiment Indicators (2023–2025)"
author: "Piotr Łukowski, Tomasz Kotliński, Jakub Gużewski"
date: "`r Sys.Date()`"
output:
  pdf_document:
    latex_engine: xelatex
    toc: true
    number_sections: true
---


# Introduction

This report aims to develop an econometric model for short-term forecasting of Bitcoin price using a set of explanatory variables related to financial markets, macroeconomic performance, and investor sentiment. The analysis is based on monthly data from 2023 to 2025, and it includes exploratory data analysis, correlation testing, stationarity verification, model estimation, and short-term prediction. Also We chose this project topic because we’re genuinely into investing some people collect stamps, we collect risk. Where others panic, we grab popcorn. For us, making sense of market chaos isn’t just academic it’s our idea of a good time.

# Data Description

The dataset contains monthly observations on the following variables:

BITCOIN – price of Bitcoin in USD,

CPI – Consumer Price Index (inflation measure),

VIX – volatility index (sentiment/fear indicator),

AMD – stock price of AMD (related to crypto mining hardware),

GOLD – gold price as a traditional safe-haven asset.

# Loading packages and libraries

```{r setup, include=FALSE}
library(tidyverse)
library(lubridate)
library(ggplot2)
library(car)
library(MASS)
library(forecast)
library(scales)
library(strucchange)
library(tseries)
library(readxl)
library(urca)
library(lmtest)
```
# Import Dataset

```{r}
df <- readxl::read_excel("C:/Users/gozia/Desktop/Piotrek/Prognozowanie i Symulacja/PIS.xlsx") %>% 
  mutate(Data = as.Date(Data))
```

# Original series

This plot visualizes the original time series of Bitcoin prices in USD over the period from January 2023 to early 2025. The data is plotted on a monthly frequency and represents untransformed values. This visualization helps to examine the raw trend and volatility in Bitcoin's market price before any differencing or transformation is applied for stationarity checks or modeling purposes.

From the graph, it is evident that the Bitcoin price exhibits significant upward movement with notable fluctuations. Especially visible are phases of rapid growth, followed by short-term corrections. The steep rise in early 2024 may indicate a market boom or a significant macroeconomic or sentiment-driven shift.

This plot serves as the foundation for further time series analysis, such as stationarity testing, differencing, and forecasting. The aim is to identify structural patterns, potential mean-reverting behavior, and prepare the data for rigorous econometric modeling.

```{r}

ggplot(df, aes(x = Data, y = BITCOIN)) + 
geom_line() + 
ggtitle("Bitcoin Prices") + 
scale_y_continuous(labels = comma) + 
theme_minimal()
```

# First difference

To prepare the Bitcoin price series for time series modeling, a first difference transformation is performed to address the presence of non-stationarity. This step is fundamental in econometric modeling, as many forecasting techniques require the underlying data to exhibit constant statistical properties over time. The plot of the first differenced series illustrates short-term changes in Bitcoin price, removing level effects and revealing underlying volatility dynamics. The data shows mean-reverting behavior around zero, which is consistent with stationarity.

The autocorrelation function ACF and partial autocorrelation function PACF plots provide further diagnostic insights. The ACF shows significant autocorrelation at lag 1, which rapidly diminishes, suggesting that recent price changes are correlated but that longer memory is weak. The PACF also shows a few short-term spikes, implying that only a limited number of autoregressive terms may be necessary to model the dynamics. These results justify the use of low-order ARIMA models and help specify their structure.

To statistically assess stationarity, two complementary tests are applied. The Augmented Dickey-Fuller ADF test fails to reject the null hypothesis of a unit root for both the original and differenced series - p-values of approximately 0.38 and 0.32, respectively. Although the p-value for the differenced series is lower, it is not below conventional significance thresholds, which weakens formal evidence of stationarity. In contrast, the Kwiatkowski–Phillips–Schmidt–Shin test, which reverses the null hypothesis, strongly supports stationarity for the differenced series - test statistic = 0.1423, well below all critical values.

The discrepancy between the two tests is not uncommon and reflects differences in testing approach. However, when interpreted jointly with visual evidence and the autocorrelation structure, the overall conclusion supports that the differenced Bitcoin series can be considered sufficiently stationary for modeling purposes. An automated ARIMA model is fit to the differenced data using the auto.arima function, producing a model specification that will be evaluated further using residual diagnostics in the next step.

```{r}
ts.plot(diff(df$BITCOIN), main = "First Difference of BITCOIN", ylab = "DPrice", col = "blue")

par(mfrow = c(1,2)) 
acf(diff(df$BITCOIN), main = "ACF of First Difference")
pacf(diff(df$BITCOIN), main = "PACF of First Difference") 
par(mfrow = c(1,1))

library(tseries) 
adf.test(df$BITCOIN)           # for original series
adf.test(diff(df$BITCOIN))    # for differentiated series

library(urca) 
summary(ur.kpss(df$BITCOIN))
summary(ur.kpss(diff(df$BITCOIN)))

library(forecast)
model_diff4 <- auto.arima(diff(df$BITCOIN))
resid_model <- residuals(model_diff4)
```

# Residual histogram

The histogram presented here illustrates the distribution of residuals obtained from the ARIMA model applied to the differenced Bitcoin price series. Analyzing residuals is a critical step in validating the adequacy and reliability of time series models, as it provides insight into whether the model captures the underlying structure of the data or leaves systematic patterns unexplained.

In this case, the residuals appear to be symmetrically distributed around zero, which is a desirable characteristic indicating that the model does not systematically over- or under-predict the target variable. The histogram approximates a bell-shaped curve, suggesting that the residuals may be normally distributed — an important assumption for constructing valid confidence intervals and performing hypothesis testing.

However, there are a few outliers present in both tails of the distribution, which may reflect occasional large prediction errors. These extreme values do not appear excessively frequent, but they merit further examination using formal normality tests and visual tools such as Q-Q plots.

Overall, the residual histogram supports the assumption of approximate normality and the absence of strong model misspecification, indicating that the ARIMA model provides a reasonably good fit to the differenced series.

```{r}
hist(resid_model, main = "Histogram of Residuals", col = "lightgray", breaks = 20)
```

# QQ-plot

The Q-Q (quantile-quantile) plot provides a graphical diagnostic for evaluating the assumption that the residuals from the ARIMA model are approximately normally distributed. In this plot, the empirical quantiles of the model residuals are plotted against the theoretical quantiles of a standard normal distribution. If the residuals are normally distributed, the points should lie approximately along the reference line shown in red. In this case, the central portion of the data aligns well with the theoretical normal distribution, suggesting that the bulk of the residuals exhibit near-normal behavior. However, deviations from the line are observed at both extremes, indicating the presence of mild heavy tails or outliers in the residual distribution. These deviations imply that the residuals may not be perfectly normal, particularly in the tails, which could slightly affect inference based on confidence intervals or prediction intervals. Nevertheless, the residuals do not exhibit extreme departures from normality, and the model may still be considered adequate for practical purposes. The Q-Q plot complements the histogram by providing a more sensitive diagnostic for tail behavior, and together, they support the conclusion that the ARIMA model residuals are approximately, though not perfectly, normally distributed.

```{r}
qqnorm(resid_model) 
qqline(resid_model, col = "red")
```

# ACF of residuals

In the presented ACF plot, the majority of the autocorrelation coefficients lie within the 95% confidence bounds (indicated by the dashed horizontal lines), suggesting that there is no statistically significant autocorrelation remaining in the residuals. This supports the assumption that the residuals resemble white noise — a fundamental requirement for model validity. Although a few minor spikes appear at higher lags, they are relatively small and do not indicate a clear violation of the independence assumption. This implies that the model has sufficiently accounted for the temporal dependencies in the data and that no strong serial structure remains unexplained.

```{r}
acf(resid_model, main = "ACF of Residuals")
```

# Combined Stationarity Diagnostics for First Difference

The top panel shows the time series of the first differences, which appears to fluctuate around a constant mean without a persistent trend. This behavior is indicative of weak stationarity, where the mean and variance remain stable over time. The bottom-left panel presents the ACF of the differenced series. The autocorrelations quickly decline and remain mostly within the 95% confidence bounds, which suggests that there is no significant serial correlation left in the transformed series. The bottom-right panel shows the PACF, which helps identify any remaining autoregressive structure. The absence of dominant spikes in both ACF and PACF plots supports the hypothesis that the differenced series behaves like white noise. Together, these plots confirm that the first difference transformation has effectively removed non-stationarity from the original series. The results validate that the transformed series is suitable for ARIMA modeling and that no strong autocorrelation structure remains to be accounted for.

```{r}
library(forecast)
tsdisplay(diff(df$BITCOIN), main = "Stationarity Diagnostics – First Difference")
```

```{r}
df$Data <- as.Date(df$Data)
df <- df[order(df$Data), ]

n <- nrow(df)
cat("Number of months:", n, "\n")

h <- 6
if (h >= n) {
  h <- n - 1
  cat("Not enough data for 6 test months, setting h to:", h, "\n")
}

btc <- ts(df$BITCOIN,
          start = c(as.numeric(format(min(df$Data), "%Y")),
                    as.numeric(format(min(df$Data), "%m"))),
          frequency = 12)

btc_train <- window(btc, end = c(time(btc)[n - h]))
btc_test <- window(btc, start = c(time(btc)[n - h + 1]))

cat("Training set has", length(btc_train), "points\n")
cat("Test set has", length(btc_test), "points\n")
```

# Arima Model Estimation Results

This section presents the output of an ARIMA(3,2,5) model fitted to the training subset of the Bitcoin time series. The model includes three autoregressive terms, two levels of differencing, and five moving average terms, selected based on prior stationarity analysis and model selection procedures.The reported coefficient estimates for the AR and MA terms indicate the dynamic structure captured by the model. While not all coefficients are individually statistically significant due to relatively large standard errors, the overall model fit can be judged by information criteria and error metrics. The residual variance (sigma\^2) is approximately 54.6 million, and the log-likelihood value is -163.24. The Akaike Information Criterion (AIC = 344.47), corrected AIC (AICc = 374.47), and Bayesian Information Criterion (BIC = 351.42) provide a basis for model comparison, with lower values indicating a better balance between model fit and complexity. Model performance is further evaluated using standard training set error metrics: ME (Mean Error): 522.15 RMSE (Root Mean Squared Error): 4926.91 MAE (Mean Absolute Error): 3132.38 MAPE (Mean Absolute Percentage Error): 7.28% MASE (Mean Absolute Scaled Error): 0.09196

These values indicate that the model captures a significant portion of the structure in the training data, with relatively low scaled and percentage errors. Additionally, the first-lag autocorrelation of residuals (ACF1 = -0.0626) suggests minimal remaining serial correlation, which supports the adequacy of the model.

Overall, the ARIMA(3,2,5) model provides a statistically coherent and reasonably accurate representation of the underlying Bitcoin price dynamics within the training period.

```{r}
model <- Arima(btc_train, order = c(3,2,5))
summary(model)
```

# Forecast Accuracy Evaluation

The table presents forecast accuracy metrics for the ARIMA(3,2,5) model across the training and test datasets. In-sample results show strong performance, with low RMSE (4926.91) and MAPE (7.28%). The MASE value of 0.09 confirms that the model outperforms a naive benchmark. Out-of-sample performance is weaker, with RMSE increasing to 12018.87 and MAPE rising to 12.16%, indicating reduced forecast accuracy. Theil’s U statistic is 0.997, suggesting that the model performs only slightly better than a naive forecast in the test set. The increase in ACF1 for forecast errors (0.37) may indicate some remaining autocorrelation.

```{r}
forecast_test <- forecast(model, h = h)
acc <- accuracy(forecast_test, btc_test)
print(acc)
```

# ARIMA Forecast Visualization with Confidence Intervals

This plot displays the in-sample fit and out-of-sample forecast generated by the ARIMA(3,2,5) model for Bitcoin prices, along with associated 80% and 95% confidence intervals. The black line represents the original observed data, while the blue line shows the model's fitted values and future forecasts. The red segment corresponds to the actual test set data, allowing direct visual comparison with forecasted values. The shaded areas depict the model's forecast uncertainty. The darker band corresponds to the 80% confidence interval, and the lighter band indicates the 95% confidence interval. These intervals capture the expected variability in future Bitcoin prices under the model assumptions. Overall, the fitted trajectory tracks the historical series closely. The forecast lies within a widening confidence region, reflecting increased uncertainty over time. The actual test values remain largely within the upper range of the forecast intervals, suggesting reasonable alignment between predicted and observed trends despite inherent market volatility.

```{r}
plot(btc, xlim = c(start(btc)[1], end(btc)[1] + length(forecast_test$mean)/6),
     ylim = range(btc, forecast_test$lower, forecast_test$upper),
     main = "ARIMA(3,2,5) Forecast with Fit and Confidence Bands",
     ylab = "BITCOIN", xlab = "Time", col = "black", lwd = 1)
time_forecast <- time(forecast_test$mean)

polygon(c(time_forecast, rev(time_forecast)),
        c(forecast_test$upper[,1], rev(forecast_test$lower[,1])),
        col = rgb(0,0,1,0.2), border = NA)

polygon(c(time_forecast, rev(time_forecast)),
        c(forecast_test$upper[,2], rev(forecast_test$lower[,2])),
        col = rgb(0,0,1,0.1), border = NA)

fitted_and_forecast <- ts(c(fitted(forecast_test$model), forecast_test$mean),
                          start = start(btc_train),
                          frequency = frequency(btc_train))
lines(fitted_and_forecast, col = "blue", lwd = 2)
lines(btc, col = "black")
lines(btc_test, col = "red", lwd = 2)

legend("topleft", legend = c("Original Data", "Test Data", "Fitted + Forecast", "95% Forecast CI"),
       col = c("black", "red", "blue", rgb(0,0,1,0.2)), lwd = c(1, 2, 2, NA), pch = c(NA, NA, NA, 15),
       pt.cex = 2, bty = "n")
```

# ETS Model Forecast and Performance

An Exponential Smoothing State Space Model (ETS) of type (M,N,N) was fitted to the Bitcoin time series using the ets() function. This model assumes multiplicative errors, no trend, and no seasonality. The smoothing parameter ALfa is estimated to be close to 1, indicating that the model places nearly full weight on the most recent observations. Model fit statistics include AIC = 504.97 and BIC = 508.50. These values are higher than those obtained from the ARIMA model, suggesting a relatively weaker in-sample fit. Training set error metrics confirm this: RMSE = 8095.12 MAE = 5500.70 MAPE = 10.17% MASE = 0.145 The ACF1 value (−0.08) indicates low autocorrelation in residuals, which supports the validity of the model.

The forecast plot shows the expected Bitcoin price trajectory over a 4-month horizon, with increasing uncertainty illustrated through expanding confidence intervals. Unlike the ARIMA forecast, the ETS forecast appears flatter, reflecting the model's strong dependence on recent levels without explicitly capturing underlying trend dynamics.

```{r}
model_ets <- ets(btc)
summary(model_ets)
forecast_ets <- forecast(model_ets, h = 4)
plot(forecast_ets, main = "Bitcoin Forecast - ETS")
```

# DESCRIBE MODELS

# Linear model of Bitcoin with CPI, VIX, AMD, GOLD

The model explains approximately 86.5% of the variance in Bitcoin prices (R² = 0.865), with an adjusted R² of 0.836, indicating a strong in-sample fit. The model is statistically significant overall (F-statistic = 30.4, p \< 0.001). Among the predictors, only the gold price is highly significant (p \< 0.001), suggesting a strong positive relationship between gold and Bitcoin prices. The coefficients for CPI, VIX, and AMD are not statistically significant at conventional levels, although VIX shows a weak negative effect. The residual standard error is 8982, and diagnostic metrics (e.g., VIF values below 2.6) indicate no serious multicollinearity concerns. All predictors show acceptable variance inflation factors, suggesting independent contributions to the model.

```{r}
model <- lm(BITCOIN ~ CPI + VIX + AMD + GOLD , data = df)
summary(model)
car::vif(model)
```

# Stationarity Tests and Variable Transformations

# Null hypothesis study

To determine the stationarity of the variables used in the modeling process, both the Augmented Dickey-Fuller (ADF) test and the Kwiatkowski–Phillips–Schmidt–Shin (KPSS) test were conducted. These tests offer complementary perspectives.The ADF test tests the null hypothesis of a unit root (i.e., non-stationarity). The KPSS test tests the null hypothesis that the series is level-stationary.

Summary of Results: BITCOIN: ADF p-value = 0.3849 → fails to reject the null ⇒ likely non-stationary KPSS = 0.8436 \> 0.739 (1% crit. value) ⇒ rejects stationarity ⇒ confirms non-stationarity

AMD: ADF p-value = 0.8883 ⇒ strong evidence of unit root KPSS = 0.5706 → close to 2.5% critical value ⇒ likely non-stationary

VIX: ADF p-value = 0.6699 ⇒ fails to reject unit root KPSS = 0.2108 → below all critical values ⇒ does not reject stationarity ⇒ inconclusive

GOLD: ADF p-value = 0.656 ⇒ fails to reject unit root KPSS = 0.7904 → above 1% critical value ⇒ rejects stationarity ⇒ likely non-stationary

CPI: ADF p-value = 0.3605 ⇒ non-stationary under ADF KPSS = 0.6555 → between 5% and 1% critical values ⇒ marginal rejection of stationarity

Conclusion: Both tests generally indicate non-stationarity in the level form of all variables, particularly for BITCOIN, AMD, GOLD, and CPI. While the result for VIX is mixed, it is safer to treat all variables as non-stationary and apply differencing before model estimation. These findings justify the use of first differences or transformations in subsequent regression and ARIMA modeling stages.

```{r}
print(adf.test(df$BITCOIN))

print(kpss.test(df$BITCOIN))

print(adf.test(df$AMD))

print(kpss.test(df$AMD))

print(adf.test(df$VIX))

print(kpss.test(df$VIX))

print(adf.test(df$GOLD))

print(kpss.test(df$GOLD))

print(adf.test(df$CPI))

print(kpss.test(df$CPI))
```

# Stationarity Testing – Differenced Variables

Summary of Results: d_BITCOIN: ADF p = 0.3191 (borderline), KPSS = 0.1423 → within 10% critical value ⇒ likely stationary

d_AMD: ADF p = 0.8883, KPSS = 0.2685 → results are weak, but KPSS does not reject stationarity ⇒ inconclusive but acceptable

d_VIX: ADF p = 0.1433, KPSS = 0.1693 → both support stationarity

d_GOLD: ADF p = 0.0849, KPSS = 0.2562 → strong support for stationarity

d_CPI: ADF p = 0.5323, KPSS = 0.3730 → both within bounds, interpreted as stationary

Interpretation: While the ADF test results are not uniformly below conventional significance thresholds, the KPSS test consistently fails to reject the null hypothesis of stationarity. Combined with visual diagnostics and model performance, these findings support the use of the first-differenced variables as stationary inputs for time series regression and forecasting models.

```{r}
df <- df %>%
  mutate(
    d_BITCOIN = c(NA, diff(BITCOIN)),
    d_AMD = c(NA, diff(AMD)),
    d_VIX = c(NA, diff(VIX)),
    d_CPI = c(NA, diff(CPI)),
    d_GOLD = c(NA, diff(GOLD))
  )

adf.test(na.omit(df$d_BITCOIN))
kpss.test(na.omit(df$d_BITCOIN))

adf.test(na.omit(df$AMD))
kpss.test(na.omit(df$d_AMD))

adf.test(na.omit(df$d_VIX))
kpss.test(na.omit(df$d_VIX))

adf.test(na.omit(df$d_GOLD))
kpss.test(na.omit(df$d_GOLD))

adf.test(na.omit(df$d_CPI))
kpss.test(na.omit(df$d_CPI))
```

# Regression Models with Differenced Variables

Model Comparison: Model 1: d_BITCOIN \~ d_AMD Very weak explanatory power (Adj. R² = –0.017), p-value = 0.4385. → d_AMD alone does not explain Bitcoin returns.

Model 2: d_BITCOIN \~ d_AMD + d_VIX Moderate improvement (Adj. R² = 0.1659), VIX is significant (p = 0.028), model nearly significant overall (p = 0.0628). → VIX negatively influences Bitcoin in the short term.

Model 3: Adds d_GOLD No further improvement; GOLD is insignificant (p = 0.6628).

Model 4: Adds d_CPI Slight improvement (Adj. R² = 0.1594), d_VIX remains significant, d_CPI not significant.

Model 5 (best balance): d_BITCOIN \~ d_AMD + d_VIX + d_CPI Highest Adj. R² (0.1947), overall model marginally significant (p = 0.0697), d_VIX remains the only consistently significant predictor.

Interpretation: Across models, d_VIX (volatility index) consistently shows a statistically significant negative relationship with Bitcoin returns, suggesting that rising market fear is associated with short-term drops in Bitcoin price. Other variables, including d_AMD, d_CPI, and d_GOLD, do not show strong or consistent effects. While explanatory power remains modest (Adj. R² \< 0.20), results indicate that investor sentiment plays a more immediate role in Bitcoin fluctuations than macro fundamentals.

```{r}
model_diff <- lm(d_BITCOIN ~ d_AMD, data = df)
summary(model_diff)

model_diff1 <- lm(d_BITCOIN ~ d_AMD + d_VIX, data = df)
summary(model_diff1)

model_diff2 <- lm(d_BITCOIN ~ d_AMD + d_VIX + d_GOLD, data = df)
summary(model_diff2)

model_diff3 <- lm(d_BITCOIN ~ d_AMD + d_VIX + d_GOLD + d_CPI, data = df)
summary(model_diff3)

model_diff4 <- lm(d_BITCOIN ~ d_AMD + d_VIX + d_CPI, data = df)
summary(model_diff4)
```

# Model Diagnostics for Differenced Regression

Durbin-Watson Test: DW = 1.9199, p = 0.426 → No evidence of significant autocorrelation in residuals, which supports model assumptions.

Breusch-Pagan Test: BP = 5.52, p = 0.138 → Homoscedasticity assumption is not rejected; residual variance is approximately constant.

Shapiro-Wilk Test for Normality: W = 0.948, p = 0.2646 → No significant departure from normality; residuals appear normally distributed.

Histogram of Residuals: The residual distribution appears somewhat irregular but loosely centered around zero. While many values fall near the center, the presence of several large deviations on both sides indicates mild dispersion and potential outliers.

Q-Q Plot: The majority of residuals align closely with the theoretical normal line. Slight deviation at tails suggests mild non-normality, though not critical.

```{r}
library(lmtest)
dwtest(model_diff4)
bptest(model_diff4)
shapiro.test(residuals(model_diff4))

hist(residuals(model_diff4), main = "Histogram of Residuals", breaks = 20)
qqnorm(residuals(model_diff4))
qqline(residuals(model_diff4), col = "red")

new_data <- tail(df[, c("d_AMD", "d_VIX", "d_CPI")], 3)
forecast_diff <- predict(model_diff4, newdata = new_data)
forecast_diff
```

# BTC price levels reconstruction

```{r}
last_price <- tail(df$BITCOIN, 1)
forecast_level <- cumsum(c(last_price, forecast_diff))[-1]
```

# 3-Month Forecast of Bitcoin Price

The table and corresponding line chart present a short-term, three-step forecast of Bitcoin price levels based on the final differenced regression model (model_diff4). Forecasted price levels were reconstructed by applying the predicted changes (DBTC) to the last known price observation. Summary of Results: In Month 1, a decline of approximately \$3785.59 is forecasted, reducing the price to \$89,771.61.

In Month 2, the model anticipates a sharp rebound, with an increase of \$14,344.22, lifting the price to \$104,115.83.

In Month 3, the price is expected to stabilize with a slight decline of \$685.36, resulting in a forecasted level of \$103,430.47.

The forecast indicates short-term volatility, with an initial drop followed by a strong upward correction and then a stabilization phase. This pattern reflects the model’s sensitivity to recent changes in d_VIX, d_CPI, and d_AMD, and suggests that sentiment and macro indicators are driving significant short-term fluctuations. However, due to the moderate explanatory power of the model, forecasts should be interpreted cautiously and supplemented with uncertainty estimates (e.g., confidence intervals or simulations).

```{r}
forecast_result <- data.frame(
  Month = c("Month 1", "Month 2", "Month 3"),
  BTC_change = round(forecast_diff, 2),
  BTC_forecast = round(forecast_level, 2)
)

print(forecast_result)

ggplot(forecast_result, aes(x = Month, y = BTC_forecast)) +
  geom_line(group = 1, color = "steelblue") +
  geom_point(size = 3) +
  ggtitle("3-Month Forecast of Bitcoin Price") +
  ylab("Predicted BTC Price (USD)") +
  theme_minimal()
```

# Sensitivity Analysis

A basic sensitivity analysis was performed to evaluate how the regression model responds to typical values of explanatory variables. Mean values of the predictors (d_AMD, d_VIX, d_CPI) were extracted from the dataset and used as input for prediction with the model_diff4. Results: Mean inputs:

d_AMD: 1.98

d_VIX: –0.089

d_CPI: –0.15

Predicted Bitcoin price change (g_base): +3062.27 USD When all predictors are set to their recent average values, the model forecasts a positive short-term change in Bitcoin price of approximately \$3062. This suggests that under typical market conditions observed during the analyzed period, Bitcoin tends to experience modest upward pressure. The result serves as a useful baseline for comparing impacts under alternative scenarios (e.g., in Monte Carlo simulations).

```{r}
X_vars <- names(coef(model_diff4))[-1]
X0 <- df %>% dplyr::select(all_of(X_vars)) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  unlist()
g_base <- predict(model_diff4, newdata = as.data.frame(t(X0)))

X0
g_base
```

# Monte Carlo simulation

```{r}
Sigma <- vcov(model_diff4)
B <- mvrnorm(n = 1000, mu = coef(model_diff4), Sigma = Sigma)

h <- 6  #Forecast horizon
```

# Explanatory variables

```{r}
X_vars <- names(coef(model_diff4))[-1]
X0 <- colMeans(df[, X_vars], na.rm = TRUE)
```

# New predictor matrix for h steps forward

```{r}
X_new_matrix <- matrix(rep(as.numeric(X0), each = h), nrow = h)
X_new_matrix <- cbind(1, X_new_matrix)
```

# Simulation

```{r}
sim_forecasts <- replicate(1000, {
  beta <- mvrnorm(1, coef(model_diff4), Sigma)
  as.numeric(X_new_matrix %*% beta)
})
```

# Statistics summarrizing the forecast

```{r}
forecast_matrix <- t(sim_forecasts)
mean_forecast <- colMeans(forecast_matrix)
lower_95 <- apply(forecast_matrix, 2, quantile, 0.025)
upper_95 <- apply(forecast_matrix, 2, quantile, 0.975)
```

# Recostruction of price levels from differences

```{r}
last_price <- tail(df$BITCOIN, 1)
forecast_price <- cumsum(c(last_price, mean_forecast))[-1]
lower_price <- cumsum(c(last_price, lower_95))[-1]
upper_price <- cumsum(c(last_price, upper_95))[-1]
```

# Monte Carlo Forecast with Fitted Model and Confidence Bands

The chart displays a comprehensive time series of Bitcoin prices, combining historical data, fitted values from the regression model (model_diff4), and a 6-step-ahead Monte Carlo simulation. The simulation incorporates uncertainty in model coefficients to generate a distribution of possible outcomes, represented by a 95% confidence interval (shaded area). Black line: Original Bitcoin price data from 2023 onward. Blue line: In-sample fitted values combined with out-of-sample forecast (levels reconstructed from differenced model). Shaded region: 95% confidence interval based on 1,000 random draws from the multivariate normal distribution of model parameters.

The forecast shows a continued upward trend in Bitcoin price, reflecting recent dynamics captured by explanatory variables (d_AMD, d_VIX, d_CPI). The confidence interval widens with time, illustrating increasing forecast uncertainty typical for multi-step projections. The fitted line tracks the historical data closely, suggesting that the model adequately captures short-run price dynamics, though minor deviations remain. This visualization underscores the probabilistic nature of forecasting and provides a realistic envelope of expected Bitcoin price paths under current market and macroeconomic conditions. The model may be especially useful for short-horizon tactical forecasts, while longer-term projections should be interpreted more cautiously due to widening uncertainty.

```{r}
# Adjusted changes

fitted_diff <- fitted(model_diff4)

# Level reconstruction with fitted

start_price <- df$BITCOIN[1]
fitted_level <- cumsum(c(start_price, fitted_diff))[-1]

# Fitted + forecast

fitted_and_forecast <- ts(c(fitted_level, forecast_price),
                          start = start(btc),
                          frequency = frequency(btc))
                      

# --- Chart ---

# Forecast time

time_start <- time(btc)[length(btc)]
forecast_time <- time_start + seq(1, h) / frequency(btc)

# Chart

plot(btc, 
     xlim = c(start(btc)[1], forecast_time[h]), 
     ylim = range(btc, lower_price, upper_price),
     main = "Monte Carlo Forecast with Fitted Model and Confidence Bands",
     ylab = "Bitcoin Price", 
     xlab = "Time", 
     col = "black", 
     lwd = 1)

# \# Confidence Band

polygon(c(forecast_time, rev(forecast_time)),
        c(upper_price, rev(lower_price)),
        col = rgb(0, 0, 1, 0.2), border = NA)
# Line fitted + forecast

lines(fitted_and_forecast, col = "blue", lwd = 2)

# Actual data

lines(btc, col = "black")

# Legend

legend("topleft", 
       legend = c("Original Data", "Fitted + Forecast", "95% CI (Monte Carlo)"),
       col = c("black", "blue", rgb(0, 0, 1, 0.2)), 
       lwd = c(1, 2, NA), 
       pch = c(NA, NA, 15), 
       pt.cex = 2, 
       bty = "n")
```

# Conclusion

This project developed and compared short-term forecasting models for Bitcoin price using a combination of financial market indicators, macroeconomic variables, and investor sentiment measures. The analysis was conducted on monthly data from 2023 to 2025 and involved a structured sequence of time series modeling, statistical testing, regression analysis, and simulation-based forecasting.

Stationarity analysis revealed that most raw time series were non-stationary, requiring first differencing to achieve stationarity and enable valid statistical inference.

ARIMA modeling provided solid in-sample fit and reasonable short-term forecasts, with moderate error metrics and well-behaved residuals. However, the model assumes past values and internal structure alone determine future prices.

ETS modeling captured trend and level patterns without requiring differencing but performed less accurately in terms of error measures compared to ARIMA.

Multiple linear regression models using differenced explanatory variables (d_AMD, d_VIX, d_CPI, etc.) provided interpretable insights. Among the predictors, d_VIX (volatility index) consistently emerged as a significant factor, negatively associated with short-term Bitcoin price movements.

Model diagnostics (Durbin-Watson, Breusch-Pagan, Shapiro-Wilk) confirmed that the final regression model met key assumptions (no autocorrelation, homoscedasticity, approximate normality).

Monte Carlo simulation was used to quantify forecast uncertainty. The model suggested upward pressure on prices, with a wide but realistic confidence envelope, reflecting growing uncertainty with forecast horizon.

Scenario-based forecasts (e.g., average input sensitivity) confirmed that under typical conditions, the model predicts a modest price increase.

Investing - risk 
No risk - no fun

