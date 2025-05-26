library(minpack.lm)
library(dplyr)
library(xtable)
library(urca)
library(forecast)
library(tseries)

Fert_data <- read.table("C:/Users/leeha/OneDrive/Desktop/Affine mortality/GLG/KORasfrVHbo.txt",header=TRUE,skip=2)
Fert_data <- Fert_data %>% mutate_at(vars(Cohort,Age,ASFR,ASFR1,ASFR2,ASFR3,ASFR4,ASFR5p),as.numeric)
Fert_data<-na.omit(Fert_data)
Fert_data <- subset(Fert_data,Age >= 15 & Age <=49)
Fert_data <- subset(Fert_data, Cohort >= 1960 & Cohort <= 2006)



# Cohort
target_data <- subset(Fert_data,Cohort==1983)

# make model
mod_asfr1 <- nlsLM(ASFR1 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
             data=target_data,start = c(C=0.6,u=30,b=3,lambda=-2),
             lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
coef(mod_asfr1)
mod_asfr2 <- nlsLM(ASFR2 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=target_data,start = c(C=0.6,u=32,b=4,lambda=-1),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
coef(mod_asfr2)
mod_asfr3 <- nlsLM(ASFR3 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=target_data,start = c(C=0.6,u=35,b=3,lambda=-7),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
coef(mod_asfr3)
mod_asfr4 <- nlsLM(ASFR4 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=target_data,start = c(C=0.6,u=35,b=3,lambda=-7),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
coef(mod_asfr4)
mod_asfr5 <- nlsLM(ASFR5p ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=target_data,start = c(C=0.6,u=35,b=3,lambda=-2),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
coef(mod_asfr5)

# 1983 cohort
plot(target_data$Age, target_data$ASFR1, type = "p", pch = 18,cex=1.5,xlim=c(15,50),
     main = "By Birth Order - Cohort Born in 1983", xlab = "Age", ylab = "ASFR")
points(target_data$Age, target_data$ASFR2, type = "p", pch = 22, cex=1.5)
points(target_data$Age, target_data$ASFR3, type = "p", pch = 17, cex=1.5)
points(target_data$Age, target_data$ASFR4, type = "p", pch = 21, cex=1.5)
points(target_data$Age, target_data$ASFR5, type = "p", pch = 20, cex=1.5)
lines(target_data$Age,fitted(mod_asfr1), lwd = 2)
lines(target_data$Age, fitted(mod_asfr2), lwd = 2)
lines(target_data$Age, fitted(mod_asfr3), lwd = 2)
lines(target_data$Age, fitted(mod_asfr4), lwd = 2)
lines(target_data$Age, fitted(mod_asfr5), lwd = 2)
legend("topright", 
       legend = c("1st","2nd","3rd","4th","5th+"),
       pch = c(18,22, 17, 21, 20))
tfr_1983<-mapply(sum,fitted(mod_asfr1),fitted(mod_asfr2),fitted(mod_asfr3),fitted(mod_asfr4),fitted(mod_asfr5))
plot(target_data$Age, target_data$ASFR, type = "p", pch = 15,
     main = "All Birth Order - Cohort Born in 1983", xlab = "Age", ylab = "ASFR")
lines(target_data$Age, tfr_1983, lwd = 2, lty = 1)
legend("topright",legend=c("Observed","Projected"),pch=c(15,NA),lty=c(NA,1))

# Period ASFR forecast

period_Fert_data <- read.table("C:/Users/leeha/OneDrive/Desktop/Affine mortality/GLG/KORasfrRRbo.txt",header=TRUE,skip=2)
head(period_Fert_data)
period_Fert_data <- period_Fert_data %>% mutate_at(vars(ASFR1,ASFR2,ASFR3,ASFR4,ASFR5p),~ replace(.,.==".",NA))
period_Fert_data <- period_Fert_data %>% mutate_at(vars(Age,ASFR1,ASFR2,ASFR3,ASFR4,ASFR5p),as.numeric)
period_Fert_data <- period_Fert_data[complete.cases(Fert_data),]

period_Fert_data <- subset(period_Fert_data,Age>= 15 & Age<=49)


#target
prd_target_data <- subset(period_Fert_data,Year==2000)
prd_target_data_1 <- subset(period_Fert_data,Year==2023)
# make model
prd_mod_asfr1 <- nlsLM(ASFR1 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=prd_target_data,start = c(C=0.6,u=30,b=3,lambda=-2),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
prd_mod_asfr1_1 <- nlsLM(ASFR1 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                       data=prd_target_data_1,start = c(C=0.6,u=30,b=3,lambda=-2),
                       lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)

coef(prd_mod_asfr1)
prd_mod_asfr2 <- nlsLM(ASFR2 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=prd_target_data,start = c(C=0.6,u=32,b=3,lambda=-2),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
prd_mod_asfr2_1 <- nlsLM(ASFR2 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                       data=prd_target_data_1,start = c(C=0.6,u=32,b=3,lambda=-1),
                       lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)

coef(prd_mod_asfr2)
prd_mod_asfr3 <- nlsLM(ASFR3 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=prd_target_data,start = c(C=0.6,u=35,b=3,lambda=-7),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
prd_mod_asfr3_1 <- nlsLM(ASFR3 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                       data=prd_target_data_1,start = c(C=0.6,u=35,b=5,lambda=-5),
                       lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)

coef(prd_mod_asfr3)
prd_mod_asfr4 <- nlsLM(ASFR4 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=prd_target_data,start = c(C=0.6,u=35,b=3,lambda=-7),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
prd_mod_asfr4_1 <- nlsLM(ASFR4 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                       data=prd_target_data_1,start = c(C=0.6,u=35,b=3,lambda=-7),
                       lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)

coef(prd_mod_asfr4)
prd_mod_asfr5 <- nlsLM(ASFR5p ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                   data=prd_target_data,start = c(C=0.6,u=35,b=3,lambda=-7),
                   lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
prd_mod_asfr5_1 <- nlsLM(ASFR5p ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                       data=prd_target_data_1,start = c(C=0.6,u=35,b=3,lambda=-2),
                       lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)

coef(prd_mod_asfr5)

#plot
plot(prd_target_data$Age, prd_target_data$ASFR1, type = "p", pch = 18,
     main = "By Birth Order - Period data observed in 2000 and 2023", xlab = "Age", ylab = "ASFR")
points(prd_target_data$Age, prd_target_data$ASFR2, type = "p", pch = 22)
points(prd_target_data$Age, prd_target_data$ASFR3, type = "p", pch = 17)
points(prd_target_data$Age, prd_target_data$ASFR4, type = "p", pch = 21)
points(prd_target_data$Age, prd_target_data$ASFR5, type = "p", pch = 20)
points(prd_target_data$Age, prd_target_data_1$ASFR1, type = "p", pch = 18,col="blue")
points(prd_target_data$Age, prd_target_data_1$ASFR2, type = "p", pch = 22,col="blue")
points(prd_target_data$Age, prd_target_data_1$ASFR3, type = "p", pch = 17,col="blue")
points(prd_target_data$Age, prd_target_data_1$ASFR4, type = "p", pch = 21,col="blue")
points(prd_target_data$Age, prd_target_data_1$ASFR5, type = "p", pch = 20,col="blue")
lines(prd_target_data$Age, fitted(prd_mod_asfr1), lwd = 2)
lines(prd_target_data$Age, fitted(prd_mod_asfr2), lwd = 2)
lines(prd_target_data$Age, fitted(prd_mod_asfr3), lwd = 2)
lines(prd_target_data$Age, fitted(prd_mod_asfr4), lwd = 2)
lines(prd_target_data$Age, fitted(prd_mod_asfr5), lwd = 2)
lines(prd_target_data$Age, fitted(prd_mod_asfr1_1), lwd = 2,col="blue")
lines(prd_target_data$Age, fitted(prd_mod_asfr2_1), lwd = 2,col="blue")
lines(prd_target_data$Age, fitted(prd_mod_asfr3_1), lwd = 2,col="blue")
lines(prd_target_data$Age, fitted(prd_mod_asfr4_1), lwd = 2,col="blue")
lines(prd_target_data$Age, fitted(prd_mod_asfr5_1), lwd = 2,col="blue")
legend("topright", 
       legend = c("1st","2nd","3rd","4th","5th+","2000","2023"),
       pch = c(18,22, 17, 21, 20,NA,NA),lty = c(NA,NA,NA,NA,NA,1,1),col = c("black","black","black","black","black","black","blue"))

# total
prd_tfr<-mapply(sum,fitted(prd_mod_asfr1),fitted(prd_mod_asfr2),fitted(prd_mod_asfr3),fitted(prd_mod_asfr4),fitted(prd_mod_asfr5))
prd_tfr_1<-mapply(sum,fitted(prd_mod_asfr1_1),fitted(prd_mod_asfr2_1),fitted(prd_mod_asfr3_1),fitted(prd_mod_asfr4_1),fitted(prd_mod_asfr5_1))

plot(prd_target_data$Age, prd_target_data$ASFR, type = "p", pch = 15,cex=0.5,
     main = "All Birth Order - Period data observed in 2000 and 2023", xlab = "Age", ylab = "ASFR")
points(prd_target_data$Age, prd_target_data_1$ASFR, type = "p", pch = 15,col="blue",cex=0.5)

lines(prd_target_data$Age, prd_tfr, lwd = 2)
lines(prd_target_data$Age, prd_tfr_1, lwd = 2, col="blue")
legend("topright",legend=c("Observed-2000","Projected-2000","Observed-2023","Projected-2023"),pch=c(15,NA,15,NA),lty=c(NA,1,NA,1),col=c("black","black","blue","blue"),cex=0.8)
plot(prd_target_data$Age, prd_target_data_1$ASFR, type = "p", pch = 15,
     main = "All Birth Order - Period data observed in 2023", xlab = "Age", ylab = "ASFR")
lines(prd_target_data$Age, prd_tfr_1, lwd = 2, col="blue")

# make function
params_fn <- function(year,initial_paras){
    #target
    prd_target_data <- subset(period_Fert_data,Year==year)
    
    # make model
    prd_mod_asfr1 <- nlsLM(ASFR1 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                           data=prd_target_data,start = initial_paras[[1]],
                           lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
    
    prd_mod_asfr2 <- nlsLM(ASFR2 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                           data=prd_target_data,start = initial_paras[[2]],
                           lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
    
    prd_mod_asfr3 <- nlsLM(ASFR3 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                           data=prd_target_data,start = initial_paras[[3]],
                           lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
    
    prd_mod_asfr4 <- nlsLM(ASFR4 ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                           data=prd_target_data,start = initial_paras[[4]],
                           lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)
    
    prd_mod_asfr5 <- nlsLM(ASFR5p ~ C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b)),
                           data=prd_target_data,start = initial_paras[[5]],
                           lower = c(-Inf,-Inf,0,-Inf),upper=c(Inf,Inf,Inf,0),trace = TRUE)

    coef(prd_mod_asfr5)
    param_list <- list(coef(prd_mod_asfr1),coef(prd_mod_asfr2),coef(prd_mod_asfr3), coef(prd_mod_asfr4),coef(prd_mod_asfr5))

    return(param_list)
}

initial_paras <- list(c(C=0.6,u=30,b=3,lambda=-2),c(C=0.6,u=32,b=3,lambda=-2),c(C=0.6,u=35,b=3,lambda=-7),c(C=0.6,u=35,b=3,lambda=-7),c(C=0.6,u=35,b=3,lambda=-7))
params_period_list <- list()
for (year in 2000:2023){
    print(year)
    temp_list <- params_fn(year,initial_paras)
    params_period_list <- append(params_period_list, temp_list)
    initial_paras  <- temp_list
    
}

# data frame
years <- rep(2000:2023, each = 5)   # 각 년도별로 5개씩 반복 (총 120)
birthorders <- rep(c("ASFR1", "ASFR2", "ASFR3", "ASFR4", "ASFR5"), times = 24)  # 24회 반복

params_matrix <- do.call(rbind, params_period_list)
params_df <- data.frame(Year = years,
                        Birthorder = birthorders,
                        params_matrix,
                        row.names = NULL)
print(params_df)

# ADF test

params <- c("C","u","b","lambda")
birthorder <- unique(birthorders)
adf_table <- data.frame(Birthorder=character(),Parameter=character(),Statistic = numeric(),Critical_Value=numeric(),Decision=character())

for(bo in birthorder){
    subset_data <- subset(params_df,Birthorder == bo)
    for (param in params){
        ts_data <- subset_data[[param]]
        #print(ts_data)
        adf_test <- ur.df(ts_data, type="drift",selectlags = "AIC")
        #print(adf_test)
        test_stat <- adf_test@teststat[,"tau2"]
        print(adf_test@cval)
        critical_value <- adf_test@cval["tau2","5pct"]
        decision <- ifelse(test_stat < critical_value,"Stationary","Random Walk")
        
        adf_table<-rbind(adf_table,data.frame(Birthorder=bo,Parameter=param,Statistic = round(test_stat,4),Critical_Value=round(critical_value,4),Decision=decision))
    }
}
adf_table

print(xtable(adf_table),include.rownames=FALSE)

# Estimating 2024's parameter

 # forecast

forecast_AR1 <- function(parameter){
    ts_data <- ts(parameter,start=2000,frequency = 1)
    
    fit <- auto.arima(ts_data)
    forecast <- forecast(fit,h=1)
    se<-(forecast$upper[2]-forecast$lower[2])/(2*1.96)
    
    return(c(forecast$mean[1],se))
}

param_2024 <-function(Order){
    target_data <- subset(params_df, Birthorder == Order)
    target_data <- target_data[order(target_data$Year), ]
    pred_C      <- forecast_AR1(target_data$C)
    pred_u      <- forecast_AR1(target_data$u)
    pred_b      <- forecast_AR1(target_data$b)
    pred_lambda <- forecast_AR1(target_data$lambda)
    print(pred_C)
    pred <- rbind(pred_C,pred_u,pred_b,pred_lambda)
    return(pred)
    
} 
# test

# not negative values
truncated_rnorm <- function(n, mean, sd, lower_bound = 0) {
    values <- rnorm(n, mean, sd)
    while (any(values < lower_bound)) {
        values[values < lower_bound] <- rnorm(sum(values < lower_bound), mean, sd)
    }
    return(values)
}

asfr_birth <- function(params_df,Order,Years, Iteration){
    if (Order == "ASFR1"){
        target_data <- subset(params_df, Birthorder==Order )
        target_data <- target_data[order(target_data$Year), ]    }
    if (Order == "ASFR2"){
        target_data <- subset(params_df, Birthorder==Order )
        target_data <- target_data[order(target_data$Year), ]    }
    if (Order == "ASFR3"){
        target_data <- subset(params_df, Birthorder==Order )
        target_data <- target_data[order(target_data$Year), ]    }
    if (Order == "ASFR4"){
        target_data <- subset(params_df, Birthorder==Order )
        target_data <- target_data[order(target_data$Year), ]    }
    if (Order == "ASFR5"){
        target_data <- subset(params_df, Birthorder==Order )
        target_data <- target_data[order(target_data$Year), ]    }
    
    C_ts_data <- ts(target_data$C,start=2000,frequency = 1)
    C_fit <- arima(C_ts_data,order=c(0,1,0))
    C_sigma <- sqrt(C_fit$sigma2)
    generated_C <- matrix(nrow=Years,ncol=Iteration)
    
    for (sim in 1:Iteration){
        C_t <- numeric(Years+1)
        C_t[1]<-tail(target_data$C,1)
        for (t in 2:(Years+1)){
            C_epsilon <- rnorm(1,mean=0,sd=C_sigma)
            C_t[t] <- C_t[t-1] + C_epsilon
        }
        generated_C[,sim] <- C_t[2:length(C_t)]
    }
    
#u
    u_ts_data <- ts(target_data$u,start=2000,frequency = 1)
    u_fit <- arima(u_ts_data,order=c(0,1,0))
    u_sigma <- sqrt(u_fit$sigma2)
    
    generated_u <-matrix(nrow=Years,ncol=Iteration)
    for (sim in 1:Iteration){
        u_t <- numeric(Years+1)
        u_t[1]<-tail(target_data$u,1)
        
        for (t in 2:(Years+1)){
            u_epsilon <- rnorm(1,mean=0,sd=u_sigma)
            u_t[t] <- u_t[t-1]+u_epsilon
        }
        generated_u[,sim] <- u_t[2:length(u_t)]
    }
    
#b
    b_ts_data <- ts(target_data$b,start=2000,frequency = 1)
    b_fit <- arima(b_ts_data,order=c(0,1,0))
    b_sigma <- sqrt(b_fit$sigma2)
    generated_b <- matrix(nrow=Years,ncol=Iteration)
    
    for (sim in 1:Iteration){
        b_t <- numeric(Years+1)
        b_t[1]<-tail(target_data$b,1)
        for (t in 2:(Years+1)){
            b_epsilon <- rnorm(1,mean=0,sd=b_sigma)
            b_t[t] <- b_t[t-1] + b_epsilon
        }
        generated_b[,sim] <- b_t[2:length(b_t)]
    }

#l
    l_ts_data <- ts(target_data$lambda,start=2000,frequency = 1)
    adjusted_l_ts_data <- log(abs(l_ts_data))
    l_fit <- arima(l_ts_data,order=c(0,1,0))
    l_sigma <- sqrt(l_fit$sigma2)
    generated_l <- matrix(nrow=Years,ncol=Iteration)
    
    for (sim in 1:Iteration){
        l_t <- numeric(Years+1)
        Adjusted_values <- numeric((Years+1))
        Adjusted_values[1]<-tail(adjusted_l_ts_data,1)
        
       
        l_t[1]<-tail(target_data$lambda,1)
        
        for (t in 2:(Years+1)){
            l_epsilon <- rnorm(1,mean=0,sd=l_sigma)
            Adjusted_values[t] <- Adjusted_values[t-1] + l_epsilon
            l_t[t] <- -exp(Adjusted_values[t])
        }
        generated_l[,sim] <- l_t[2:length(l_t)]
    }

    # ASFR
    asfr_fnc <- function(C,u,b,lambda){
        Age <- 15:49
        asfr_list<-lapply(Age,function(Age){
            C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age+1-u)/b))-C *pgamma(lambda^{-2},lambda^{-2} * exp(lambda * (Age-u)/b))})
        return(asfr_list)
    }
    # first year
    Age <- 15:49
    
    Itr_list <-list()

    for ( i in 1:Iteration){
        asfr_mat<-matrix(ncol=Years,nrow=length(Age))
        for ( year in 1:Years){
            asfr_mat[,year]<-unlist(asfr_fnc(generated_C[year,i],generated_u[year,i],generated_b[year,i],generated_l[year,i]))
        }
        Itr_list[[i]]<-asfr_mat
    }

    return(Itr_list)    
}

tfr_function <- function(params_df,Years,Iteration){
    tfr_list <- list()
    for (i in 1:Iteration){
        tfr_list[[i]]<-asfr_birth(params_df,"ASFR1",Years, Iteration)[[i]]+asfr_birth(params_df,"ASFR2",Years, Iteration)[[i]]+asfr_birth(params_df,"ASFR3",Years, Iteration)[[i]]+asfr_birth(params_df,"ASFR4",Years, Iteration)[[i]]+asfr_birth(params_df,"ASFR5",Years, Iteration)[[i]]
        print(i)
        }
    return(tfr_list)
}

TFR_list<-tfr_function(params_df,Years=100,Iteration = 200)

