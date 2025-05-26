library(readxl)
cpi_data <- read_excel("C:/Users/leeha/OneDrive/Desktop/Affine mortality/wage/cpi.xlsx")
head(cpi_data)
cpi_year <- as.numeric(cpi_data$시점)
cpi_rate <- cpi_data$총지수 / 100
# from 1992
inflation <- ts(cpi_rate[27:length(cpi_rate)],start=cpi_year[27])
inflation_model <- auto.arima(inflation)
summary(inflation_model)
plot(inflation)
theta <- inflation_model$coef[[1]]

Years<-100
Iteration<-100
generated_cpi <- matrix(nrow=Iteration,ncol=Years)
#Cpi
cpi_sigma <- sqrt(inflation_model$sigma2)
for (sim in 1:Iteration){
    cpi_t <- numeric(Years+1)
    cpi_t[1]<-tail(inflation,1)
    cpi_epsilon <- rnorm(Years+1,mean=0,sd=cpi_sigma)
    for (t in 2:(Years+1)){
        
        cpi_t[t] <- cpi_t[t-1] + cpi_epsilon[t] + theta*cpi_epsilon[t-1]
    }
    generated_cpi[sim,] <- cpi_t[2:length(cpi_t)]
}
