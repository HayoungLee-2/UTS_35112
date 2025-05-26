library(dplyr)
library(fmsb)
library(pyramid)
library(scales)
library(stringr)
Pop_Fund <- function(Forecast_year=100,Iteration=100,Premium_rate=0.09,Replacement_rate=0.3,option){

    
    wx <- 101
    tp <- fp <- matrix(rep(0, wx*(Forecast_year+1)), wx, (Forecast_year+1))
    Kpopl <- read.table("C:/Users/leeha/OneDrive/Desktop/Affine mortality/Affine mortality/Population.txt",header=TRUE,skip=2)
    
    Kpopl <- Kpopl %>% filter(str_detect(Age,"^[0-9]+$"))%>% mutate_at(vars(Year,Age,Female,Male,Total),as.double)%>%filter(Age >= 0 & Age <= 100)
    # initial year is 2023
    tp[,1] <- Kpopl[Kpopl$Year==2023,][,"Total"]
    fp[,1] <- Kpopl[Kpopl$Year==2023,][,"Female"]
    
    # population matrix
    pop_data <- matrix(nrow = (Forecast_year+1), ncol = Iteration)
    # Dependency matrix
    dependecy_data <- matrix(nrow = (Forecast_year+1), ncol = Iteration)
    # Fund matrix
    fund_data <- matrix(nrow = (Forecast_year+1), ncol = Iteration)
    # Payer change matrix
    payer_forecast_mat <- matrix(nrow = Iteration, ncol= Forecast_year-1)
    lifeexp_mat<-matrix(nrow=Iteration,ncol=(Forecast_year-1))
    Result_mat<-matrix(nrow=Iteration,ncol=(Forecast_year-1))
    generated_cpi <- matrix(nrow=Iteration,ncol=Forecast_year)
    #Cpi
    
    cpi_sigma <- sqrt(inflation_model$sigma2)
    for (sim in 1:Iteration){
        cpi_t <- numeric(Forecast_year+1)
        cpi_t[1]<-tail(inflation,1)
        cpi_epsilon <- rnorm(Forecast_year+1,mean=0,sd=cpi_sigma)
        for (t in 2:(Forecast_year+1)){
            
            cpi_t[t] <- cpi_t[t-1] + cpi_epsilon[t] + theta*cpi_epsilon[t-1]
        }
        generated_cpi[sim,] <- cpi_t[2:length(cpi_t)]
    }
    for (itr in 1:Iteration){
        
        for (j in 2:(Forecast_year+1)) {
            asfr <- c(rep(0, 15), TFR_list[[itr]][,j-1], rep(0, wx-50))
            
            baby <- sum(fp[, j-1]*asfr)
            babym <- as.integer(baby*1.06/2.06+0.5)
            babyf <- as.integer(baby*1/2.06+0.5)
            tp[1, j] <- as.integer(baby)
            fp[1, j] <- babyf
            
            sfx <- Generated_fm[[itr]][1:(wx-1),j-1]
            stx <- Generated_t[[itr]][1:(wx-1),j-1]
            
            tp[2:wx, j] <- as.integer(tp[1:(wx-1), j-1]*stx)
            fp[2:wx, j] <- as.integer(fp[1:(wx-1), j-1]*sfx)
        }
        
        
        pop_data[, itr] <- colSums(tp)
        
        Dependcy_Ratio <- colSums(tp[66:101,])/colSums(tp[16:65,])*100
        dependecy_data[, itr] <- Dependcy_Ratio
        # pension pay rate 73.3%
        # pension take rate 51.2%
        
        # Fund
        Fund_forecast <- c()
        Fund_forecast <- append(Fund_forecast,fund_2023)
        for ( i in 1:Forecast_year){
            C <-12*Premium_rate*wage_years[i]*sum(tp[(18+1):(59+1),i+1])*73.9/100
            E <-12*Replacement_rate*wage_years[i]*sum(tp[(65+1):(100+1),i+1])*51.2/100
            Fund_forecast<-append(Fund_forecast,(fund_2023+C-E)*(1+avg_investment/100))
        }
        fund_data[, itr] <- Fund_forecast
        
        if ( option =="Macroslide"){
            # Apply macro slide
            Itr_payer<-sapply(1:Forecast_year, function(i) sum(tp[(18+1):(59+1),i+1]))
            payer_forecast_mat[itr,]<-diff(Itr_payer)/head(Itr_payer,-1)
            # life expectancy increase rate
            
            Itr_lifeexp<- sapply(1:Forecast_year, function(j) sum(Generated_tm[[itr]][,j]))
            lifeexp_mat[itr,]<-diff(Itr_lifeexp)/head(Itr_lifeexp,-1)
            #print(length(payer_forecast_mat[itr,]))
            #print(length(lifeexp_mat[itr,]))
            
            Adjustment_list <-payer_forecast_mat[itr,]+lifeexp_mat[itr,]
            #print(length(Adjustment_list))
            Result_mat[itr,]<-unlist(mapply(function(a,b){
                if ((a>0)&(a>b)){
                    (a-b)}
                else if (a<0){
                    a
                }
                else{
                    0
                }
            },generated_cpi[itr,2:length(generated_cpi[itr,])],Adjustment_list))
            
            
            
            
            # Fund
            Fund_forecast <- c()
            Fund_forecast <- append(Fund_forecast,fund_2023)
            for ( i in 1:Forecast_year){
                if (i==1){
                    C <-12*Premium_rate*wage_years[i]*sum(tp[(18+1):(59+1),i+1])*73.3/100
                    E <-12*Replacement_rate*wage_years[i]*sum(tp[(65+1):(100+1),i+1])*51.2/100
                    Fund_forecast<-append(Fund_forecast,(fund_2023+C-E)*(1+avg_investment/100))
                    
                }
                else{
                    C <-12*Premium_rate*wage_years[i]*sum(tp[(18+1):(59+1),i+1])*73.3/100
                    E <-(1+Result_mat[itr,i-1])*12*Replacement_rate*wage_years[i]*sum(tp[(65+1):(100+1),i+1])*51.2/100
                    Fund_forecast<-append(Fund_forecast,(fund_2023+C-E)*(1+avg_investment/100))
                    
                }
                            }
            fund_data[, itr] <- Fund_forecast
            #print(fund_data)
            
        }
                
        
        }
    
    # Df change
    Years_list <- 2023:(2023+Forecast_year)
    pop_df <- as.data.frame(t(pop_data))
    colnames(pop_df) <- as.character(Years_list)
    pop_df$iteration <- paste0("sim", 1:Iteration)
    
    #transform
    df_long <- pivot_longer(pop_df, cols = -iteration, names_to = "year", values_to = "population")
    df_long$year <- as.integer(df_long$year)
    
        
    
    # summary
    summary_df <- df_long %>%
        group_by(year) %>%
        summarise(
            median = median(population),
            lower = quantile(population, 0.025),
            upper = quantile(population, 0.975)
        )
        
    
    if (option=="Population"){
        # Pop
        
        return(ggplot() +
            
            geom_line(data = df_long, aes(x = year, y = population, group = iteration), 
                      color = "grey80", alpha = 0.5) +
            
            
            geom_ribbon(data = summary_df, aes(x = year, ymin = lower, ymax = upper), 
                        fill = "lightblue", alpha = 0.4) +
            
            
            geom_line(data = summary_df, aes(x = year, y = median), 
                      color = "blue", size = 1.2) +
            
            scale_y_continuous(labels = comma) +
            labs(
                title = "Population Forecast with 95% Confidence Interval",
                x = "Year",
                y = "Population"
            ) +
            theme_minimal() )
    }
      
    
    if (option=="Dependency Ratio"){
        dependency_df <- as.data.frame(t(dependecy_data))
        colnames(dependency_df) <- as.character(Years_list)
        dependency_df$iteration <- paste0("sim", 1:Iteration)
        
        df_long_dep <- pivot_longer(dependency_df, cols = -iteration, names_to = "year", values_to = "dependencyratio")
        df_long_dep$year <- as.integer(df_long_dep$year)
        
        summary_df_dep <- df_long_dep %>%
            group_by(year) %>%
            summarise(
                median = median(dependencyratio,na.rm = TRUE),
                lower = quantile(dependencyratio, 0.025,na.rm = TRUE),
                upper = quantile(dependencyratio, 0.975,na.rm = TRUE)
            )
        # Dep
        return(ggplot() +
            
            geom_line(data = df_long_dep, aes(x = year, y = dependencyratio, group = iteration), 
                      color = "grey80", alpha = 0.5) +
            
            
            geom_ribbon(data = summary_df_dep, aes(x = year, ymin = lower, ymax = upper), 
                        fill = "lightblue", alpha = 0.4) +
            
            
            geom_line(data = summary_df_dep, aes(x = year, y = median), 
                      color = "blue", size = 1.2) +
            
            scale_y_continuous(labels = comma) +
            labs(
                title = "Dependecy Ratio Forecast with 95% Confidence Interval",
                x = "Year",
                y = "Dependecy Ratio"
            ) +
            theme_minimal() )
        
    }

    
    if (option == "Fund" | option == "Macroslide"){
        fund_df <- as.data.frame(t(fund_data))
        colnames(fund_df) <- as.character(Years_list)
        fund_df$iteration <- paste0("sim", 1:Iteration)
        df_long_fund <- pivot_longer(fund_df, cols = -iteration, names_to = "year", values_to = "Fund")
        df_long_fund$year <- as.integer(df_long_fund$year)
        #print(df_long_fund,n=8000)
        summary_df_fund <- df_long_fund %>%
            group_by(year) %>%
            summarise(
                median = median(Fund),
                lower = quantile(Fund, 0.025),
                upper = quantile(Fund, 0.975)
            )
        #print(fund_df[,-ncol(fund_df)])
        zero_years <- apply(fund_df[,-ncol(fund_df)],1,function(x){Years_list[min(which(x <= 0))]})
        #print(zero_years)
        max_years <- apply(fund_df[,-ncol(fund_df)],1,function(x){Years_list[which.max(x)]})
        #Fund
        
        return(ggplot() +
            
            geom_line(data = df_long_fund, aes(x = year, y = Fund, group = iteration), 
                      color = "grey80", alpha = 0.5) +
            
            geom_vline(xintercept = median(zero_years,na.rm=TRUE),color = "orange", linetype = "dashed", )+
            annotate("rect",xmin=quantile(zero_years,0.025,na.rm=TRUE),xmax=quantile(zero_years,0.975,na.rm=TRUE),ymin=-Inf,ymax=Inf,fill="lightpink",alpha=0.3)+
            #geom_vline(xintercept = quantile(zero_years,0.025), color = "red",size = 0.8 ,linetype = "dotted") +
            #geom_vline(xintercept = quantile(zero_years,0.975), color = "red",size = 0.8, linetype = "dotted") +
            geom_hline(yintercept = 0,color="orange", linetype = "dashed", size = 1)+
            
            geom_vline(xintercept = median(max_years),color = "darkgreen", linetype = "dashed", )+
            
            geom_ribbon(data = summary_df_fund, aes(x = year, ymin = lower, ymax = upper), 
                        fill = "lightblue", alpha = 0.4) +
            
            scale_x_continuous(limits = c(2023, 2023+Forecast_year)) +
            geom_line(data = summary_df_fund, aes(x = year, y = median), 
                      color = "blue", size = 1.2) +
            
            scale_y_continuous(labels = comma) +
            labs(
                title = "Fund Forecast with 95% Confidence Interval",
                subtitle = paste0("Median depletion year = ", median(zero_years,na.rm=TRUE)," Median Maximum fund year = ",median(max_years)),
                x = "Year",
                y = "Fund Balance"
            ) +
            theme_minimal() )
        
    }
 
}
Pop_Fund(option="Population",Forecast_year = 50)
Pop_Fund(option="Dependency Ratio",Forecast_year = 50)
Pop_Fund(option="Fund")
# Macro slide
Pop_Fund(option="Macroslide")
# Premium rate increase
Pop_Fund(Forecast_year = 100,Premium_rate=0.13,option="Fund")










