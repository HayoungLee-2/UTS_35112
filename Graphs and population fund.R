library(ggplot2)

Forecast_year=80
Iteration=200
Premium_rate=0.09
Start_age = 65
Replacement_rate=0.4
avg_investment = 4
fund_2023 <- 950.344 * 10^12
avg_growth<-3.7


Simu_fund <- function(Premium_rate=0.09,Replacement_rate=0.4 ,Start_age=65, avg_growth=3.7, avg_investment=4,OPTION,detail=FALSE){
    fund_2023 <- 950.344 * 10^12
    wage_years <-c()
    for (yrs in 1:Forecast_year){
        wage_years<-append(wage_years,wage_2023 * (1+avg_growth/100)^yrs)
    }
    
    Kpopl_raw <- read.table("C:/Users/leeha/OneDrive/Desktop/Affine mortality/Affine mortality/Population.txt",header=TRUE,skip=2)
    wx <- 110
    tp <- fp <- matrix(rep(0, wx*(Forecast_year+1)), wx, (Forecast_year+1))
    
    Kpopl <- Kpopl_raw %>% filter(str_detect(Age,"^[0-9]+$"))%>% mutate_at(vars(Year,Age,Female,Male,Total),as.double)%>%filter(Age >= 0 & Age <= 109)
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
            babym <- baby*1.06/2.06+0.5
            babyf <- baby*1/2.06+0.5
            tp[1, j] <- baby
            fp[1, j] <- babyf
            
            stx <- Generated_t[[itr]][1:(wx-1),j-1]
            
            tp[2:wx, j] <- tp[1:(wx-1), j-1]*stx
            fp[2:wx, j] <- fp[1:(wx-1), j-1]*stx
        }
        pop_data[, itr] <- colSums(tp)
        
        Dependcy_Ratio <- colSums(tp[66:110,])/colSums(tp[16:65,])*100
        dependecy_data[, itr] <- Dependcy_Ratio
        # pension pay rate 73.3%
        # pension take rate 51.2%
        
        if (OPTION=="fund")
        {
            Fund_forecast <- numeric(Forecast_year+1)
            Fund_forecast[1] <- fund_2023
            for ( i in 1:Forecast_year){
                C <-12*Premium_rate*wage_years[i]*sum(tp[(18+1):(59+1),i+1])*73.9/100
                E <-12*Replacement_rate*wage_years[i]*sum(tp[(Start_age+1):(109+1),i+1])*51.2/100
                Fund_forecast[i+1]<-(Fund_forecast[i]+C - E)*(1+avg_investment/100)
            }
            fund_data[, itr] <- Fund_forecast
            
        }
                
        if ( OPTION =="Macroslide"){
            # Apply macro slide
            Itr_payer<-sapply(1:Forecast_year, function(i) sum(tp[(18+1):(59+1),i+1]))
            payer_forecast_mat[itr,]<-diff(Itr_payer)/head(Itr_payer,-1)
            # life expectancy increase rate
            
            Itr_lifeexp<- sapply(1:Forecast_year, function(j) sum(Generated_t[[itr]][,j]))
            lifeexp_mat[itr,]<-diff(Itr_lifeexp)/head(Itr_lifeexp,-1)
            
            
            Adjustment_list <- -payer_forecast_mat[itr,]+lifeexp_mat[itr,]
            
            Result_mat[itr,]<-pmax(0,generated_cpi[itr,2:length(generated_cpi[itr,])]-Adjustment_list)
            
            # Fund
            Fund_forecast <- numeric(Forecast_year+1)
            Fund_forecast[1] <- fund_2023
            for ( i in 1:Forecast_year){
                if (i==1){
                    C <-12*Premium_rate*wage_years[i]*sum(tp[(18+1):(59+1),i+1])*73.3/100
                    E <-12*Replacement_rate*wage_years[i]*sum(tp[(Start_age+1):(109+1),i+1])*51.2/100
                    Fund_forecast[i + 1] <- (Fund_forecast[i] + C - E) * (1 + avg_investment / 100)
                }
                else{
                    C <-12*Premium_rate*wage_years[i]*sum(tp[(18+1):(59+1),i+1])*73.3/100
                    E <-12*(1-Result_mat[itr,i-1])*Replacement_rate*wage_years[i]*sum(tp[(Start_age+1):(109+1),i+1])*51.2/100
                    Fund_forecast[i + 1] <- (Fund_forecast[i] + C - E) * (1 + avg_investment / 100)     
                }
            }
            fund_data[, itr] <- Fund_forecast
            
        }
        
    }
    if (OPTION == "Macroslide" & detail==TRUE) {
        return(list(
            fund_data = fund_data,
            Result_mat = Result_mat,
            generated_cpi = generated_cpi,
            payer_forecast_mat = payer_forecast_mat,
            lifeexp_mat = lifeexp_mat
        ))
    }
    if (OPTION == "fund" & detail==TRUE) {
        return(fund_data)
    }
    if (OPTION=="pop"){
        Years_list <- 2023:(2023+Forecast_year)
        pop_df <- as.data.frame(t(pop_data))
        colnames(pop_df) <- as.character(Years_list)
        pop_df$iteration <- paste0("sim", 1:Iteration)
        df_long <- pivot_longer(pop_df, cols = -iteration, names_to = "year", values_to = "population")
        df_long$year <- as.integer(df_long$year)
        summary_df <- df_long %>%
            group_by(year) %>%
            summarise(
                median = median(population),
                lower = quantile(population, 0.025),
                upper = quantile(population, 0.975)
            )
        return(ggplot()+
            geom_line(data = df_long, aes(x = year, y = population, group = iteration),
                      color = "grey80", alpha = 0.5)+
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
            theme_minimal()
        )
        summary_df_2070 <- summary_df %>% filter(year == 2070)
        return(summary_df_2070)
    }
    if (OPTION=="dep"){
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
        df_long_dep_50 <- df_long_dep %>% filter(year <= min(year) + 49)
        summary_df_dep_50 <- summary_df_dep %>% filter(year<=min(year) + 49)
        
        summary_df_dep_2070 <- summary_df_dep %>% filter(year == 2070)
        return(summary_df_dep_2070)
        Dep
        return(ggplot() +

            geom_line(data = df_long_dep_50, aes(x = year, y = dependencyratio, group = iteration),
                      color = "grey80", alpha = 0.5) +


            geom_ribbon(data = summary_df_dep_50, aes(x = year, ymin = lower, ymax = upper),
                        fill = "lightblue", alpha = 0.4) +


            geom_line(data = summary_df_dep_50, aes(x = year, y = median),
                      color = "blue", size = 1.2) +

                scale_y_continuous(
                    labels = scales::comma,
                    breaks = seq(10, 150, by = 20)  # ← 여기 추가!
                )  +
            labs(
                title = "Old-age Dependecy Ratio Forecast with 95% Confidence Interval",
                x = "Year",
                y = "Ratio"
            ) +
            theme_minimal()
        )
    }
    if(OPTION=="fund"|OPTION=="Macroslide"){
        fund_df <- as.data.frame(t(fund_data))
        colnames(fund_df) <- as.character(Years_list)
        fund_df$iteration <- paste0("sim", 1:Iteration)
        df_long_fund <- pivot_longer(fund_df, cols = -iteration, names_to = "year", values_to = "Fund")
        df_long_fund$year <- as.integer(df_long_fund$year)
        summary_df_fund <- df_long_fund %>%
            group_by(year) %>%
            summarise(
                median = median(Fund),
                lower = quantile(Fund, 0.025),
                upper = quantile(Fund, 0.975)
            )
        zero_years <- apply(fund_df[,-ncol(fund_df)],1,function(x){Years_list[min(which(x <= 0))]})
        max_years <- apply(fund_df[,-ncol(fund_df)],1,function(x){Years_list[which.max(x)]})
        return(ggplot() +
            
            geom_line(data = df_long_fund, aes(x = year, y = Fund, group = iteration), 
                      color = "grey80", alpha = 0.5) +coord_cartesian(ylim = c(-5*10^15,5 * 10^15 ))
            +geom_vline(xintercept = median(zero_years,na.rm=TRUE),color = "orange", linetype = "dashed", )+
            annotate("rect",xmin=quantile(zero_years,0.025,na.rm=TRUE),xmax=quantile(zero_years,0.975,na.rm=TRUE),ymin=-Inf,ymax=Inf,fill="lightpink",alpha=0.3)+
            geom_vline(xintercept = quantile(zero_years,0.025,na.rm=TRUE), color = "red",size = 0.8 ,linetype = "dotted") +
            geom_vline(xintercept = quantile(zero_years,0.975,na.rm=TRUE), color = "red",size = 0.8, linetype = "dotted") +
            geom_hline(yintercept = 0,color="orange", linetype = "dashed", size = 1)+
            
            geom_vline(xintercept = median(max_years),color = "darkgreen", linetype = "dashed", )+
            
            geom_ribbon(data = summary_df_fund, aes(x = year, ymin = lower, ymax = upper), 
                        fill = "lightblue", alpha = 0.4) +
            
            scale_x_continuous(limits = c(2023, 2023+Forecast_year)) +
            geom_line(data = summary_df_fund, aes(x = year, y = median), 
                      color = "blue", size = 1.2) +
            labs(
                title = "National Pension Fund Simulation",
                subtitle = paste0(
                    "Median Depletion Year: ", round(median(zero_years, na.rm = TRUE)), 
                    " (95% CI: ", round(quantile(zero_years, 0.025, na.rm = TRUE)), " - ", round(quantile(zero_years, 0.975, na.rm = TRUE)), ")",
                    "\nMedian Maximum Year: ", round(median(max_years, na.rm = TRUE)), 
                    " (95% CI: ", round(quantile(max_years, 0.025, na.rm = TRUE)), " - ", round(quantile(max_years, 0.975, na.rm = TRUE)), ")"
                ),
                x = "Year",
                y = "Fund Balance (KRW)"
            ) +
            theme_minimal() )
    }
    
    
    
    if (OPTION == "heatmap") {
        fund_df <- as.data.frame(t(fund_data))
        colnames(fund_df) <- as.character(Years_list)
        fund_df$iteration <- paste0("sim", 1:Iteration)
        df_long_fund <- pivot_longer(fund_df, cols = -iteration, names_to = "year", values_to = "Fund")
        df_long_fund$year <- as.integer(df_long_fund$year)
        summary_df_fund <- df_long_fund %>%
            group_by(year) %>%
            summarise(
                median = median(Fund),
                lower = quantile(Fund, 0.025),
                upper = quantile(Fund, 0.975)
            )
        zero_years <- apply(fund_df[,-ncol(fund_df)],1,function(x){Years_list[min(which(x <= 0))]})

        return(median(zero_years, na.rm = TRUE))
    }
    
    }
    

# Sensitivity Analysis for economic variables
avg_growth_vec <- c(2, 3, 4, 5)
avg_investment_vec <- c(2, 3, 4, 5)


heatmap_df <- expand.grid(avg_growth = avg_growth_vec, avg_investment = avg_investment_vec)
heatmap_df$depletion_median <- NA


for (i in 1:nrow(heatmap_df)) {
    g <- heatmap_df$avg_growth[i]
    inv <- heatmap_df$avg_investment[i]
    
    heatmap_df$depletion_median[i] <- Simu_fund(
        Premium_rate = 0.09,
        Replacement_rate=0.4,
        Start_age = 65,
        avg_growth = g,
        avg_investment = inv,
        OPTION = "heatmap"
    )
}

# heatmap 
ggplot(heatmap_df, aes(x = factor(avg_growth), y = factor(avg_investment), fill = depletion_median)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(depletion_median, 0)), color = "black", size = 4, na.rm = TRUE) +
    
    scale_fill_gradientn(
        colours = c("#deebf7", "#9ecae1", "#3182bd", "#08519c"),
        na.value = "grey90"
    )   + labs(
        title = "Median Depletion Year by Wage Growth and Investment Rate",
        x = "Average Wage Growth (%)",
        y = "Average Investment Return (%)",
        fill = "Depletion Year"
    ) +
    theme_minimal()

# Premium_rate × Replacement_rate × Start_age
premium_vec <- c(0.09, 0.11, 0.13, 0.15)
replace_vec <- c(0.35, 0.4, 0.45, 0.5)
start_age_vec <- c(65, 66, 67)


heatmap_df3 <- expand.grid(
    Premium_rate = premium_vec,
    Replacement_rate = replace_vec,
    Start_age = start_age_vec
)
heatmap_df3$depletion_median <- NA


for (i in 1:nrow(heatmap_df3)) {
    prem <- heatmap_df3$Premium_rate[i]
    rep <- heatmap_df3$Replacement_rate[i]
    sa <- heatmap_df3$Start_age[i]
    
    heatmap_df3$depletion_median[i] <- Simu_fund(
        Premium_rate = prem,
        Start_age = sa,
        avg_growth = 3.5,
        avg_investment = 4,
        Replacement_rate = rep,
        OPTION = "heatmap"
    )
}



ggplot(heatmap_df3, aes(x = factor(Premium_rate), y = factor(Replacement_rate), fill = depletion_median)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(depletion_median, 0)), color = "black", size = 3.5, na.rm = TRUE) +
    scale_fill_gradientn(
        colours = c("#deebf7", "#9ecae1", "#3182bd", "#08519c"),
        na.value = "grey90"
    ) +
    facet_wrap(~ Start_age) +
    labs(
        title = "Median Depletion Year by Contribution, Replacement and Start Age",
        x = "Contribution Rate",
        y = "Replacement Rate",
        fill = "Depletion Year"
    ) +
    theme_minimal()


# Macro Slide Effect
result <- Simu_fund(OPTION="Macroslide",detail = TRUE)
Result_mat <- result$Result_mat
generated_cpi <- result$generated_cpi
payer_forecast <- result$payer_forecast_mat
lifeexp_forecast <- result$lifeexp_mat

# adjustment
Adjustment_list <- -colMeans(payer_forecast) + colMeans(lifeexp_forecast)
CPI_vector <- colMeans(generated_cpi)
CPI_vector<-CPI_vector[2:length(CPI_vector)]
Years <- 2025:(2023 + Forecast_year)

# Trigger
triggered <- CPI_vector > Adjustment_list
slide_effect <- pmax(0, CPI_vector - Adjustment_list)
cumulative_effect <- cumsum(slide_effect)

df1 <- data.frame(
    Year = Years,
    CPI = CPI_vector,
    Adjustment = Adjustment_list,
    Trigger = triggered
)
df1$SlideEffect <- slide_effect
df1$Triggered <- df1$SlideEffect > 0

ggplot(df1, aes(x = Year)) +
    geom_line(aes(y = CPI, color = "CPI"), size = 1.2) +
    geom_line(aes(y = Adjustment, color = "Adjustment"), size = 1.2) +
    geom_point(
        data = df1[df1$Triggered, ], 
        aes(y = CPI), 
        shape = 21, 
        fill = "#ff5733",    
        color = "#b30000",     
        size = 2,
        stroke = 1,
        alpha = 1
    ) +
    geom_ribbon(
        data = df1[df1$Triggered, ],
        aes(ymin = Adjustment, ymax = CPI), 
        fill = "#c6dbef", 
        alpha = 0.4
    ) +
    scale_color_manual(values = c("CPI" = "#08519c", "Adjustment" = "#3182bd")) +
    labs(
        title = "Macro Slide Trigger Points",
        y = "CPI/Rate",
        color = NULL
    ) +
    theme_minimal()

no_slide_result <- Simu_fund(OPTION = "fund", detail = TRUE)
fund_no_slide <- no_slide_result
fund_slide <- result$fund_data
get_depletion_years <- function(fund_mat, years) {
    apply(fund_mat, 2, function(f) {
        depletion_idx <- which(f <= 0)
        if (length(depletion_idx) == 0) return(NA)
        years[min(depletion_idx)]
    })
}

depletion_noslide <- get_depletion_years(fund_no_slide, Years_list)
depletion_slide   <- get_depletion_years(fund_slide, Years_list)

get_depletion_summary <- function(depletion_years) {
    data.frame(Mean = mean(depletion_years, na.rm = TRUE),
               Median = median(depletion_years, na.rm = TRUE),
               Lower_95 = quantile(depletion_years, 0.025, na.rm = TRUE),
               Upper_95 = quantile(depletion_years, 0.975, na.rm = TRUE)
    )
    
}


summary_noslide <- get_depletion_summary(depletion_noslide)
summary_slide   <- get_depletion_summary(depletion_slide)


comparison_summary <- rbind(
    NoMacroslide = summary_noslide,
    WithMacroslide = summary_slide
)


print(round(comparison_summary, 1))
write.csv(round(comparison_summary, 1),"C:/Users/leeha/OneDrive/Desktop/Affine mortality/fig/slide.csv")

 

