library(devtools)
library(AffineMortality)
library(tidyr)
library(dplyr)
RatestoAvg<-function(mu) 
{
    mu_bar <- mu
    for (row in 1:nrow(mu)) {
        for (col in 1:ncol(mu)) {
            mu_bar[row, col] <- mean(mu[1:row, col])
        }
    }
    return(mu_bar)
}
# bring data
mxt_data <- read.table("C:/Users/leeha/OneDrive/Desktop/Affine mortality/Affine mortality/Mx_1x1.txt",header=TRUE,skip=2)

# data processing
mxt_data <- mxt_data %>% mutate_at(vars(Year,Age,Female,Male,Total),as.numeric)

mxt_total <- mxt_data %>%filter(Age >= 40 & Age <= 109)%>% select(Age,Year,Total)%>% pivot_wider(names_from = Year,values_from = Total)
mxt_total<-as.data.frame(mxt_total)
rownames(mxt_total)<-mxt_total$Age
mxt_total<-mxt_total[,-1]

mxt_female<-mxt_data %>%filter(Age >= 40 & Age <= 109)%>% select(Age,Year,Female)%>% pivot_wider(names_from = Year,values_from = Female)
mxt_female<-as.data.frame(mxt_female)
rownames(mxt_female)<-mxt_female$Age
mxt_female<-mxt_female[,-1]

mxt_male<-mxt_data %>% filter(Age >= 40 & Age <= 109)%>%select(Age,Year,Male)%>% pivot_wider(names_from = Year,values_from = Male)
mxt_male<-as.data.frame(mxt_male)
rownames(mxt_male)<-mxt_male$Age
mxt_male<-mxt_male[,-1]
# 0 to 39
mxt_0_39<-mxt_data %>%filter(Age >= 0 & Age <= 39)%>% select(Age,Year,Total)%>% pivot_wider(names_from = Year,values_from = Total)
mxt_0_39<-as.data.frame(mxt_0_39)
rownames(mxt_0_39)<-mxt_total$Age
mxt_0_39<-mxt_0_39[,-1]
mu_bar_young <- RatestoAvg(mxt_0_39)
Sx_young<- exp(-mu_bar_young[,"2023"])


mu_bar_tot<-RatestoAvg(mxt_total)
mu_bar_fem<-RatestoAvg(mxt_female)
mu_bar_man <- RatestoAvg(mxt_male)

# Total 
#Blackburn-Sherris model(BS) - independent
t_bsd_3fi<-AffineMortality::affine_fit(model="BS",fact_dep = FALSE, n_factors=3, data=mu_bar_tot,
                       st_val=sv_default$BSi,max_iter=200,tolerance=0.1)
#Blackburn-Sherris model(BS) - three dependent
t_bsd_3f<-AffineMortality::affine_fit(model="BS",fact_dep = TRUE, n_factors=3, data=mu_bar_tot,
                      st_val=sv_default$BSd,max_iter=200,tolerance=0.1)

#Arbigrage Free Nelson Siegel (AFNS) - independent
t_AFNSi<-AffineMortality::affine_fit(model="AFNS",fact_dep = FALSE, n_factors=3, data=mu_bar_tot,
                     st_val=sv_default$AFNSi,max_iter=200,tolerance=0.1)
#Arbigrage Free Nelson Siegel (AFNS) - dependent
t_AFNS<-AffineMortality::affine_fit(model="AFNS",fact_dep = TRUE, n_factors=3, data=mu_bar_tot,
                    st_val=sv_default$AFNSd,max_iter=200,tolerance=0.1)

t_bsd_3fi$AIC
t_bsd_3fi$BIC

t_bsd_3f$AIC
t_bsd_3f$BIC

t_AFNSi$AIC
t_AFNSi$BIC

t_AFNS$AIC
t_AFNS$BIC

t_fitted_bsdi <- mubar_hat(model="BS", fact_dep=FALSE,n_factors=3,parameters=t_bsd_3fi$fit$par_est,data=mu_bar_tot)
t_fitted_bsd <- mubar_hat(model="BS", fact_dep=TRUE,n_factors=3,parameters=t_bsd_3f$fit$par_est,data=mu_bar_tot)
t_fitted_afnsi <- mubar_hat(model="AFNS", fact_dep=FALSE,n_factors=3,parameters=t_AFNSi$fit$par_est,data=mu_bar_tot)
t_fitted_afns <- mubar_hat(model="AFNS", fact_dep=TRUE,n_factors=3,parameters=t_AFNS$fit$par_est,data=mu_bar_tot)

#RMSE
sum((mu_bar_tot-t_fitted_bsdi)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot))
sum((mu_bar_tot-t_fitted_bsd)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot))
sum((mu_bar_tot-t_fitted_afnsi)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot))
sum((mu_bar_tot-t_fitted_afns)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot))

future_year=10
bsdi_3f_proj <- affine_project(model="BS",fact_dep = FALSE,n_factors=3,parameters=t_bsd_3fi$fit$par_est,data=mu_bar_tot,years_proj=future_year)
bsd_3f_proj <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=t_bsd_3f$fit$par_est,data=mu_bar_tot,years_proj=future_year)
afnsi_3f_proj <- AffineMortality::affine_project(model="AFNS",fact_dep = FALSE,n_factors=3,parameters=t_AFNSi$fit$par_est,data=mu_bar_tot,years_proj=future_year)
afns_3f_proj <- affine_project(model="AFNS",fact_dep = TRUE,n_factors=3,parameters=t_AFNS$fit$par_est,data=mu_bar_tot,years_proj=future_year)
plot(rownames(mu_bar_tot),bsdi_3f_proj,type = "l",ylim=c(0,2),ylab="S(t)",xlab="Age",lwd=2,main = paste("Year",2023+future_year))
lines(rownames(mu_bar_tot),bsd_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darkkhaki",lwd=2)
lines(rownames(mu_bar_tot),afnsi_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darkmagenta",lwd=2)
lines(rownames(mu_bar_tot),afns_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darksalmon",lwd=2)
legend("bottomleft",legend=c("BS-ind-3","BS-dep-3","AFNS-ind","AFNS-dep"),cex=0.7,lty = c(1,1,1,1),lwd=c(2,2,2,2),col=c("black","darkkhaki","darkmagenta","darksalmon"))
plot(rownames(mu_bar_tot),afnsi_3f_proj,type = "l",ylim=c(0,2),ylab="S(t)",xlab="Age",lwd=2,main = paste("Year",2023+future_year))


t_survival_curve <- function(future_year){
    bsdi_3f_proj <- affine_project(model="BS",fact_dep = FALSE,n_factors=3,parameters=t_bsd_3fi$fit$par_est,data=mu_bar_tot,years_proj=future_year)
    bsd_3f_proj <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=t_bsd_3f$fit$par_est,data=mu_bar_tot,years_proj=future_year)
    afnsi_3f_proj <- affine_project(model="AFNS",fact_dep = FALSE,n_factors=3,parameters=t_AFNSi$fit$par_est,data=mu_bar_tot,years_proj=future_year)
    afns_3f_proj <- affine_project(model="AFNS",fact_dep = TRUE,n_factors=3,parameters=t_AFNS$fit$par_est,data=mu_bar_tot,years_proj=future_year)
    plot(rownames(mu_bar_tot),bsdi_3f_proj,type = "l",ylim=c(0,2),ylab="S(t)",xlab="Age",lwd=2,main = paste("Year",2023+future_year))
    lines(rownames(mu_bar_tot),bsd_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darkkhaki",lwd=2)
    lines(rownames(mu_bar_tot),afnsi_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darkmagenta",lwd=2)
    lines(rownames(mu_bar_tot),afns_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darksalmon",lwd=2)
    legend("bottomleft",legend=c("BS-ind-3","BS-dep-3","AFNS-ind","AFNS-dep"),cex=0.7,lty = c(1,1,1,1),lwd=c(2,2,2,2),col=c("black","darkkhaki","darkmagenta","darksalmon"))
    
}
par(mfrow=c(1,3),oma=c(0,0,2,0))
t_survival_curve(10)
t_survival_curve(20)
t_survival_curve(30)
mtext("Predicted Survival Curve",outer=TRUE,cex=1.2,line=0.1)

prob_neg_mu(model="BS", fact_dep=TRUE, n_factors=3,
            parameters=t_bsd_3f$fit$par_est, data=mu_bar_tot, years_proj=30,
            n_simulations=1000)
prob_neg_mu(model="BS", fact_dep=FALSE, n_factors=3,
            parameters=t_bsd_3fi$fit$par_est, data=mu_bar_tot, years_proj=30,
            n_simulations=1000)
prob_neg_mu(model="AFNS", fact_dep=TRUE, n_factors=3,
            parameters=t_AFNS$fit$par_est, data=mu_bar_tot, years_proj=30,
            n_simulations=1000)
prob_neg_mu(model="AFNS", fact_dep=FALSE, n_factors=3,
            parameters=t_AFNSi$fit$par_est, data=mu_bar_tot, years_proj=30,
            n_simulations=1000)

par(mfrow=c(1,1),oma=c(0,0,0,0))

prj <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=t_bsd_3f$fit$par_est,data=mu_bar_tot,1)
plot(rownames(mu_bar_tot),prj,col=rainbow(7)[1],type = "l",ylab="S(x)",xlab="Age",lwd=2,main="S(x) By BS 3 Dep. From 2024 to 2053")
year<-c(5,10,20,30)
for (k in 1:4){
    prj <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=t_bsd_3f$fit$par_est,data=mu_bar_tot,year[k])
    lines(rownames(mu_bar_tot),prj,col=rainbow(7)[k+3],type = "l",ylab="S(x)",xlab="Age",lwd=2)
    
}
legend("bottomleft",legend=c("1 yr","5 yr","10 yr","20 yr","30 yr"),col=c(rainbow(7)[1],rainbow(7)[4],rainbow(7)[5],rainbow(7)[6],rainbow(7)[7]),lty = c(1,1,1,1,1,1),lwd=c(2,2,2,2,2))

mbsd_3f_proj1 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=1)
mbsd_3f_proj10 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=10)
mbsd_3f_proj30 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=30)
mbsd_3f_proj50 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=50)
mbsd_3f_proj70 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=70)
mbsd_3f_proj100 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=100)

plot(rownames(mu_bar_tot),mbsd_3f_proj1,type = "l",xlim=c(45,100),ylim=c(0,1.2),ylab="S(t)",xlab="Age",lwd=2,main = "2024 ~ 2124 Male Survival Curve",col="red")
title(sub="BS 3-factor Dependent Model")
lines(rownames(mu_bar_man),mbsd_3f_proj10,type = "l",ylab="S(t)",xlab="Age",col="orange",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj30,type = "l",ylab="S(t)",xlab="Age",col="yellow",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj50,type = "l",ylab="S(t)",xlab="Age",col="green",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj70,type = "l",ylab="S(t)",xlab="Age",col="blue",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj100,type = "l",ylab="S(t)",xlab="Age",col="darkslateblue",lwd=2)
legend("bottomleft",legend=c("1 yr","10 yr","30 yr","50 yr","70 yr","100 yr"),col=c("red","orange","yellow","green","blue","darkslateblue"),lty = c(1,1,1,1,1,1))


# Women case
#Blackburn-Sherris model(BS) - independent
fm_bsd_3fi<-AffineMortality::affine_fit(model="BS",fact_dep = FALSE, n_factors=3, data=mu_bar_fem,
           st_val=sv_default$BSi,max_iter=200,tolerance=0.1)

#Blackburn-Sherris model(BS) - three dependent
fm_bsd_3f<-AffineMortality::affine_fit(model="BS",fact_dep = TRUE, n_factors=3, data=mu_bar_fem,
                      st_val=sv_default$BSd,max_iter=200,tolerance=0.1)

#Arbigrage Free Nelson Siegel (AFNS) - independent
fm_AFNSi<-AffineMortality::affine_fit(model="AFNS",fact_dep = FALSE, n_factors=3, data=mu_bar_fem,
                      st_val=sv_default$AFNSi,max_iter=200,tolerance=0.1)
#Arbigrage Free Nelson Siegel (AFNS) - dependent
fm_AFNS<-AffineMortality::affine_fit(model="AFNS",fact_dep = TRUE, n_factors=3, data=mu_bar_fem,
                     st_val=sv_default$AFNSd,max_iter=200,tolerance=0.1)

fm_bsd_3fi$AIC
fm_bsd_3fi$BIC

fm_bsd_3f$AIC
fm_bsd_3f$BIC

fm_AFNSi$AIC
fm_AFNSi$BIC

fm_AFNS$AIC
fm_AFNS$BIC

fm_fitted_bsdi <- mubar_hat(model="BS", fact_dep=FALSE,n_factors=3,parameters=fm_bsd_3fi$fit$par_est,data=mu_bar_fem)
fm_fitted_bsd <- mubar_hat(model="BS", fact_dep=TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_fem)
fm_fitted_afnsi <- mubar_hat(model="AFNS", fact_dep=FALSE,n_factors=3,parameters=fm_AFNSi$fit$par_est,data=mu_bar_fem)
fm_fitted_afns <- mubar_hat(model="AFNS", fact_dep=TRUE,n_factors=3,parameters=fm_AFNS$fit$par_est,data=mu_bar_fem)

#RMSE
sum((mu_bar_fem-fm_fitted_bsdi)^2)/(nrow(mu_bar_fem)*ncol(mu_bar_fem))
sum((mu_bar_fem-fm_fitted_bsd)^2)/(nrow(mu_bar_fem)*ncol(mu_bar_fem))
sum((mu_bar_fem-fm_fitted_afnsi)^2)/(nrow(mu_bar_fem)*ncol(mu_bar_fem))
sum((mu_bar_fem-fm_fitted_afns)^2)/(nrow(mu_bar_fem)*ncol(mu_bar_fem))

# S(t) forecasts

# future survival curve
survival_curve <- function(future_year){
    bsdi_3f_proj <- affine_project(model="BS",fact_dep = FALSE,n_factors=3,parameters=fm_bsd_3fi$fit$par_est,data=mu_bar_fem,years_proj=future_year)
    bsd_3f_proj <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_fem,years_proj=future_year)
    afnsi_3f_proj <- affine_project(model="AFNS",fact_dep = FALSE,n_factors=3,parameters=fm_AFNSi$fit$par_est,data=mu_bar_fem,years_proj=future_year)
    afns_3f_proj <- affine_project(model="AFNS",fact_dep = TRUE,n_factors=3,parameters=fm_AFNS$fit$par_est,data=mu_bar_fem,years_proj=future_year)
    plot(rownames(mu_bar_fem),bsdi_3f_proj,type = "l",ylim=c(0,2),ylab="S(t)",xlab="Age",lwd=2,main = paste("Year",2023+future_year))
    lines(rownames(mu_bar_fem),bsd_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darkkhaki",lwd=2)
    lines(rownames(mu_bar_fem),afnsi_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darkmagenta",lwd=2)
    lines(rownames(mu_bar_fem),afns_3f_proj,type = "l",ylab="S(t)",xlab="Age",col="darksalmon",lwd=2)
    legend("bottomleft",legend=c("BS-ind-3","BS-dep-3","AFNS-ind","AFNS-dep"),cex=0.7,lty = c(1,1,1,1),lwd=c(2,2,2,2),col=c("black","darkkhaki","darkmagenta","darksalmon"))

}
par(mfrow=c(1,3),oma=c(0,0,2,0))
survival_curve(10)
survival_curve(30)
survival_curve(50)
mtext("Predicted Female Survival Curve",outer=TRUE,cex=1.2,line=0.3)

prob_neg_mu(model="BS", fact_dep=TRUE, n_factors=3,
            parameters=fm_bsd_3f$fit$par_est, data=mu_bar_fem, years_proj=10,
            n_simulations=1000)
prob_neg_mu(model="BS", fact_dep=FALSE, n_factors=3,
            parameters=fm_bsd_3fi$fit$par_est, data=mu_bar_fem, years_proj=10,
            n_simulations=1000)
prob_neg_mu(model="AFNS", fact_dep=TRUE, n_factors=3,
            parameters=fm_AFNS$fit$par_est, data=mu_bar_fem, years_proj=10,
            n_simulations=1000)
prob_neg_mu(model="AFNS", fact_dep=FALSE, n_factors=3,
            parameters=fm_AFNSi$fit$par_est, data=mu_bar_fem, years_proj=10,
            n_simulations=1000)

# men case
#Blackburn-Sherris model(BS) - independent
m_bsd_3fi<-AffineMortality::affine_fit(model="BS",fact_dep = FALSE, n_factors=3, data=mu_bar_man,
                       st_val=fm_bsd_3fi$fit$par_est,max_iter=200,tolerance=0.1)
m_bsd_4fi<-AffineMortality::affine_fit(model="BS",fact_dep = FALSE, n_factors=4, data=mu_bar_man,
                      st_val=fm_bsd_4fi$fit$par_est,max_iter=200,tolerance=0.1)

#Blackburn-Sherris model(BS) - three dependent
m_bsd_3f<-AffineMortality::affine_fit(model="BS",fact_dep = TRUE, n_factors=3, data=mu_bar_man,
                      st_val=sv_default$BSd,max_iter=15,tolerance=0.01)
notpositive_m_parameters <- m_bsd_3f$fit$par_est

m_bsd_3f<-AffineMortality::affine_fit(model="BS",fact_dep = TRUE, n_factors=3, data=mu_bar_man,
                     st_val=fm_bsd_3f$fit$par_est,max_iter=10,tolerance=0.1)

#Arbitrage Free Nelson Siegel (AFNS) - independent
m_AFNSi<-AffineMortality::affine_fit(model="AFNS",fact_dep = FALSE, n_factors=3, data=mu_bar_man,
                     st_val=fm_AFNSi$fit$par_est,max_iter=200,tolerance=0.1)
#Arbitrage Free Nelson Siegel (AFNS) - dependent
m_AFNS<-AffineMortality::affine_fit(model="AFNS",fact_dep = TRUE, n_factors=3, data=mu_bar_man,
                    st_val=fm_AFNS$fit$par_est,max_iter=200,tolerance=0.1)


m_bsd_3fi$AIC
m_bsd_3fi$BIC

m_bsd_3f$AIC
m_bsd_3f$BIC

m_AFNSi$AIC
m_AFNSi$BIC

m_AFNS$AIC
m_AFNS$BIC


#fitting
m_fitted_bsd3fi <- mubar_hat(model="BS", fact_dep=FALSE,n_factors=3,parameters=m_bsd_3fi$fit$par_est,data=mu_bar_man)
m_fitted_bsd3f <- mubar_hat(model="BS", fact_dep=TRUE,n_factors=3,parameters=m_bsd_3f$fit$par_est,data=mu_bar_man)
m_fitted_afnsi <- mubar_hat(model="AFNS", fact_dep=FALSE,n_factors=3,parameters=m_AFNSi$fit$par_est,data=mu_bar_man)
m_fitted_afns <- mubar_hat(model="AFNS", fact_dep=TRUE,n_factors=3,parameters=m_AFNS$fit$par_est,data=mu_bar_man)


#RMSE
sum((mu_bar_man-m_fitted_bsd3fi)^2)/(nrow(mu_bar_man)*ncol(mu_bar_man))
sum((mu_bar_man-m_fitted_bsd3f)^2)/(nrow(mu_bar_man)*ncol(mu_bar_man))
sum((mu_bar_man-m_fitted_afnsi)^2)/(nrow(mu_bar_man)*ncol(mu_bar_man))
sum((mu_bar_man-m_fitted_afns)^2)/(nrow(mu_bar_man)*ncol(mu_bar_man))

models <- c("BS ind.3","BS dep.","AFNS ind","AFNS dep")
aic_values_t <- c(t_bsd_3fi$AIC,t_bsd_3f$AIC,t_AFNSi$AIC,t_AFNS$AIC)
bic_values_t <-c(t_bsd_3fi$BIC,t_bsd_3f$BIC,t_AFNSi$BIC,t_AFNS$BIC)

rmse_values_t<-c(sum((mu_bar_tot-t_fitted_bsdi)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot)),
                  sum((mu_bar_tot-t_fitted_bsd)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot))
                  ,sum((mu_bar_tot-t_fitted_afnsi)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot))
                  ,sum((mu_bar_tot-t_fitted_afns)^2)/(nrow(mu_bar_tot)*ncol(mu_bar_tot))
)
comp_t <- data.frame(Model=models,AIC=aic_values_t,BIC=bic_values_t,RMSE=rmse_values_t)
comp_t
write.csv(comp_t,"C:/Users/leeha/OneDrive/Desktop/Affine mortality/fig/comp_t.csv")

#project
bsd_3f_proj_man <- affine_project(model="BS",fact_dep = FALSE,n_factors=3,parameters=m_bsd_3fi$fit$par_est,data=mu_bar_man,years_proj=1)
plot(rownames(mu_bar_man),bsd_3f_proj_man,type = "l",ylab="S(t)",xlab="Age")

bsd_3f_proj_man_depen <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=1)
plot(rownames(mu_bar_man),bsd_3f_proj_man_depen,type = "l",ylab="S(t)",xlab="Age")

afnsi_3f_proj_man <- affine_project(model="AFNS",fact_dep = FALSE,n_factors=3,parameters=m_AFNSi$fit$par_est,data=mu_bar_man,years_proj=1)
plot(rownames(mu_bar_man),afnsi_3f_proj_man,type = "l",ylab="S(t)",xlab="Age")

afns_3f_proj_man <- affine_project(model="AFNS",fact_dep = TRUE,n_factors=3,parameters=m_AFNS$fit$par_est,data=mu_bar_man,years_proj=1)
plot(rownames(mu_bar_man),afns_3f_proj_man,type = "l",ylab="S(t)",xlab="Age")

plot(rownames(mu_bar_man),bsd_3f_proj_man,type = "l",ylim=c(0,1.2),ylab="S(t)",xlab="Age",lwd=2,main = "Year 2024 - Predicted Male Survival Curve")
lines(rownames(mu_bar_man),bsd_3f_proj_man_depen,type = "l",ylab="S(t)",xlab="Age",col="darkkhaki",lwd=2)
lines(rownames(mu_bar_man),afnsi_3f_proj_man,type = "l",ylab="S(t)",xlab="Age",col="darkmagenta",lwd=2)
lines(rownames(mu_bar_man),afns_3f_proj_man,type = "l",ylab="S(t)",xlab="Age",col="darksalmon",lwd=2)
legend("bottomleft",legend=c("BS-ind-3","BS-ind-4","BS-dep-3","AFNS-ind","AFNS-dep"),lty = c(1,1,1,1,1),lwd=c(2,2,2,2,2),col=c("black","cornflowerblue","darkkhaki","darkmagenta","darksalmon"))


survival_curve_man <- function(year){
    bsd_3f_proj_man <- affine_project(model="BS",fact_dep = FALSE,n_factors=3,parameters=m_bsd_3fi$fit$par_est,data=mu_bar_man,years_proj=year)
    bsd_3f_proj_man_depen <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=m_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=year)
    afnsi_3f_proj_man <- affine_project(model="AFNS",fact_dep = FALSE,n_factors=3,parameters=m_AFNSi$fit$par_est,data=mu_bar_man,years_proj=year)
    afns_3f_proj_man <- affine_project(model="AFNS",fact_dep = TRUE,n_factors=3,parameters=m_AFNS$fit$par_est,data=mu_bar_man,years_proj=year)
    plot(rownames(mu_bar_man),bsd_3f_proj_man,type = "l",ylim=c(0,1.2),ylab="S(t)",xlab="Age",lwd=2,main = paste("Year", 2023+year))
    lines(rownames(mu_bar_man),bsd_4f_proj_man,type = "l",ylab="S(t)",xlab="Age",col="cornflowerblue",lwd=2)
    lines(rownames(mu_bar_man),bsd_3f_proj_man_depen,type = "l",ylab="S(t)",xlab="Age",col="darkkhaki",lwd=2)
    lines(rownames(mu_bar_man),afnsi_3f_proj_man,type = "l",ylab="S(t)",xlab="Age",col="darkmagenta",lwd=2)
    lines(rownames(mu_bar_man),afns_3f_proj_man,type = "l",ylab="S(t)",xlab="Age",col="darksalmon",lwd=2)
    legend("bottomleft",legend=c("BS-ind-3","BS-ind-4","BS-dep-3","AFNS-ind","AFNS-dep"),cex=0.6,lty = c(1,1,1,1,1),lwd=c(2,2,2,2,2),col=c("black","cornflowerblue","darkkhaki","darkmagenta","darksalmon"))
    
}


par(mfrow=c(2,3),oma=c(0,0,4,0))
survival_curve_man(1)
survival_curve_man(10)
survival_curve_man(30)
survival_curve_man(50)
survival_curve_man(70)
survival_curve_man(100)
mtext("Predicted Male Survival Curve",outer=TRUE,cex=1.5,line=1)

par(mfrow=c(1,1),oma=c(0,0,0,0))
mbsd_3f_proj1 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=1)
mbsd_3f_proj10 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=10)
mbsd_3f_proj30 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=30)
mbsd_3f_proj50 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=50)
mbsd_3f_proj70 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=70)
mbsd_3f_proj100 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_man,years_proj=100)

plot(rownames(mu_bar_man),mbsd_3f_proj1,type = "l",xlim=c(45,100),ylim=c(0,1.2),ylab="S(t)",xlab="Age",lwd=2,main = "2024 ~ 2124 Male Survival Curve",col="red")
title(sub="BS 3-factor Dependent Model")
lines(rownames(mu_bar_man),mbsd_3f_proj10,type = "l",ylab="S(t)",xlab="Age",col="orange",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj30,type = "l",ylab="S(t)",xlab="Age",col="yellow",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj50,type = "l",ylab="S(t)",xlab="Age",col="green",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj70,type = "l",ylab="S(t)",xlab="Age",col="blue",lwd=2)
lines(rownames(mu_bar_man),mbsd_3f_proj100,type = "l",ylab="S(t)",xlab="Age",col="darkslateblue",lwd=2)
legend("bottomleft",legend=c("1 yr","10 yr","30 yr","50 yr","70 yr","100 yr"),col=c("red","orange","yellow","green","blue","darkslateblue"),lty = c(1,1,1,1,1,1))
#BS_3f_plot
bsd_3f_proj30 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_fem,years_proj=30)
bsd_3f_proj50 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_fem,years_proj=50)
bsd_3f_proj100 <- affine_project(model="BS",fact_dep = TRUE,n_factors=3,parameters=fm_bsd_3f$fit$par_est,data=mu_bar_fem,years_proj=100)
plot(rownames(mu_bar_fem),bsd_3f_proj30,type = "l",ylim=c(0,2),ylab="S(t)",xlab="Age",lwd=2,main = "Year 2024 - Predicted Female Survival Curve")
lines(rownames(mu_bar_fem),bsd_3f_proj50,type = "l",ylab="S(t)",xlab="Age",col="cornflowerblue",lwd=2)
lines(rownames(mu_bar_fem),bsd_3f_proj100,type = "l",ylab="S(t)",xlab="Age",col="darkkhaki",lwd=2)

# mape
?MAPE_age
plot(MAPE_age(mu_bar_fem,fm_fitted_bsd),type = "l",ylim = c(0,2))
lines(MAPE_age(mu_bar_fem,fm_fitted_bsdi),col="red")
lines(MAPE_age(mu_bar_fem,fm_fitted_bsdi4),col="purple")

lines(MAPE_age(mu_bar_fem,fm_fitted_afns),col="blue")
lines(MAPE_age(mu_bar_fem,fm_fitted_afnsi),col="green")

plot(MAPE_age(mu_bar_man,m_fitted_bsd3f),type = "l",ylim = c(0,2))
lines(MAPE_age(mu_bar_man,m_fitted_bsd3fi),col="red")
lines(MAPE_age(mu_bar_man,m_fitted_afns),col="blue")
lines(MAPE_age(mu_bar_man,m_fitted_afnsi),col="green")

fm_std_red<-std_res(model="BS",fact_dep = TRUE, n_factors = 3,parameters =fm_bsd_3f$fit$par_est,data=mu_bar_fem)
heatmap_res(residuals = fm_std_red)

fm_std_redi<-std_res(model="BS",fact_dep = FALSE, n_factors = 3,parameters =fm_bsd_3fi$fit$par_est,data=mu_bar_fem)
heatmap_res(residuals = fm_std_redi)

fm_std_redi4<-std_res(model="BS",fact_dep = FALSE, n_factors = 4,parameters =fm_bsd_4fi$fit$par_est,data=mu_bar_fem)
heatmap_res(residuals = fm_std_redi4)

m_std_red<-std_res(model="BS",fact_dep = TRUE, n_factors = 3,parameters =m_bsd_3f$fit$par_est,data=mu_bar_man)
heatmap_res(residuals = m_std_red)

m_std_redi<-std_res(model="BS",fact_dep = FALSE, n_factors = 3,parameters =m_bsd_3fi$fit$par_est,data=mu_bar_man)
heatmap_res(residuals = m_std_redi)

