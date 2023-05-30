#install package
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")

#library package
library(TwoSampleMR)

#exposure data
coffee_dataset<-read.csv("Table_1_coffee_intake_id.csv",header=T)
id<-coffee_dataset$coffee_id
dir.create("exposure_data/intake/")
exposure<-c()
for (n in 1:2){
  dat<-extract_instruments(id[n],clump=F)
  if (is.null(dat)==FALSE){
    dat<-clump_data(dat,clump_kb = 1000,clump_r2 = 0.01,pop="EUR")
    dat<-add_rsq(dat)
    dat$F<-dat$rsq.exposure*(dat$effective_n.exposure-2)/(1-dat$rsq.exposure)*1
    dat<-dat[dat$F>10,]
    exposure<-rbind(exposure,dat)
  }
}
save(exposure,file="exposure_data/intake/coffee_intake_exposure.Rdata")

##disease outcome
disease_list<-read.csv("diseases.csv")
disease_id<-c(disease_list$pcos,disease_list$endometriosis,disease_list$ovarian_cyst,disease_list$ovarian_cancer)
disease_id<-disease_id[disease_id != ""]
outcome <- extract_outcome_data(snps = exposure$SNP, outcomes = disease_id, proxies = FALSE)
save(outcome,file="outcome_data/intake/intake_disease_outcome.Rdata")

#harmonise data
data<-harmonise_data(exposure_dat = exposure,
                     outcome_dat = outcome,action = 2)
#perform 2smr 
res<-mr(data)

#calculate odds ratio
res<-generate_odds_ratios(res)

#sig mr filter
sig<-res[res$pval<0.05,]

#filter sig data2 for sensitive analysis
data_filter<-data[data$id.outcome %in% sig$id.outcome,]
mr_het<-mr_heterogeneity(data_filter, method_list=c("mr_egger_regression"))
mr_ple<-mr_pleiotropy_test(data_filter)
res_loo<-mr_leaveoneout(data_filter)
res_single<-mr_singlesnp(data_filter)

#calculate odds ratio
res_single <- generate_odds_ratios(res_single)

dir.create("mr")
#save results
write.table(sig,file="mr/intake_mr_sig.csv",sep=",",row.names = F)
write.table(res,file="mr/intake_mr.csv",sep=",",row.names = F)
write.table(mr_het,file = "mr/intake_hterogeneity.csv",sep=",",row.names = F)
write.table(mr_ple,file = "mr/intake_pleiotropy.csv",sep=",",row.names = F)
write.table(res_loo,file = "mr/intake_leaveoneout.csv",sep=",",row.names = F)
write.table(res_single,file = "mr/intake_single.csv",sep=",",row.names = F)
