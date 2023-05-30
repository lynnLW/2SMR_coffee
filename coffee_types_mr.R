#library package
library(TwoSampleMR)

#exposure data
coffee_dataset2<-read.csv("Table_2_coffee_type_id.csv",header=T)
id2<-coffee_dataset2$coffee_id
dir.create("exposure_data/types/",recursive = T)
exposure2<-c()
for (n in 1:12){
  dat<-extract_instruments(id2[n],clump=F)
  if (is.null(dat)==FALSE){
    dat<-clump_data(dat,clump_kb = 1000,clump_r2 = 0.01,pop="EUR")
    dat<-add_rsq(dat)
    dat$F<-dat$rsq.exposure*(dat$effective_n.exposure-2)/(1-dat$rsq.exposure)*1
    dat<-dat[dat$F>10,]
    exposure2<-rbind(exposure2,dat)
  }
}
save(exposure2,file="exposure_data/types/coffee_types_exposure.Rdata")

#extract outcome
dir.create("outcome_data/types/",recursive = T)
disease_list<-read.csv("diseases.csv")
disease_id<-c(disease_list$endometriosis,disease_list$ovarian_cyst,disease_list$ovarian_cancer)
disease_id<-disease_id[disease_id != ""]
outcome2 <- extract_outcome_data(snps = exposure2$SNP, outcomes = disease_id,proxies = FALSE)
save(outcome2,file="outcome_data/types/types_disease_outcome.Rdata")

#harmonise data
data2<-harmonise_data(exposure_dat = exposure2,
                        outcome_dat = outcome2,action = 2)
#perform 2smr 
res2<-mr(data2)

#calculate odds ratio
res2<-generate_odds_ratios(res2)

#filter sig mr result
sig2<-res2[res2$pval<0.05,]

#filter sig data2 for sensitive analysis
data2_filter<-data2[data2$id.outcome %in% sig2$id.outcome,]
mr_het<-mr_heterogeneity(data2_filter, method_list=c("mr_egger_regression"))
mr_ple<-mr_pleiotropy_test(data2_filter)
res_loo<-mr_leaveoneout(data2_filter)
res_single<-mr_singlesnp(data2_filter)
    
#calculate odds ratio
res_single <- generate_odds_ratios(res_single)

dir.create("mr")
#save results
write.table(sig2,file="mr/types_mr_sig.csv",sep=",",row.names = F)
write.table(res2,file="mr/types_mr.csv",sep=",",row.names = F)
write.table(mr_het,file = "mr/types_hterogeneity.csv",sep=",",row.names = F)
write.table(mr_ple,file = "mr/types_pleiotropy.csv",sep=",",row.names = F)
write.table(res_loo,file = "mr/types_leaveoneout.csv",sep=",",row.names = F)
write.table(res_single,file = "mr/types_single.csv",sep=",",row.names = F)






