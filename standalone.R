setwd("/Users/niek/repos/ecgenetics/")
source("global.R")
source("helpers.R")
source("helpers.heatmap.R")
source("helpers.regionalplot.R")
source("get_nearest_gene.r")

library(ggpubr)

# df.static <- df.static[df.static$f %in% c(1,2),]
# save(df.static,file=f_static)
rsids=c("rs1045716","rs145425222")
rsids=c("rs201106462")
rsids=c("rs540538026","rs367896724","rs201106462","rs1045716")
#input<-list()
#input$rsids="rs888578 rs888579 rs888580 rs888627 rs889339 rs889495 rs889802 rs890053 rs890054 rs890055 rs8902 rs890326 rs890327 rs8909 rs890955 rs890959 rs890960 rs890972 rs890973 rs890974 rs890975 rs890976 rs890977 rs890978 rs890979 rs890980 rs890981 rs890982 rs890983 rs890984 rs890985 rs890986 rs891145 rs891173 rs891174 rs891175 rs891176 rs891206 rs891735 rs891877 rs891878 rs891879 rs892201 rs892202 rs892203 rs892204 rs892209 rs892225 rs892291 rs892292 rs892294 rs892295 rs892296 rs892518 rs892520 rs892522 rs892523 rs892524 rs892525 rs892526 rs892527 rs892528 rs892529 rs892530 rs892665 rs893235 rs893237 rs893538 rs893539 rs893540 rs893541 rs893542 rs893543 rs893657 rs893658 rs893659 rs893660 rs893661 rs893870 rs893871 rs893874 rs893875 rs893876 rs893877 rs893878 rs893879 rs894534 rs894628 rs894633 rs895341 rs895342 rs895343 rs895344 rs895345 rs895346 rs895347 rs895656 rs895657 rs895658 rs896522 rs896523 rs896525 rs896526 rs896527 rs896528 rs896529 rs897078 rs897443 rs897444 rs897445 rs897855 rs897856 rs897857 rs897858 rs897860 rs898049 rs898051 rs898053 rs898057 rs898058 rs898060 rs898061 rs898062 rs898671 rs898748 rs898749 rs898750 rs898751 rs899125 rs899195 rs899385 rs899387 rs899388 rs899389 rs899390 rs899391 rs899392 rs899393 rs899397 rs899399 rs899402 rs899405 rs899406 rs899408 rs899409 rs899411"
#input$rsids="rs10797060\nrs12143842\nrs59852339\nrs116381680\nrs11580561\nrs531706\nrs12569204\nrs1200130\nrs67779898\nrs1591734\nrs575784488\nrs340849\nrs6735214\nrs7605146\nrs907683\nrs35394392\nrs141501074\nrs760390566\nrs76521806\nrs1403767\nrs1706003\nrs6533694\nrs61760320\n4:155350842_GT_G\nrs260399\nrs61367953\nrs78144291\nrs730506\nrs28763977\nrs663474\nrs11153730\nrs10457469"

#input$rsids="rs540538026,rs367896724 ,rs201106462"
#rsids=unlist(strsplit( input$rsid ,"\n|[, ]+"))
#rsids=unlist(strsplit( rsids ,"\n|[, ]+"))
rsids=trimws(rsids)

input = "rs140985773,rs202536,rs8117600,rs72484093,rs202535,rs201441911,rs202534,rs202533,rs202532,rs79012224,rs202531,rs756734607,rs202530,rs76531410,rs202529,rs112052313,rs117457337,rs62186956,rs778107293,rs17770828,rs11697150,rs202528,20:1671665_AT_A,rs202527,rs59979485,rs151003361,rs6110951,rs202526,rs117725692,rs62189067,rs113392269,rs149173145,rs543610294,rs200730035,rs112362679,rs202525,rs202524,rs11696890,rs202523,rs202522,rs202521,rs202520,rs2244191,rs7269065,rs7274477,rs76021634,rs2273329,rs2273328,rs202487,rs73069315,rs59238791,rs202486,rs143370970,rs202485,rs202484,rs548669,rs548666,rs202483,rs202482,rs202481,rs202480,rs202479,rs202478,rs5839915,20:1680109_GA_G,rs74764447,rs149538314,rs8124886,rs6043867,rs757096628,rs73069323,rs73069324,rs73069326,rs73069328,rs6111024,rs1438101,rs1438099,rs1438098,rs6111031,rs56198409,rs6111034,rs79020972,rs202477,rs11480839,rs6111041,rs6111046,20:1683531_TAGAC_T,rs202475,rs138813509,rs138436223,rs143719429,rs202474,rs6105571,rs6105572,rs12480227,rs11906069,rs202473,rs6105573,rs6105574,rs6111054,rs141638318,rs112799593,rs6111055,rs1595460,rs1595459,rs202472,rs139258893,rs6111058,rs692847,rs659348,rs692891,rs497134,rs648051,rs148617538,rs2247997,20:1686447_CT_C,rs6105579,rs528256,rs6135734,rs532097,rs629878,rs629867,rs6135736,rs629812,rs480743,rs480769,rs482716,20:1687979_TA_T,rs117047021,20:1688002_AAAAC_A,rs2164863,rs615007,rs59719880,rs513119,rs7266946,rs513201,rs514813,rs514965,rs515883,rs515973,rs601652,rs6111080,rs517578,rs600303,rs599837,rs5839916,rs2317999,rs2318000,rs2318001,rs2318002,rs543245,rs588555,rs544110,rs735994,rs735993,rs743082,rs549733,rs585046,rs200393694,rs572620,rs577238,rs578127,rs679310,rs578970,rs493466,rs73071317,rs11699191,rs498139,rs118020298,rs75794831,rs202493,20:1694183_CT_C,rs4633984,rs202494,rs202495,rs202496,rs202497,rs202498,rs202499,rs202500,rs202501,rs202502,rs4814498,rs202503,rs59027828,rs112882613,rs111553886,rs202504,rs202505,rs7509334,rs139268472,rs11480959"
input="rs540538026,rs367896724 ,rs201106462"
input="rs2234962,rs28581202,6:122113614_CT_C,rs55851300,rs1763604"
input="rs55851300, rs1051375, rs34081637,rs12541595 , rs10774625, rs10774625,rs806322,rs4765663"
input="rs776293589"
input="rs115286434"
input = "TTN"
input = "KCNH2"
input = "SCN10A"
input ="rs55851300" # kcnh2
input ="rs34081637" # CACNA1G
input = "rs11067264"
input = "BCAT1"
input = "20:62000000-62300000"
input = "rs104571" #wrong
input="rs2234962,rs1763604 "
input="20:40000000-50000000"
query <- process_user_input(input,mapping.proteincoding,standalone=T)
tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)
data <- extract_multiple_variants(tabix_query,
                                  f.data_p=paste0(datadir,"/tophits_data/unadjusted.logP.outfile.tsv.gz.tophits.gz"),
                                  f.data_beta=paste0(datadir,"/tophits_data/unadjusted.BETA.outfile.tsv.gz.tophits.gz"),
                                  f.data_se=paste0(datadir,"/tophits_data/unadjusted.SE.outfile.tsv.gz.tophits.gz"),
                                  f.data.index=paste0(datadir,"/tophits_data/unadjusted.logP.outfile.index.tsv.gz.tophits.gz")
                                  )
# input="rs10399793,rs530867301"
# query <- process_user_input(input,mapping.proteincoding)
# tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)
# 
# data <- extract_multiple_variants(tabix_query,
#                                   f.data_p=paste0(datadir,"/full_data_combined/unadjusted.logP.outfile.tsv.gz"),
#                                   f.data_beta="",
#                                   f.data_se="",
#                                   f.data.index=paste0(datadir,"/full_data_combined/unadjusted.logP.outfile.tsv.gz")
# )

i=data$max.snp.index[1]
ecg_plot <- make_ecg_plot(vct_snp_p=data$df_snp_p[i,],
                          vct_snp_beta=data$df_snp_beta[i,],
                          vct_snp_se=data$df_snp_se[i,],
                          vct_snp_info=data$df_snp_info[i,],
                          df_ecg_stats=df_ecg_unadjusted,plot_adjusted_means = T,
                          invert=FALSE)
ecg_plot


dftable1 <- fread("/Users/niek/Dropbox/Gwasshared_2/ECG_morphology/table1.tsv")

data$df_snp_info[data$df_snp_info$SNP %in% "rs57926290",]
regional_plot(data,LD =T)


fwrite(
  data.frame(
           uniqid=data$df_snp_info$uniqid,
          SNP=data$df_snp_info$SNP,
          CHR=data$df_snp_info$CHR,
          BP=data$df_snp_info$BP,
          GENPOS=0,
          ALLELE1=data$df_snp_info$EFAL,
          ALLELE0=data$df_snp_info$NEFAL,
          A1FREQ=data$df_snp_info$EAF,
          INFO=data$df_snp_info$INFO,
          BETA=data$df_snp_beta[,data$max.snp.index[2]],
          SE=data$df_snp_se[,data$max.snp.index[2]],
          pval=10^-abs((data$df_snp_p[,data$max.snp.index[2]]))
),
file="/Users/niek/Dropbox/Gwasshared_2/ECG_morphology/st-depression/sumstats.KCNB1.ECG.142ms.tsv",quote = F,row.names = F,sep="\t")

for (v in  unique(dftable1$SNP)) {
   input=v
   query <- process_user_input(input,mapping.proteincoding)
   tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)
   data <- extract_multiple_variants(tabix_query,
                                     f.data_p="unadjusted.logP.outfile.tsv.chr.gz.tophits.gz",
                                     f.data_beta="",
                                     f.data_se="",
                                     f.data.index="unadjusted.logP.outfile.index.tsv.gz.tophits.gz"
   )
   pdf(paste0("/Users/niek/Dropbox/Gwasshared_2/ECG_morphology/Figures/allvariants_ecgenetics/",unique(dftable1[dftable1$SNP %in% v,]$`Nearest gene (or any within 10kb)`),".",v,".pdf"))
   ecg_plot <- make_ecg_plot(vct_snp_p=data$df_snp_p[1,],
                             vct_snp_beta=data$df_snp_beta,
                             vct_snp_se=data$df_snp_se,
                             vct_snp_info=data$df_snp_info[1,],
                             df_ecg_stats=df_ecg_unadjusted,
                             invert=FALSE)   
   print(ecg_plot)
   dev.off()
}



make_ecg_plot <- function(vct_snp_p,vct_snp_beta="",vct_snp_se="",vct_snp_info,df_ecg_stats=df_ecg_unadjusted,invert=FALSE){
   #print(nrow(vct_snp_p))
   dfecg <- cbind(as.data.frame(df_ecg_stats),
                  t(vct_snp_p)
   )
   names(dfecg)[5:ncol(dfecg)] <- vct_snp_info$SNP
   dfecg = melt(dfecg, id.vars=names(dfecg)[1:4])
   names(dfecg)[5:6] <- c("SNP","logP")
   dfecg$logP <- as.numeric(dfecg$logP)
   if(invert==TRUE) dfecg$logP <- -1*dfecg$logP
   
   if(is.data.frame(vct_snp_beta) & is.data.frame(vct_snp_se) & length(vct_snp_p)==500 ) {
      dfecg$BETA <- as.numeric(vct_snp_beta)
      if(invert==TRUE) dfecg$BETA <- -1*dfecg$BETA
      dfecg$cil <-dfecg$BETA - 1.96*as.numeric(vct_snp_se)
      dfecg$ciu <-dfecg$BETA + 1.96*as.numeric(vct_snp_se)
      
      # 
      BETAf=max(abs(c(dfecg$cil,  dfecg$ciu)) )/max(abs(as.numeric(dfecg$logP)),na.rm=TRUE)
      dfecg$BETA  <- dfecg$BETA /BETAf
      dfecg$cil  <- dfecg$cil /BETAf
      dfecg$ciu  <- dfecg$ciu /BETAf
      
      # 
   }
   
   # dfecg <- cbind(as.data.frame(df_ecg_stats),
   #                vct_snp_p=as.numeric(as.character(unlist(vct_snp_p))))
   
   dfecg$time <- 2*(1:500-250)
   dfecg$SNP <- unname(unlist(as.data.frame(sapply( vct_snp_info$SNP, function(x) rep(x,500)))))
   
   Pf=max(abs(dfecg$mean0))/max(abs(as.numeric(dfecg$logP)),na.rm=TRUE)
   Pfdefault=max(abs(dfecg$mean0))/8
   Pf = min(Pfdefault,Pf)
   
   
   #print(Pf)
   #print(Pfdefault)
   dfecg$mean0  <- dfecg$mean0 /Pf
   #print(dfecg[order(dfecg$time),])
   #### FOR MULTIPLE LINES IN ONE PLOT, DONT KNOW YET HOW THIS IS DONE BEST WITH USER INPUT. 
   ecgplot <- ggplot()+
     geom_line(data=dfecg, aes(x = time, y= logP ,group=SNP,colour = SNP,
                               text = paste0('Time(ms): ',time,
                                             '<br>P-value: ',formatC(10^(-abs(logP)),digits=3),
                                             '<br>-Log10(P-value): ',abs(logP) ) ))  +
     geom_line(data=dfecg,aes(x = time,y=mean0),linetype="dotted",colour="black") +
         scale_color_manual(name="Groups",values=c("red", "blue","#00AFBB", "#E7B800", "#FC4E07"))+
         labs(title=paste(vct_snp_info$SNP,collapse = ", " ), subtitle=vct_snp_info$Gene)
   # 
   # # 
   # ecgplot <- ggplot(dfecg,aes(x=time))+
   #    geom_line(data=dfecg,colour="red", aes(x = time, y= logP ,
   #                                           text = paste0('Time(ms): ',time,
   #                                                         '<br>P-value: ',formatC(10^(-abs(logP)),digits=3),
   #                                                         '<br>-Log10(P-value): ',abs(logP) ),
   #                                           group=1  ))  +
   #    geom_line(data=dfecg,aes(x = time,y=mean0),linetype="dotted",colour="black") +
   #    scale_color_manual(name="Groups",values=c("red", "blue"))+
   #    ylab("Signed -log10(P-value)") + xlab("Time (ms)") +
   #    labs(title=vct_snp_info$SNP, subtitle=vct_snp_info$Gene)
   
   if(is.data.frame(vct_snp_beta) & is.data.frame(vct_snp_se) & length(vct_snp_p)==500 ) {
      ecgplot <- ecgplot +  geom_ribbon(aes(ymin=cil, ymax=ciu,x = time), linetype=2, alpha=0.1) + 
         geom_line(colour="grey",aes(y=BETA)) + 
         scale_y_continuous(
            
            sec.axis = sec_axis(~ . *BETAf, name = "Effect (per standard deviation) and CI")
         )
   }
   
   #ecgplot
   
   return(ecgplot) 
}

data$df_snp_p[1,] <- -1*data$df_snp_p[1,]

make_ecg_plot(vct_snp_p=data$df_snp_p,
                          vct_snp_beta=data$df_snp_beta,
                          vct_snp_se=data$df_snp_se,
                          vct_snp_info=data$df_snp_info,
                          df_ecg_stats=df_ecg_unadjusted,
                          invert=FALSE)



hm <- make_heatmap_plot(df_snp_p = data$df_snp_p,
                        vct_snp_info = data$df_snp_info,
                  df_ecg_unadjusted, pvaltreshold = 0.00000001
                  )

str(hm$x )

rp <- regional_plot(data,for_plotly = F,LD=FALSE)
rp
rply <- regional_plot(data,for_plotly = T,LD=FALSE)
p1 <- ggplotly(rply$plt_regional,source="source_rp")
attr(p1,"source")<-"asd"
p2 <- ggplotly(rply$plt_genes,source="Asd")
subplot(p1, p2, nrows = 2, margin = 0.09, heights = c(0.8, 0.2))

attr(ggplotly(rply$plt_regional,source="asd"),"source")

 ensure_one <- function(plots, attr) {
   attrs <- lapply(plots, "[", attr)
   for (i in seq_along(attrs)) {
     if (!identical(attrs[[1]], attrs[[i]])) {
       warning("Can only have one: ", attr, call. = FALSE)
       break
     }
   }
   attrs[[length(attrs)]][[1]]
 }
 
pdf(paste0("/Users/niekverw/Downloads/ecg_plot",input,".pdf"),height = 5,width=6)
ecg_plot
dev.off()
paste(unname(data$df_snp_p[1,]),collapse=',' )

dput(data$df_snp_p[100,])
#df.static.rsid , df.static.pos






# rs17162854
# rs2503715
# rs60965469
# rs846111
# rs578347
# rs113712703
# rs1763604
# rs2298632
# rs61775418
# rs116015634
# rs1707334
# rs112288630
# rs590041
# rs985783
# rs34370185
# rs2474370
# rs72677052
# rs6702939
# rs143221246
# rs72681760
# rs584096
# rs12142815
# rs776293589
# rs115286434
# rs4073778
# rs12034347
# rs553620963
# rs1870940
# rs4745
# rs10908504
# rs10797060
# rs12143842
# rs59852339
# rs116381680
# rs10918192
# rs11580561
# rs531706
# rs12569204
# rs1200130
# rs67779898
# rs1591734
# rs57831591
# rs76120745
# rs575784488
# rs1104859
# rs4951365
# rs863724
# rs35232605
# rs568925773
# rs340849
# rs1532501
# rs115282653
# 1:225664878_AAGAG_A
# rs35937215
# rs12117603
# rs10193517
# rs6730558
# rs1047392
# rs13020178
# rs1368888
# rs56036545
# rs13026826
# rs58470117
# rs2121266
# rs4671222
# rs7609203
# rs759664
# rs11677371
# 2:84085386_AGT_A
# rs2732861
# rs17407658
# rs16860332
# rs115142950
# rs16866397
# rs955738
# rs11681740
# rs6735214
# rs7605146
# 2:213211032_GC_G
# rs907683
# rs35394392
# rs10165767
# rs2595508
# rs141501074
# rs142536918
# rs7650482
# rs4685090
# rs13092272
# rs760390566
# rs55644439
# rs6772646
# rs62241772
# rs41312411
# rs572474770
# rs571037730
# rs76521806
# rs6801957
# rs11921457
# rs1986599
# rs113855341
# rs62256002
# rs6774470
# rs332371
# rs114564188
# rs4677137
# rs766206762
# rs72620089
# 3:111574251_TCTTC_T
# rs1254386
# rs1403767
# rs4683632
# rs11923657
# rs11928162
# rs759623919
# rs6441207
# rs55853354
# rs4894803
# rs545261249
# rs375818588
# rs1706003
# rs2292804
# rs768695
# rs199745043
# rs60063444
# rs1906599
# rs6533694
# rs139446049
# rs61760320
# rs6846721
# 4:155350842_GT_G
# rs13127767
# rs2062591
# rs6824178
# rs76383255
# rs2333469
# rs260399
# rs3828591
# rs111869282
# rs7714577
# rs114419215
# rs111778721
# 5:82594418_CT_C
# rs35658696
# rs890924
# rs335214
# 5:132428677_ATT_A
# rs1476096
# rs11748963
# 5:145398089_ACAC_A
# rs34289553
# rs10054375
# rs10076602
# rs763696876
# rs11745586
# rs12652726
# rs61367953
# rs9502009
# rs72825038
# rs775204216
# rs547691407
# rs9396858
# rs78144291
# rs13211463
# rs9467373
# rs301395
# rs772081778
# rs3130582
# rs730506
# rs533119559
# rs4478445
# rs28763977
# rs59365541
# rs210977
# rs663474
# rs376054711
# rs11153730
# rs12212870
# rs28581202
# 6:122113614_CT_C
# rs12209641
# rs6940945
# rs10457469
# rs11961815
# rs378273
# rs573680085
# rs77170649
# rs6570952
# rs13208950
# rs6914036
# rs10464373
# rs6461703
# 7:35281105_TG_T
# rs11764098
# rs343057
# rs372543528
# rs74524365
# rs11771799
# rs13245899
# rs2520258
# rs3128384
# rs1997571
# 7:116890450_AT_A
# rs56216811
# rs35680304
# rs117104597
# rs55851300
# rs143016484
# rs4517166
# 8:5828558_GGA_G
# rs7820146
# rs73195040
# rs61494024
# rs62505473
# rs7833518
# rs571690420
# rs11776287
# rs1562986
# rs2445704
# rs72671655
# rs147543128
# rs1461990
# rs16898508
# rs4871397
# rs72711245
# rs34866937
# rs13264524
# rs11784619
# 9:33272391_CA_C
# rs770250170
# rs748802
# rs1415221
# rs34382015
# rs10119967
# rs74552897
# rs7866966
# rs7030930
# rs10980770
# rs7850041
# rs12001868
# rs11257655
# rs72786007
# rs788076
# rs12779954
# rs211315
# rs12777211
# rs1733724
# rs10740131
# rs1037985
# 10:73832438_CCT_C
# rs2177843
# rs6480713
# rs80341670
# rs1658323
# rs1769758
# rs3802662
# rs117524646
# rs10736059
# rs117929802
# rs11187838
# rs10882661
# rs11189038
# rs72845807
# rs12257568
# rs1977821
# rs10787270
# rs770537744
# rs7918405
# rs740746
# rs2234962
# rs1078806
# rs2074238
# rs117236107
# rs2625322
# rs7118294
# rs3729989
# rs7479922
# 11:60533354_CA_C
# rs174536
# rs36085849
# rs638138
# rs777300465
# 11:75923599_CA_C
# rs10898646
# 11:95034279_TA_T
# rs671333
# rs4935884
# rs11220270
# rs11221500
# 12:1633920_TGGG_T
# rs11062103
# rs11609729
# rs1051375
# rs10849144
# rs7307613
# rs34219605
# rs11049136
# 12:31620079_AG_A
# rs7315695
# rs148885520
# 12:51197281_TA_T
# rs113671070
# 12:57156410_AAAAC_A
# rs1149003
# 12:110625059_GTTGTTT_G
# rs141066394
# rs722726
# rs883079
# rs61931005
# rs188348223
# rs11067264
# rs7132327
# rs142576422
# rs373465055
# rs7317804
# rs11616720
# rs7993628
# 13:49534769_CT_C
# rs1965800
# rs1886512
# rs2769884
# rs111403728
# rs6571769
# rs365990
# 14:23929566_AT_A
# rs562161427
# rs2108706
# rs17180489
# rs12880339
# rs11626401
# 14:96872394_AT_A
# rs34509141
# rs78496852
# rs3131560
# rs74018128
# rs7496259
# rs2019185
# rs3825977
# rs4777245
# rs7164883
# 15:75187470_CT_C
# rs200201338
# rs12595786
# rs28513927
# 15:94705104_GTGTA_G
# rs4966020
# rs1045716
# rs6600238
# rs12924536
# rs735052
# rs2080512
# rs1344533
# rs30226
# rs3915425
# rs12448482
# rs4238961
# rs143721880
# rs4784939
# rs150728296
# rs62044987
# rs72791317
# rs35559328
# rs6565060
# 16:88527325_AGGGGC_A
# rs4614769
# rs12945974
# 17:12642266_TTA_T
# rs4925095
# rs8077048
# 17:28267098_CT_C
# rs1088450
# rs772490616
# rs17608766
# rs34081637
# rs138998340
# rs9909004
# rs8067476
# rs236586
# rs312755
# rs28360848
# rs12051698
# rs766320
# rs35236676
# rs12454640
# rs11665344
# rs356907
# rs476348
# rs2848901
# rs1126225
# rs9952412
# rs689049
# rs200181073
# rs113394178
# rs10853935
# rs10409940
# rs9630873
# rs35828463
# rs191904180
# rs12610825
# rs6084574
# rs11699612
# rs369307656
# rs2025096
# rs4811601
# rs59499156
# rs8115156
# rs310614
# rs3787112
# rs2830965
# rs762834393
# rs1999323
# rs41312993
# rs1805128
# rs2186370
# rs133885
