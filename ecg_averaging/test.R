#install_github("vankesteren/rpeaks")
library(rpeaks)
source("/Users/niek/repos/itc-main-repo/Niek/ECG/exerciseecg_data/ECG_averaging/E-ECG_saecg.functions.R")

fxml <- "/Users/niek/data/ecg_exercise/1000131_6025_0_0.xml"
samplingrate=500


dfxml <- getECGsignal(fxml)
vctecg <- rowSums(dfxml)
y=vctecg[30000:39000]
plot(y,pch=".")

# REMOVE LOWER FREQUENCIES FROM SIGNAL
# APPLY FOURIER TRANSFORM
fresult = fft(as.numeric(unlist(y)))
# Once in frequency domain, revert the first n values to 0
fresult[1 : round(length(fresult)*5/samplingrate)]=0;
# Also revert the last n values to 0
fresult[(length(fresult) - round(length(fresult)*5/samplingrate)):length(fresult)]=0;
# Translate signal back to ECG data using inverse Fourier transform
y_r=Re(ifft(fresult));

#####


ir<- rpeaks_pan_tompkins(y,samplingrate)*samplingrate
plot(y,type="l") + abline (v=ir,col="red")


window.finetune=25 #(+/- window within the total-window (defined above) that will be used to finetune). 

dfy=saecg_unadjusted(QRSpeaks=ir,VctECG=y_r,tsc=FALSE,samplerate=500,shift_ms=-125,window_ms=250)
dfy_outliers <- removeoutliers(dfy,return="Outlierindex",omethod="pearson",iterative=TRUE)

find.max <- function(target,ref,istart,iend){
  c <- unlist(lapply(-(istart-1):(istart),function(i) cor(target[(istart+i):(iend+i)],ref[istart:iend])))
  return(which.max(c))
}
i.template <- which.max(rowSums(cor(dfy))) # define template according to the maxium sum of R. 
ir_new <- c()
istart=round((nrow(dfy)/2))-window.finetune
iend=round((nrow(dfy)/2))+window.finetune
for (i in 1:ncol(dfy)){
  i.max <- find.max(target=dfy[,i],ref=dfy[,i.template],istart,iend)
  ir_new <- c(ir_new,(ir[i]-(istart-i.max) ))

}

#plot(dfy[,i.template]) #
shift=median(1:nrow(dfy))-which.max(dfy[,i.template]) # shift template to the middle and adjust all R's. 
ir_new <- ir_new - shift # <-- 
plot(y[0:1000],type="l") + abline (v=ir_new,col="red")



#### CHECK. 
dfy_new=saecg_unadjusted(QRSpeaks=ir_new,VctECG=y_r,tsc=FALSE,samplerate=500,shift_ms=-125,window_ms=250)
dfy_outliers <- removeoutliers(dfy_new,return="Outlierindex",omethod="pearson",iterative=TRUE) # no outliers anymore. 
plot(dfy_new[,2])
pheatmap::pheatmap(cor(dfy_new))
pheatmap::pheatmap(as.matrix(dfy_new),cluster_cols = T,cluster_rows = F)



