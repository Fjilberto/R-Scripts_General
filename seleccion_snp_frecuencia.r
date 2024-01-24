library(dartR)
library(hierfstat)
library(stringr)

#### Los datos DEBEN estar en formato genind!!!

basic_datos = basic.stats(datos, diploid = TRUE)
freq<-basic_datos$pop.freq
freq2<-as.data.frame(freq)
freq2<-t(freq2)
freq2 <- cbind(rownames(freq2), data.frame(freq2, row.names=NULL))
freq3 = freq2[seq(3, nrow(freq2), 3),]

# En las siguientes lineas cambiar los nombres de columnas según corresponda a los grupos de interes.
# Los grupos en este ejemplo son F y ;

colnames(freq3)<-c("snp","F1","F2","M1","M2")
freq3$F1<-as.numeric(freq3$F1)
freq3$F2<-as.numeric(freq3$F2)
freq3$M1<-as.numeric(freq3$M1)
freq3$M2<-as.numeric(freq3$M2)

freq3$diff1<-abs(freq3$M1-freq3$F1)
freq3$diff2<-abs(freq3$M2-freq3$F2)

freq3$sumdiff<-freq3$diff1+freq3$diff2
freq3$sum<-freq3$F1+freq3$M1+freq3$F2+freq3$M2

#### Selección de SNPs contrastantes entre grupos

freq3<-subset(freq3, sum!=4) # Eliminar SNP monomorficos
seleccion<-subset(freq3, sumdiff==2) # Seleccionar SNP contrastantes
seleccion$snp<-sub(".", "", seleccion$snp)
seleccion$snp<-str_c(str_sub(seleccion$snp,1,-6))

write.table(seleccion,file = "SNP_contrastantes.txt")

#### Selección de SNPs fijos en un grupo diversos en el otro

# Seleccionar F 1/0

seleccion2f<-subset(freq3,freq3$F1==1|freq3$F1==0)
seleccion2f<-subset(seleccion2f,0.3<seleccion2f$M1 & seleccion2f$M1<0.7)

# Seleccionar M 1/0

seleccion2m<-subset(freq3,freq3$M1==1|freq3$M1==0)
seleccion2m<-subset(seleccion2m,0.3<seleccion2m$F1 & seleccion2m$F1<0.7)

# Unir dataframes

seleccion_fijodiverso <- merge(seleccion2f, seleccion2m, 
                       by = c("snp","F1","F2","M1","M2","diff1",
                              "diff2","sumdiff","sum"), all = TRUE)

seleccion_fijodiverso$snp<-sub(".", "", seleccion_fijodiverso$snp)
seleccion_fijodiverso$snp<-str_c(str_sub(seleccion_fijodiverso$snp,1,-6))

write.table(seleccionsexo,file = "SNP_fijodiverso.txt")
