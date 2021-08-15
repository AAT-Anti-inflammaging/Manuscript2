library("tidyverse")

dat <- read.csv("fly.csv", header = T)
dat <- dat%>%mutate(GB=paste("M:",F,"x","F:",M),
                    hAAT=hAAT2*0.4)%>%mutate(AnimalhAAT=hAAT/No,WeighthAAT=hAAT/Weight)

dat$GB <- factor(dat$GB, levels = c("M: Background x F: Lsp2-Gal4", "M: UAS-hAAT x F: Lsp2-Gal4",
                                    "M: UAS-D6-hAAT (G257) x F: Lsp2-Gal4","M: Background x F: Yolk-Gal4",
                                    "M: UAS-hAAT x F: Yolk-Gal4"))
dat$Sex <- factor(dat$Sex, levels = c("Male","Female","-"))

dat%>%filter(F=="Background")%>%summarise(max(AnimalhAAT),max(WeighthAAT)) # LLOQ


dat2 <- dat%>%mutate(cnt=1)%>%group_by(Note,GB,Sex)%>%
  summarise(hAAT_animal=mean(AnimalhAAT),sd_animal=sd(AnimalhAAT),
            hAAT_weight=mean(WeighthAAT),sd_weight=sd(WeighthAAT),
            N=sum(cnt))%>%mutate(se_animal=sd_animal/sqrt(N),se_weight=sd_weight/sqrt(N))%>%
  mutate(ymax_animal=hAAT_animal+se_animal, ymax_weight=hAAT_weight+se_weight)

ggplot(dat2%>%filter(Note=="Larvae"))+aes(x=GB, y=hAAT_animal, color=GB)+
  geom_bar(stat = "identity", fill=NA, width=0.5, size=1)+
  geom_hline(yintercept = 2.55755, linetype="dashed", size=1)+
  scale_color_manual(values = c("orange","dodgerblue1","dodgerblue4","orange3","dodgerblue3"))+
  scale_y_continuous(limits = c(0,30),breaks = seq(0,30,5))+
  labs(x="", y="Human AAT (ng/larvae)",title="Larvae")+
  geom_errorbar(aes(ymin=hAAT_animal-se_animal, ymax=hAAT_animal+se_animal), width=0.2,size=1)+
  theme_bw()+theme(text=element_text(size=18), axis.text.x = element_text(angle=30,hjust = 1),
                   legend.position = "")+theme(axis.text.x = element_blank())
ggsave("larvae_animal.png",width=15,height = 12,units = "cm",dpi = 600)

ggplot(dat2%>%filter(Note=="Adults"))+aes(x=GB, y=hAAT_animal, color=GB)+
  geom_bar(stat = "identity", fill=NA, width=0.5, size=1)+
  geom_hline(yintercept = 2.55755, linetype="dashed", size=1)+
  facet_grid(~Sex,scales = "free_x")+
  scale_y_continuous(limits = c(0,30),breaks = seq(0,30,5))+
  scale_color_manual(values = c("orange","dodgerblue1","dodgerblue4","orange3","dodgerblue3"))+
  labs(x="", y="Human AAT (ng/fly)",title = "Adult")+
  geom_errorbar(aes(ymin=hAAT_animal-se_animal, ymax=hAAT_animal+se_animal), width=0.2,size=1)+
  theme_bw()+theme(text=element_text(size=18), axis.text.x = element_text(angle=30,hjust = 1),
                   legend.position = "")+theme(axis.text.x = element_blank())
ggsave("adult_animal.png",width=15,height = 12,units = "cm",dpi = 600)

ggplot(dat2)+aes(x=GB, y=hAAT_animal, color=GB)+
  geom_bar(stat = "identity", fill=NA, width=0.5, size=1)+
  geom_hline(yintercept = 2.55755, linetype="dashed", size=1)+
  facet_grid(~Note+Sex,scales = "free_x")+
  scale_y_continuous(limits = c(0,30),breaks = seq(0,30,5))+
  scale_color_manual(values = c("orange","dodgerblue1","dodgerblue4","orange3","dodgerblue3"))+
  labs(x="", y="Human AAT [ng/fly]",title = "Human AAT Level")+
  geom_errorbar(aes(ymin=hAAT_animal-se_animal, ymax=hAAT_animal+se_animal), width=0.2,size=1)+
  theme_bw()+theme(text=element_text(size=18), axis.text.x = element_text(angle=30,hjust = 1),
                   legend.position = "")+theme(axis.text.x = element_blank())
ggsave("animal.png",width=18,height = 12,units = "cm",dpi = 600)



ggplot(dat2%>%filter(Note=="Larvae"))+aes(x=GB, y=hAAT_weight, color=GB)+
  geom_bar(stat = "identity", fill=NA, width=0.5, size=1)+
  geom_hline(yintercept = 3.499304, linetype="dashed", size=1)+
  scale_y_continuous(limits = c(0,30),breaks = seq(0,30,5))+
  scale_color_manual(values = c("orange","dodgerblue1","dodgerblue4","orange3","dodgerblue3"))+
  labs(x="", y="Human AAT (ng/mg)",title="Larvae")+
  geom_errorbar(aes(ymin=hAAT_weight-se_weight, ymax=hAAT_weight+se_weight), width=0.2,size=1)+
  theme_bw()+theme(text=element_text(size=18), axis.text.x = element_text(angle=30,hjust = 1),
                   legend.position = "")+theme(axis.text.x = element_blank())
ggsave("larvae_weight.png",width=15,height = 12,units = "cm",dpi = 600)

ggplot(dat2%>%filter(Note=="Adults"))+aes(x=GB, y=hAAT_weight, color=GB)+
  geom_bar(stat = "identity", fill=NA, width=0.5, size=1)+
  geom_hline(yintercept = 3.499304, linetype="dashed", size=1)+
  facet_grid(~Sex,scales = "free_x")+
  scale_y_continuous(limits = c(0,30),breaks = seq(0,30,5))+
  scale_color_manual(values = c("orange","dodgerblue1","dodgerblue4","orange3","dodgerblue3"))+
  labs(x="", y="Human AAT (ng/mg)", title = "Adult")+
  geom_errorbar(aes(ymin=hAAT_weight-se_weight, ymax=hAAT_weight+se_weight), width=0.2,size=1)+
  theme_bw()+theme(text=element_text(size=18), axis.text.x = element_text(angle=30,hjust = 1),
                   legend.position = "")+theme(axis.text.x = element_blank())
ggsave("adult_weight.png",width=15,height = 12,units = "cm",dpi = 600)


ggplot(dat2)+aes(x=GB, y=hAAT_weight, color=GB)+
  geom_bar(stat = "identity", fill=NA, width=0.5, size=1)+
  geom_hline(yintercept = 3.499304, linetype="dashed", size=1)+
  facet_grid(~Note+Sex,scales = "free_x")+
  scale_y_continuous(limits = c(0,30),breaks = seq(0,30,5))+
  scale_color_manual(values = c("orange","dodgerblue1","dodgerblue4","orange3","dodgerblue3"))+
  labs(x="", y="Human AAT [ng/mg]", title = "Human AAT Level")+
  geom_errorbar(aes(ymin=hAAT_weight-se_weight, ymax=hAAT_weight+se_weight), width=0.2,size=1)+
  theme_bw()+theme(text=element_text(size=18), axis.text.x = element_text(angle=30,hjust = 1),
                   legend.position = "")+theme(axis.text.x = element_blank())
ggsave("weight.png",width=18,height = 12,units = "cm",dpi = 600)
