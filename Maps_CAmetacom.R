# Create maps of the sites used per year:


library(ggmap)
library(OpenStreetMap)

CA <- get_map(location=c(lon=-120, lat=37), zoom=6, maptype="hybrid",
              source="google")

# Create a wide shot of CA bay area
CAmap <- ggmap(CA, extent="device") +
  labs(x="",y="")+
  theme(axis.text=element_blank(), axis.ticks=element_blank())

quartz(height=4, width=4)
print(CAmap)

setwd("~/Documents/Thesis Research/CA metacoms/Manuscript/Figures/Maps")
png("CAmap.png", height=600, width=600, res=200)
print(CAmap)
dev.off()

# 2009:
CA <- get_map(location=c(lon=-121.8, lat=37.6), zoom=9, maptype="satellite",
              source="google")

map_2009 <- ggmap(CA, extent="device") +
  geom_point(aes(x=Long, y=Lat), data=data.frame(Xcov_2009),
             color="red", size=1.5)+
  labs(x="",y="")+
  theme(axis.text=element_blank(), axis.ticks=element_blank())
            
quartz(height=4, width=4)
print(map_2009)

setwd("~/Documents/Thesis Research/CA metacoms/Manuscript/Figures/Maps")
png("2009map.png", height=600, width=600, res=200)
print(map_2009)
dev.off()

# 2010:

map_2010 <- ggmap(CA, extent="device") +
  geom_point(aes(x=Long, y=Lat), data=data.frame(Xcov_2010),
             color="red", size=1.5)+
  labs(x="",y="")+
  theme(axis.text=element_blank(), axis.ticks=element_blank())

quartz(height=4, width=4)
print(map_2010)

png("2010map.png", height=600, width=600, res=200)
print(map_2010)
dev.off()

# 2011:

map_2011 <- ggmap(CA, extent="device") +
  geom_point(aes(x=Long, y=Lat), data=data.frame(Xcov_2011),
             color="red", size=1.5)+
  labs(x="",y="")+
  theme(axis.text=element_blank(), axis.ticks=element_blank())

quartz(height=4, width=4)
print(map_2011)

png("2011map.png", height=600, width=600, res=200)
print(map_2011)
dev.off()

# 2012:

map_2012 <- ggmap(CA, extent="device") +
  geom_point(aes(x=Long, y=Lat), data=data.frame(Xcov_2012),
             color="red", size=1.5)+
  labs(x="",y="")+
  theme(axis.text=element_blank(), axis.ticks=element_blank())

quartz(height=4, width=4)
print(map_2012)

png("2012map.png", height=600, width=600, res=200)
print(map_2012)
dev.off()