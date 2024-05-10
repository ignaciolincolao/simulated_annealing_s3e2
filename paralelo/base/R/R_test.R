#install.packages(c("ggplot2", "sf", "ggmap"))
library(ggplot2)

samples <- read.csv("./samples.dat", sep=" ", dec=".", header = FALSE)
observations <- read.csv("observations.dat", sep=" ", dec=".", header = FALSE)

observations$V4 <- observations$V3*(1/10000000) + observations$V2
observations[which.min(observations$V4),]

valRange <- function(a,b, c){
  return (a + c * ( b - a ));
}

samples$tempInit <-valRange(100,10000,samples$V2) 
samples$coolingRate <-valRange(0.9,0.999,samples$V3) 
samples$len1 <-valRange(1,100,samples$V4) 
samples$len2 <-valRange(1,100,samples$V5) 
samples$n_block <-32* as.integer(valRange(1,32,samples$V6))
samples$n_thread <-32*as.integer(valRange(1,32,samples$V7))
samples$id <- 1:nrow(samples)
observations$id <- 1:nrow(observations)
observations <- setNames(observations, c("V1","z","count","lz","id"))

data <- merge(samples, observations, by="id")


observations[which.min(data$lz),]

# data$tempInit<-factor(data$tempInit)
# data$coolingRate<-factor(data$coolingRate)
# data$len1<-factor(data$len1)
# data$len2<-factor(data$len2)
# data$n_block<-factor(data$n_block)
# data$n_thread<-factor(data$n_thread)

data$end_n_block <- c(data$n_block[-1], NA)
data$end_n_thread <- c(data$n_thread[-1], NA)

ggplot(data, aes(x=n_block, y=n_thread)) +
  geom_segment(aes(xend = end_n_block, yend = end_n_thread, color = z,alpha=id/max(id)),
               arrow = arrow(type = "closed"))


ggplot(data, aes(x=n_block, y=n_thread, color=log(z))) +
  geom_point(size=4) +
  scale_color_gradient(low="blue", high="red") +
  labs(title="Scatter plot with 3 variables",
       x="n_block",
       y="n_thread",
       color="Z") +
  theme_minimal()


modelo_mult <- lm(lz ~ tempInit + coolingRate + len1 + len2 + n_block + n_thread, data=data)
summary(modelo_mult) 
par(mfrow=c(2,2))  
plot(modelo_mult)




log(0.05)
log(0.03)



p_i = p_max - (p_max -p_0) * exp(-k * i)



# Cargar la librería ggplot2
library(ggplot2)

# Definir los parámetros
p_max <- 0.8     # Probabilidad máxima objetivo
p_0 <- 0.01     # Probabilidad inicial
k <- 1e-5     # Parámetro de ajuste de la convergencia
iteraciones <- 1:200000  # Rango de iteraciones

# Calcular los valores de p_i para cada iteración
p_i <- p_max - (p_max - p_0) * exp(-k * iteraciones)

# Crear un data frame para ggplot
datos <- data.frame(Iteracion = iteraciones, Probabilidad = p_i)

# Graficar
ggplot(datos, aes(x = Iteracion, y = Probabilidad)) +
  geom_line() + 
  theme_minimal() +
  ggtitle("Convergencia de la Probabilidad") +
  xlab("Iteración") +
  ylab("Probabilidad p_i")

p_i[50000]
# Definir los parámetros



# Cargar la librería ggplot2
library(ggplot2)
p <- 0.4
n <- 10
iteraciones <- 1:n
# Calcular los valores de p_i para cada iteración
p_i <- p * (1-p)**(iteraciones-1)
p_i[length(p_i)]
promedio <-sum(p_i)
p_i <- p_i/promedio
p_i
for (i in 2:n) {
  p_i[i] <- p_i[i]+p_i[i-1] 
}
p_i[n] = 1
# Crear un data frame para ggplot
datos <- data.frame(Iteracion = iteraciones, Probabilidad = p_i)
# Graficar
ggplot(datos, aes(x = Iteracion, y = Probabilidad)) +
  geom_line() + 
  theme_minimal() +
  ggtitle("Convergencia de la Probabilidad") +
  xlab("Iteración") +
  ylab("Probabilidad p_i")
p_i



senoPersonalizado <- function(x, minVal, maxVal) {
  # Amplitud es la mitad de la diferencia entre el valor máximo y mínimo
  amplitud <- (maxVal - minVal) / 2
  
  # El centro es el punto medio entre el máximo y el mínimo
  centro <- (maxVal + minVal) / 2
  
  # Ajustar la función seno para que oscile entre minVal y maxVal
  # x es la iteración actual, que afecta la fase del seno
  return(amplitud * sin(x) + centro)
}

# Valores máximos y mínimos deseados para la función seno
minVal <- 0.01
maxVal <- 0.9

# Simulación de las iteraciones con diferentes valores de 'x'
x <- seq(0, 2 * pi, length.out = 10000)

y <- sapply(x, senoPersonalizado, minVal = minVal, maxVal = maxVal)

# Graficar el resultado
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(data.frame(x, y), ggplot2::aes(x, y)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Función Seno Personalizada",
                  x = "x",
                  y = "senoPersonalizado") +
    ggplot2::theme_minimal()
} else {
  plot(x, y, type = 'l', main = "Función Seno Personalizada", xlab = "x", ylab = "senoPersonalizado")
}




library(ggplot2)

set.seed(13243) # Para reproducibilidad

# Punto central aproximado (e.g., Santiago de Chile)
lat_central <- -33.45
lon_central <- -70.65

# Conversión de 5 km a grados
km_en_grados_lat <- 5 / 111 # Aproximadamente 0.045 grados
km_en_grados_lon <- 5 / (111 * cos(lat_central * pi / 180)) # Ajuste por coseno de la latitud

n_colegios <- 10
n_estudiantes <- sum(sample(50:100, n_colegios, replace = TRUE)) # Total de estudiantes basado en capacidad de colegios

# Crear dataframe de colegios dentro de 5km del punto central
colegios <- data.frame(
  ID = 1:n_colegios,
  X = lon_central + runif(n_colegios, min = -km_en_grados_lon, max = km_en_grados_lon),
  Y = lat_central + runif(n_colegios, min = -km_en_grados_lat, max = km_en_grados_lat),
  num_alumnos = rep(0, n_colegios), # Inicializar en 0, se actualizará después
  num_sep = round(runif(n_colegios, min = 10, max = 30))
)

# Generar estudiantes dentro del área de 5km
alumnos <- data.frame(
  ID_COL = integer(n_estudiantes), # Se llenará más adelante
  X = lon_central + runif(n_estudiantes, min = -km_en_grados_lon, max = km_en_grados_lon),
  Y = lat_central + runif(n_estudiantes, min = -km_en_grados_lat, max = km_en_grados_lat),
  sep = sample(0:1, n_estudiantes, replace = TRUE, prob = c(0.7, 0.3)) # 30% de probabilidad de ser SEP
)

# Asignar estudiantes a colegios de manera equitativa y ajustar num_alumnos en colegios
alumnos$ID_COL <- rep(1:n_colegios, length.out = n_estudiantes)
alumnos <- alumnos[sample(nrow(alumnos)), ] # Mezclar alumnos para distribución aleatoria

# Actualizar num_alumnos en colegios basado en la asignación
for (i in 1:n_colegios) {
  colegios$num_alumnos[i] <- sum(alumnos$ID_COL == i)
}

# Gráfico con formas distintas para cada colegio
ggplot(data = colegios, aes(x = X, y = Y)) +
  geom_point(color = "red", size = 4, alpha = 0.6) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  labs(title = "Distribución de Colegios y Alumnos en un Área de 5 km",
       x = "Longitud",
       y = "Latitud",
       color = "Colegio",
       shape = "Colegio") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  # Añadir estudiantes de cada colegio con diferentes formas
  geom_point(data = alumnos, aes(x = X, y = Y, shape = factor(ID_COL), color = factor(sep)), alpha = 0.5) +
  scale_color_manual(values = c('0' = 'blue', '1' = 'green')) +
  guides(shape = guide_legend(title = "ID Colegio"), color = guide_legend(title = "Estudiante SEP"))

print(head(colegios))
print(head(alumnos))

# Exportar el dataframe de colegios a CSV sin nombres de columnas
write.table(colegios, "colegios_test.txt", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)

# Exportar el dataframe de alumnos a CSV sin nombres de columnas
write.table(alumnos, "alumnos_test.txt", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)


####################################


# Cargar la librería dplyr para las funciones de manipulación de datos
library(dplyr)

# Define los vectores de posibles valores para n_block y n_thread
n_block_values <- c(1,2,4,6,8,32,64,128,256,512,1024)
n_thread_values <- c(32,64,128,256,512,1024)

# Crea todas las posibles combinaciones de n_block y n_thread
combinations <- expand.grid(n_block = n_block_values, n_thread = n_thread_values)

# Calcula el producto de n_block y n_thread y crea un nuevo dataframe con los resultados deseados
desired_combinations <- combinations %>%
  mutate(product = n_block * n_thread) %>%
  filter(product %in% c(32, 64, 256, 512, 1024,2048, 4096)) %>%
  arrange(product)

# Muestra el dataframe resultante
print(desired_combinations)




############################### Analizar Resultados

library(readr)
library(dplyr)
library(ggplot2)
library(car)

datos <- read_csv("seed_iteration.dat")

datos <-read_csv("seed_iteration_2.dat")

datos <-read_csv("seed_iteration_3.dat")

datos <- read_csv("ultimo.csv")

datos <- read_csv("finalfinal.csv")

datos <- rbind(datos,datos2)

datos <- cbind(datos, block_x_thread=c((datos$n_block*datos$n_thread)))

datos <- select(datos, z, it, n_block,n_thread,block_x_thread, temp, coolingRate, len1, len2)

modelo <- lm(z ~ block_x_thread, data = datos)
modelo <- lm(z ~ n_block + n_thread, data = datos)
library(lmtest)

modelo <- lm((z) ~ (n_block) + (n_thread), data = datos)
summary(modelo)
bptest(modelo) # Modelo homoscedástico de hilos y bloques en cada réplica
plot(modelo)
shapiro.test(modelo$residuals)


modelo <- lm(log(z) ~ log(n_block) + log(n_thread), data = datos)
summary(modelo)
bptest(modelo) # Modelo heteroscedástico de hilos y bloques en cada réplica
plot(modelo)
shapiro.test(modelo$residuals)

# Ajustar un modelo lineal con las variables que son significativas
modelo <- lm((z) ~ (n_block) + (n_thread), data = datos)

block_grid <- seq(min(datos$n_block), max(datos$n_block), length.out = 100)
thread_grid <- seq(min(datos$n_thread), max(datos$n_thread), length.out = 100)
grid <- expand.grid(n_block = block_grid, n_thread = thread_grid)

# Predecimos 'z' usando el modelo y nuestro nuevo grid
grid$z <- predict(modelo, newdata = grid)

ggplot(data=grid)+
  geom_contour_filled(aes(x=n_block,y=n_thread,z=z))+
  scale_x_continuous(trans= "log")+
  scale_y_continuous(trans= "log")
  

modelo <- lm(log(z) ~ log(block_x_thread), data = datos)
summary(modelo) # Valor del exponente
bptest(modelo) # Modelo homoscedástico de hilos totales en cada réplica
plot(modelo)
shapiro.test(modelo$residuals)

library(mgcv)
modelo<-gam((z)~s((n_block),k=6)+s((n_thread),k=6),data=datos)
summary(modelo)
plot(modelo)

library(mgcv)
modelo<-gam(log(z)~s(log(block_x_thread),k=4),data=datos)
summary(modelo)
plot(modelo)

library(mgcv)
modelo<-gam((z)~s((block_x_thread),k=6),data=datos)
summary(modelo)
plot(modelo) #BINGO!!! Esto es lo que necesitamos

library(mgcv)
modelo<-gam((z)~te(n_block,n_thread,k=6),data=datos)
summary(modelo)
plot(modelo)

modelo <- lm(z ~ log(n_block) + log(n_thread) + log(len2), data = datos)
summary(modelo)
vif(modelo)


modelo<-lm(temp~block_x_thread, data=datos)
summary(modelo)

par(mfrow=c(2,2))
plot(modelo)

library(mgcv)
modelo <- gam(z~ s(log(block_x_thread),k=4) + s((it)) + s(coolingRate) + s(len1) + s(len2), data = datos)
summary(modelo)
plot(modelo)

datos$block_x_thread_f<-factor(datos$block_x_thread)

ggplot(datos)+
  geom_boxplot(aes(x=block_x_thread_f,y=it))

ggplot(datos)+
  geom_boxplot(aes(x=block_x_thread_f,y=z))


# Gráfico de n_block vs. z
ggplot(data = datos, aes(x = n_block, y = z)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Gráfico de n_block vs. z",
       x = "n_block",
       y = "z")

# Gráfico de n_thread vs. z
ggplot(data = datos, aes(x = n_thread, y = z)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Gráfico de n_thread vs. z",
       x = "n_thread",
       y = "z")


# Gráfico de n_thread vs. z
ggplot(data = datos, aes(x = block_x_thread, y = z)) +
  geom_point() +
  scale_x_continuous(trans= "log")+
  scale_y_continuous(trans= "log")+
  geom_smooth(method = "lm", col = "red",se=T) +
  labs(title = "Gráfico de total threads vs. z",
       x = "n_thread",
       y = "z") # Gráfico para publicación

# Gráfico de n_thread vs. z
ggplot(data = datos, aes(x = n_block, y = z)) +
  geom_point() +
  scale_x_continuous(trans= "log")+
  scale_y_continuous(trans= "log")+
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Gráfico de n_block vs. z",
       x = "n_block",
       y = "z")

# Gráfico de n_thread vs. z
ggplot(data = datos, aes(x = n_thread, y = z)) +
  geom_point() +
  scale_x_continuous(trans= "log")+
  scale_y_continuous(trans= "log")+
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Gráfico de n_thread vs. z",
       x = "n_thread",
       y = "z")




# Cargar la librería GGally para gráficos de pares mejorados
#install.packages("GGally")
library(GGally)

# Selecciona solo las variables significativas y z
significant_vars <- datos[c("z", "it", "coolingRate", "len1", "len2","block_x_thread")]

# Crea gráficos de pares
ggpairs(significant_vars)



library(ggplot2)
### grafico importanteee!!!!!
# Crear un gráfico de dispersión para 'it' y 'z'
ggplot(data = datos, aes(x = it, y = z,col=factor(block_x_thread))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'it'",
       y = "'z'")
  #scale_x_continuous(trans= "log")+
  #scale_y_continuous(trans= "log")

datos2<-datos[which(datos$block_x_thread==4096),]
datos2<-datos2[-which.min(datos2$z),]
ggplot(data = datos2, 
       aes(x = it, y = z,col=factor(n_block))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'it'",
       y = "'z'")

ggplot(data = datos2, 
       aes(x = it, y = z,col=factor(n_thread))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'it'",
       y = "'z'")
#scale_x_continuous(trans= "log")+
#scale_y_continuous(trans= "log")

modelo<-lm(z~it*factor(n_block),datos2)
summary(modelo)

modelo<-lm(z~it*factor(n_thread),datos2)
summary(modelo)

datos3<-datos[which(datos$block_x_thread==4096),]
modelo<-lm(z~it*factor(n_block),datos3)
summary(modelo)

modelo<-lm(z~it*factor(n_thread),datos3)
summary(modelo)

#Regresiones robustas
library(MASS)
library(sfsmisc)

modelo<-rlm(z~it*factor(n_block),datos3)
summary(modelo)
f.robftest(modelo, var = "factor(n_block)128")

ggplot(data = datos3, aes(x = it, y = z,col=factor(n_thread))) +
  geom_point() +
stat_smooth(method="rlm",
            fullrange=TRUE)

ggplot(data = datos, aes(x = it, y = z,col=factor(block_x_thread))) +
  geom_point() +
  stat_smooth(method="rlm",
              fullrange=TRUE)


modelo<-rlm(z~it*factor(n_thread),datos3)
summary(modelo)
f.robftest(modelo, var = "factor(n_block)128")

ggplot(data = datos, aes(x = it, y = z,col=factor(block_x_thread))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'it'",
       y = "'z'")
#scale_x_continuous(trans= "log")+
#scale_y_continuous(trans= "log")

ggplot(data = datos, aes(x = coolingRate, y = z)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'coolingRate'",
       y = "'z'")
ggplot(data = datos, aes(x = len1, y = z)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'len1'",
       y = "'z'")
ggplot(data = datos, aes(x = len2, y = z)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'len2'",
       y = "'z'")
ggplot(data = datos, aes(x = block_x_thread, y = z)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relación entre 'it' y 'z'",
       x = "'len2'",
       y = "'z'")

ggplot(data = datos, aes(x = block_x_thread, y = it)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relación entre 'it' y 'blocxthread'",
       x = "'block_x_thread'",
       y = "'it'")
ggplot(data = datos, aes(x = n_block, y = it)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relación entre 'it' y 'block'",
       x = "'block'",
       y = "'it'")
ggplot(data = datos, aes(x = n_thread, y = it)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relación entre 'it' y 'thread'",
       x = "'thread'",
       y = "'it'")



#install.packages("plotly")
library(plotly)

# Ajustar un modelo lineal con las variables que son significativas
modelo <- lm(log(z) ~ log(n_block) + log(n_thread), data = datos)

block_grid <- seq(min(datos$n_block), max(datos$n_block), length.out = 100)
thread_grid <- seq(min(datos$n_thread), max(datos$n_thread), length.out = 100)
grid <- expand.grid(n_block = block_grid, n_thread = thread_grid)

# Predecimos 'z' usando el modelo y nuestro nuevo grid
grid$z <- predict(modelo, newdata = grid)


# Creamos el gráfico de contorno 3D
fig <- plot_ly(x = ~grid$n_block, y = ~grid$n_thread, z = ~grid$z, type = "contour")

# Añadir características adicionales si es necesario
fig <- fig %>% layout(title = "3D Contour Plot of z with n_block and n_thread",
                      scene = list(
                        xaxis = list(title = "block"),
                        yaxis = list(title = "thread"),
                        zaxis = list(title = "z")))
#fig <- layout(fig,xaxis=list(type="log"),yaxis=list(type="log"))

# Mostrar el gráfico
fig # BINGO!!! El otro gráfico. Hay una caida exponencial de z en función de bloques e hilos

library(ggplot2)

info <- read.csv("./info-graphics.txt", sep=",", dec=".", header = FALSE)
position_solution <- read.csv("./info-graphicsBestSolution.txt", sep=",", dec=".", header = FALSE)
colnames(info) <- c("it",'dn','d','s','c','z','temp')
info$difz<- c(0, abs(diff(x = info$z)))


ggplot(data=info[1:9500,],aes(x=it, y=z))+
  geom_line()
  #scale_x_continuous(trans = 'log')+
  #scale_y_continuous(trans = 'log')

ggplot(data=info[1:1000,],aes(x=it, y=z))+
  geom_line()

ggplot(data=info,aes(x=it, y=difz))+
  geom_line()+
  scale_x_continuous(trans = 'log')+
  scale_y_continuous(trans = 'log')

ggplot(data=info,aes(x=it, y=difz))+
  geom_smooth()+
  scale_x_continuous(trans = 'log')+
  scale_y_continuous(trans = 'log')

###################################################################
## V1

cromosoma <- read.csv("./alumnos_utm.txt", sep=",", dec=".", header = FALSE)
colnames(cromosoma)<- c('id','x','y','sep')
cromosoma1 <- position_solution[1,]

cromosoma2 <- position_solution[2,]
escuelas <- read.csv("./colegios_utm.txt", sep=",", dec=".", header = FALSE)
colnames(escuelas) <- c('id','x','y','total','sep')
escuelas<-escuelas[order(escuelas$id),]
cromosoma<-cromosoma[order(cromosoma$id),]




library(tidyverse)
hull <- cromosoma %>% group_by(id) %>% 
      slice(chull(x, y))
hull$id<-factor(hull$id)
cromosoma$id<-factor(cromosoma$id)
escuelas$id<-factor(escuelas$id)


for(i in hull$id){
pdf(paste(i,".pdf",sep=""))
print(ggplot() + geom_polygon(data = hull[which(hull$id==i),], alpha = 0.2, 
                    aes(x,y,fill = id,colour = id))+
        geom_point(data=cromosoma[which(cromosoma$id==i),],aes(x,y,col=id))+
        geom_point(data=escuelas[which(escuelas$id==i),],aes(x,y,col=id),pch=3))
dev.off()
}


#####################################
### V2

###################################################################


position_solution <- read.csv("./info-graphicsBestSolution.txt", sep=",", dec=".", header = FALSE)
cromosoma1 <- position_solution[1,]
cromosoma2 <- as.numeric(position_solution[2,])


cromosoma <- read.csv("./alumnos_utm.txt", sep=",", dec=".", header = FALSE)
colnames(cromosoma)<- c('id','x','y','sep')
cromosoma$id <- cromosoma2[-length(cromosoma2)]+1

escuelas <- read.csv("./colegios_utm.txt", sep=",", dec=".", header = FALSE)
colnames(escuelas) <- c('id','x','y','total','sep')
escuelas$id <- seq_along(escuelas[,1])


escuelas<-escuelas[order(escuelas$id),]
cromosoma<-cromosoma[order(cromosoma$id),]




library(tidyverse)
hull <- cromosoma %>% group_by(id) %>% 
  slice(chull(x, y))


hull$id<-factor(hull$id)
cromosoma$id<-factor(cromosoma$id)
escuelas$id<-factor(escuelas$id)
cromosoma$sep <-factor(cromosoma$sep)


for(i in hull$id){
  pdf(paste(i,".pdf",sep=""))
  print(ggplot() + geom_polygon(data = hull[which(hull$id==i),], alpha = 0.2, 
                                aes(x,y,fill = id,colour = id))+
          geom_point(data=cromosoma[which(cromosoma$id==i),],aes(x,y,col=sep+1))+
          geom_point(data=escuelas[which(escuelas$id==i),],aes(x,y,col=id), shape = 22, size = 3, stroke = 1,color="red"))
  dev.off()
}


i<-28
ggplot()+
  geom_point(data=cromosoma,aes(x,y,col=id))+
  geom_point(data=escuelas,aes(x,y,col=id), shape = 22, size = 3, stroke = 1)


ggplot() + geom_polygon(data = hull, alpha = 1, 
                        aes(x,y,color = id))+
  geom_point(data=cromosoma,aes(x,y,col=sep,fill=id),shape=21,alpha = 0.2)+
  geom_point(data=escuelas,aes(x,y,col=id), shape = 22, size = 3, stroke = 1)





plot.esc<-function(cromosoma,escuelas,color=T){
  #  browser()
  require(RColorBrewer)
  cromosoma_orig<-cromosoma
  cualesEsc<-as.numeric(names(table(cromosoma))) #identificar escuelas abiertas
  #  cromosoma<-cromosomas[[1]]
  cuales<-cromosoma %in% escuelas #Filtrar estudiantes según escuelas
  cromosoma<-factor(cromosoma[cuales],cualesEsc)
  vulni<-factor(vuln[cuales],0:1) #Ver estudiantes vulnerables
  #  posEstEsc<-as.data.frame(posEst[cromosomas[[1]] %in% escuelas,])
  posEstEsc<-as.data.frame(posEst[cromosoma_orig %in% escuelas,])
  posEstEsc<-cbind(posEstEsc,Esc=cromosoma,Vuln=vulni)
  posEsci<-cbind(as.data.frame(posEsc[cualesEsc,]),id=factor(cualesEsc,cualesEsc))
  posEsci<-posEsci[posEsci$id %in% escuelas,]
  posEsci<-cbind(posEsci,Si=unlist(lapply(cualesEsc,function(i) Sesc(cromosoma_orig,i))))
  
  rayos<-data.frame(X=0,Y=0,Esc=factor(unlist(lapply(1:nest,function(i) return(c(cromosoma_orig[i],cromosoma_orig[i])))),cualesEsc))
  
  for(i in unique(posEsci$id)){
    cualesEstEsc<-which(posEstEsc$Esc==i)
    rayos[2*cualesEstEsc-1,1:2]<-posEsci[i,1:2]
    rayos[2*cualesEstEsc,1:2]<-posEstEsc[cualesEstEsc,1:2]
    if(i==unique(posEsci$id)[1]){
      findhull<-posEstEsc[cualesEstEsc[chull(posEstEsc$X[cualesEstEsc],
                                             posEstEsc$Y[cualesEstEsc])],]
    } else {
      findhull<-rbind(findhull,posEstEsc[cualesEstEsc[chull(posEstEsc$X[cualesEstEsc],
                                                            posEstEsc$Y[cualesEstEsc])],])
    }
  }
  #browser()
  posEstEsc$Esc<-factor(posEstEsc$Esc,escuelas)
  posEsci$id<-factor(posEsci$id,escuelas)
  myColors <- colorRampPalette(brewer.pal(8, "Dark2"))(nesc)
  names(myColors) <- levels(posEsci$id)
  
  p<-ggplot()
  p<-p+geom_polygon(data=findhull,aes(x=X,y=Y,fill=Esc),alpha=0.1)
  p<-p+geom_polygon(data=findhull,aes(x=X,y=Y,col=Esc),alpha=0)
  #p<-p+geom_path(data=rayos,aes(x=X,y=Y,color=Esc))
  p<-p+geom_point(data=posEsci,aes(x=X,y=Y,size=Si,fill=id),shape=25,color="black")
  p<-p+scale_size(limits = c(0,1),breaks=c(0,0.3,0.6,1),labels=c("Nula","Alta","Aberrante","Absoluta"))
  #p<-p+scale_color_manual(values=c("lightred","pink","lightblue","lightgreen","orange"))
  #p<-p+scale_fill_manual(values=c("lightred","pink","lightblue","lightgreen","orange"))
  p<-p+scale_color_manual(values=myColors)
  p<-p+scale_fill_manual(values=myColors)
  p<-p+geom_point(data=posEstEsc,aes(x=X,y=Y,color=Esc,shape=Vuln),size=3)
  #p<-p+theme(legend.position="bottom",legend.box="horizontal")
  p<-p+annotate("text",label=paste("S=",round(S(cromosoma),3),sep=""),
                x=min(posEstEsc$X),y=min(posEstEsc$Y),size=3,hjust = 0)
  if(!color) p<-p+scale_color_grey()+scale_fill_grey()
  p
}


