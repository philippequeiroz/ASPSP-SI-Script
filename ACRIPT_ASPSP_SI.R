##### SCRIPT_ASPSP_SI #####

###### CDA (Canonical Discriminant Analysis) ######

install.packages("candisc")   # apenas se ainda não estiver instalado
install.packages("ggplot2")   # para visualização bonita
library(candisc)
library(ggplot2)

# --- 1. Importar seus dados ---
# Carregar os dados
dados <- read.table("clipboard", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

dados$Species <- as.factor(dados$Species)

manova_model <- manova(cbind(d13C, d15N) ~ Species, data = dados)
summary(manova_model)

cda_model <- candisc(manova_model, term = "Species")
summary(cda_model)

# Dados para o gráfico
cda_scores <- as.data.frame(cda_model$scores)
cda_scores$Species <- dados$Species

# Plot com ggplot2
ggplot(cda_scores, aes(x = Can1, y = Can2, color = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.68, size = 1) +  # elipses de confiança (68%)
  labs(title = "CDA: δ¹³C e δ¹⁵N por espécie",
       x = "Canônica 1", y = "Canônica 2") +
  theme_minimal()
30
##### nÍVEL tRÓFICO E PLOTS #####


# Pacotes necessários
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(ellipse)

# --- 1. Importar seus dados ---
# Carregar os dados
dados <- read.table("clipboard", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# --- 2. Separar predadores e presas ---
predadores <- dados %>% filter(Species %in% c("E.bipinnulata", "T.albacares", "A.solandri"))
presas <- dados %>% filter(Species %in% c("FlyingFish", "Shrimp", "Squid"))

# --- 3. Resumo das presas com TL específico ---
baseline <- presas %>%
  group_by(Species) %>%
  summarise(mean_d15N = mean(d15N, na.rm = TRUE),
            mean_d13C = mean(d13C, na.rm = TRUE),
            TL = mean(TP_baseline, na.rm = TRUE))

# --- 4. Definir fracionamento isotópico ---
delta15N <- 1.9
delta13C <- 1.8

# --- 5. Função para calcular proporções de cada presa ---
calc_mix <- function(pred_d15N, baseline_df) {
  diff <- abs(pred_d15N - baseline_df$mean_d15N)
  prop <- 1 / diff
  prop <- prop / sum(prop)
  return(prop)
}

# --- 6. Calcular TP ponderado para cada predador ---
resultado <- predadores %>%
  rowwise() %>%
  mutate(
    prop_FlyingFish = calc_mix(d15N, baseline)[baseline$Species == "FlyingFish"],
    prop_Shrimp      = calc_mix(d15N, baseline)[baseline$Species == "Shrimp"],
    prop_Squid       = calc_mix(d15N, baseline)[baseline$Species == "Squid"],
    TP_mix = prop_FlyingFish * (baseline$TL[baseline$Species=="FlyingFish"] + (d15N - baseline$mean_d15N[baseline$Species=="FlyingFish"]) / delta15N) +
      prop_Shrimp      * (baseline$TL[baseline$Species=="Shrimp"]      + (d15N - baseline$mean_d15N[baseline$Species=="Shrimp"]) / delta15N) +
      prop_Squid       * (baseline$TL[baseline$Species=="Squid"]       + (d15N - baseline$mean_d15N[baseline$Species=="Squid"]) / delta15N)
  ) %>%
  ungroup()

# --- 7. Visualizar resultado ---
head(resultado[, c("Species", "d15N", "TP_mix")])

# Ver quantos indivíduos há no resultado
nrow(resultado)

# Resumo do TP por espécie
resumo_TP <- resultado %>%
  group_by(Species) %>%
  summarise(
    n_individuos = n(),
    mean_TP = mean(TP_mix, na.rm = TRUE),
    sd_TP = sd(TP_mix, na.rm = TRUE)
  )

# Visualizar o resumo
resumo_TP


# --- 8. Gráfico de TP por espécie ---
ggplot(resultado, aes(x = Species, y = TP_mix, fill = Species)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Nível Trófico estimado (TP) dos predadores",
       y = "TP estimado", x = "Espécie") +
  scale_fill_brewer(palette = "Set2")


###### OUTRO GRÁFICO ######

# Criar resumo por espécie
resumo <- resultado %>%
  group_by(Species) %>%
  summarise(
    mean_TP = mean(TP_mix, na.rm = TRUE),
    sd_TP = sd(TP_mix, na.rm = TRUE)
  )

# Gráfico estilo Figure 2a
ggplot() +
  # Pontos individuais dos predadores
  geom_jitter(data = resultado, aes(x = Species, y = TP_mix, color = Species),
              width = 0.15, alpha = 0.7, size = 2) +
  # Média de cada espécie
  geom_point(data = resumo, aes(x = Species, y = mean_TP),
             color = "black", size = 3, shape = 18) +
  # Desvio padrão
  geom_errorbar(data = resumo, aes(x = Species, ymin = mean_TP - sd_TP, ymax = mean_TP + sd_TP),
                width = 0.2, color = "black") +
  theme_minimal() +
  labs(title = "Nível trófico estimado dos predadores",
       x = "Espécie", y = "TP estimado") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "none")



######## Gráfico de densidade por espécie

ggplot(resultado, aes(x = TP_mix, fill = Species)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribuição do nível trófico dos predadores",
       x = "TP estimado", y = "Densidade") +
  scale_fill_brewer(palette = "Set2")



######### Gráfico de violino (violin plot)

ggplot(resultado, aes(x = Species, y = TP_mix, fill = Species)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2, color = "black") +
  geom_point(data = resumo, aes(x = Species, y = mean_TP), size = 3, shape = 18, color = "red") +
  theme_minimal() +
  labs(title = "Nível trófico estimado (TP) por espécie",
       x = "Espécie", y = "TP estimado") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")


####### Gráfico de barras com erro padrão ou SD

ggplot(resumo, aes(x = Species, y = mean_TP, fill = Species)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_TP - sd_TP, ymax = mean_TP + sd_TP), width = 0.2, color = "black") +
  theme_minimal() +
  labs(title = "Média do nível trófico por espécie",
       x = "Espécie", y = "TP estimado") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")


####### Gráfico de dispersão TP vs δ¹⁵N ou δ¹³C

ggplot(resultado, aes(x = d15N, y = TP_mix, color = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "TP estimado vs δ15N", x = expression(delta^15*N~"\u2030"), y = "TP estimado") +
  scale_color_brewer(palette = "Set2")



####### Scatter plot 3D ou plot “biplot” colorido

ggplot(resultado, aes(x = d13C, y = d15N, color = TP_mix, size = TP_mix)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "δ13C vs δ15N com TP estimado",
       x = expression(delta^13*C~"\u2030"),
       y = expression(delta^15*N~"\u2030"),
       color = "TP estimado",
       size = "TP estimado")



####### Plot 2D com elipses por espécie


# Exemplo simples
ggplot(resultado, aes(x = d13C, y = d15N, color = Species)) +
  geom_point(size = 2, alpha = 0.6) +
  stat_ellipse(level = 0.95, aes(fill = Species), alpha = 0.2, geom = "polygon") +
  theme_minimal() +
  labs(title = "Isotopic niche com TP sobreposto",
       x = expression(delta^13*C~"\u2030"),
       y = expression(delta^15*N~"\u2030"))



# Gráfico δ13C vs δ15N com TP como cor
ggplot(resultado, aes(x = d13C, y = d15N, color = TP_mix)) +
  geom_point(aes(shape = Species), size = 3, alpha = 0.8) +  # Diferencia espécies por shape
  scale_color_viridis_c(option = "D", name = "TP estimado") + # Paleta contínua
  theme_minimal(base_size = 14) +
  labs(
    title = "δ13C vs δ15N com nível trófico (TP) estimado",
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030")
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )



###### Heatmap ou escala de cores do TP

# Supondo que 'resultado' já contenha: Species, d13C, d15N, TP_mix

# Gráfico δ13C vs δ15N com tamanho pelo TP e cor pela espécie
ggplot(resultado, aes(x = d13C, y = d15N, color = Species, size = TP_mix)) +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +   # Paleta para espécies
  theme_minimal(base_size = 14) +
  labs(
    title = "δ13C vs δ15N com nível trófico (TP) estimado",
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    size = "TP estimado"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


##### Scatter plot 3D ou plot “biplot” colorido + densidades

# Gráfico base: δ13C vs δ15N, cor = Species, tamanho = TP
grafico_base <- ggplot(resultado, aes(x = d13C, y = d15N, color = Species, size = TP_mix)) +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "δ13C vs δ15N com nível trófico (TP) estimado",
    x = expression(delta^13*C~"\u2030"),
    y = expression(delta^15*N~"\u2030"),
    size = "TP estimado"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Adiciona distribuições marginais
ggMarginal(
  grafico_base, 
  type = "density",      # Tipo de distribuição: density
  groupColour = TRUE,    # Mantém cores por espécie
  groupFill = TRUE       # Preenche área das densidades por espécie
)


# SALVAR

# Supondo que 'grafico_base' já foi criado
# Exemplo de salvamento em TIF de alta qualidade
tiff("delta13C_delta15N_TP.tif", 
     width = 5000, height = 3000,  # dimensões em pixels
     res = 600,                    # resolução em dpi
     compression = "lzw")          # compressão sem perda

# Plota o gráfico com marginais
ggMarginal(
  grafico_base,
  type = "density",
  groupColour = TRUE,
  groupFill = TRUE
)

dev.off()  # Finaliza a escrita do arquivo





ggplot(resultado, aes(x = d13C, y = TP_mix, color = Species, size = d15N)) +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "δ13C vs Nível Trófico (TP) estimado",
    x = expression(delta^13*C~"\u2030"),
    y = "TP estimado",
    color = "Espécie",
    size = expression(delta^15*N~"\u2030")
  )





ggplot(resultado, aes(x = d13C, y = TP_mix, color = Species, size = d15N)) +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "δ13C vs Nível Trófico (TP) estimado",
    x = expression(delta^13*C~"\u2030"),
    y = "TP estimado",
    color = "Espécie",
    size = expression(delta^15*N~"\u2030")
  )

library(ggplot2)
library(dplyr)

ggplot(resultado, aes(x = d13C, y = TP_mix, color = Species, size = TP_mix)) +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +   # cores por espécie
  theme_minimal(base_size = 14) +
  labs(
    title = "δ13C vs Nível Trófico (TP) estimado",
    x = expression(delta^13*C~"\u2030"),
    y = "TP estimado",
    color = "Espécie",
    size = "TP estimado"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

31
######## BIPLOT COM DENSIDADES #########


# Instalar pacotes se necessário
# install.packages("ggplot2")
# install.packages("dplyr")

# Carregar pacotes
library(ggplot2)
library(dplyr)

# Importar os dados
dados <- read.table("clipboard", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


##### SÓ O BIPLOT ######


# Checar se as colunas têm os nomes corretos
head(dados)

# Criar o gráfico semelhante à Figura 2a
ggplot(dados, aes(x = d13C, y = d15N, color = Species, fill = Species)) +
  geom_point(size = 3, alpha = 0.8, shape = 21, stroke = 1) +
  stat_ellipse(geom = "polygon", alpha = 0.2) +  # Ellipses para cada grupo
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    text = element_text(size = 14)
  ) +
  labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030")),
    title = "Isótopos estáveis por espécie"
  )



###### BIPLOT + DENSIDADE #######


install.packages("ggplot2")
install.packages("ggExtra")
install.packages("dplyr")

library(ggplot2)
library(ggExtra)
library(dplyr)

# Filtrar apenas grupos com pelo menos 3 pontos (necessário para as elipses)
dados_filtrados <- dados %>%
  group_by(Species) %>%
  filter(n() >= 3) %>%
  ungroup()

# Gráfico base
grafico_base <- ggplot(dados_filtrados, aes(x = d13C, y = d15N, color = Species, fill = Species)) +
  geom_point(size = 3, alpha = 0.8, shape = 21, stroke = 1) +
  stat_ellipse(geom = "polygon", alpha = 0.2, color = NA) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    text = element_text(size = 14)
  ) +
  labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030")),
    title = "Nicho isotópico das espécies (δ¹³C vs δ¹⁵N)"
  )

# Adicionar densidades marginais
ggMarginal(
  grafico_base,
  type = "density",
  groupColour = TRUE,
  groupFill = TRUE
)



######## BIPLOT + DENSIDADE SEM ELIPSES ########


# Pacotes necessários
library(ggplot2)
library(ggExtra)
library(dplyr)

# Calcular médias e desvios padrão por espécie
media_sd <- dados %>%
  group_by(Species) %>%
  summarise(
    mean_d13C = mean(d13C, na.rm = TRUE),
    mean_d15N = mean(d15N, na.rm = TRUE),
    sd_d13C = sd(d13C, na.rm = TRUE),
    sd_d15N = sd(d15N, na.rm = TRUE)
  )

# Gráfico base com pontos
grafico_base <- ggplot(dados, aes(x = d13C, y = d15N, color = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  
  # Adiciona cruzes (média ± 1 SD)
  geom_errorbar(data = media_sd, 
                aes(x = mean_d13C, ymin = mean_d15N - sd_d15N, ymax = mean_d15N + sd_d15N, color = Species),
                width = 0.1, inherit.aes = FALSE) +
  geom_errorbarh(data = media_sd, 
                 aes(y = mean_d15N, xmin = mean_d13C - sd_d13C, xmax = mean_d13C + sd_d13C, color = Species),
                 height = 0.1, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = expression(delta^13*C~"\u2030"),
       y = expression(delta^15*N~"\u2030"),
       title = "Isotopic niche with mean ± SD per group")

# Adiciona distribuições marginais
ggMarginal(grafico_base, type = "density", groupColour = TRUE, groupFill = TRUE)

32
######## NICHO ISOTÓPICO #######

##### SIBER #####

install.packages("SIBER")
library(SIBER)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

# Importar os dados
siber_data <- read.table("clipboard", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Preparar os dados
siber_data <- dados %>%
  select(d13C, d15N, Species) %>%
  rename(iso1 = d13C, iso2 = d15N, group = Species) %>%
  mutate(community = 1)  # todas as espécies na mesma comunidade

# 3. Criar objeto SIBER
siber_object <- createSiberObject(siber_data)

# 4. Plot simples (opcional)
plotSiberObject(siber_object)

# 5. Definir priors
priors <- list(
  R = diag(2),
  k = 2,
  mu = c(0, 0),
  V = diag(2) * 100,
  tau.mu = 0.001  # necessário para evitar erro do JAGS
)

# 6. Parâmetros MCMC
parms <- list(n.iter = 10000, n.burnin = 1000, n.thin = 10, n.chains = 2)

# 7. Rodar modelo bayesiano
ellipses <- siberMVN(siber_object, parms, priors)

# 8. Calcular SEA.B (área elíptica bayesiana)
SEA.B <- siberEllipses(ellipses)  # matriz com 2000 linhas (amostras) x N grupos (espécies)

# 9. Converter para formato tidy (long)
SEA.df <- as.data.frame(SEA.B)
colnames(SEA.df) <- levels(as.factor(siber_data$group))  # nomes das espécies
SEA_long <- SEA.df %>%
  pivot_longer(cols = everything(), names_to = "Species", values_to = "SEA_B")

# 10. Plot da SEA.B com boxplot
ggplot(SEA_long, aes(x = Species, y = SEA_B, fill = Species)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Área elíptica bayesiana (SEA.B)",
       x = "Espécie", y = expression("SEA"[B]~"(‰²)")) +
  scale_fill_brewer(palette = "Set2")





################################ nicheROVER ################################
rm(list = ls())

library(nicheROVER)
library(dplyr)

# ----------------------------- Importação ---------------------------------
ray <- read.table("clipboard", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ------------------------- Preparar dados ---------------------------------
ray <- ray %>%
  mutate(
    d13C = as.numeric(d13C),
    d15N = as.numeric(d15N)
  ) %>%
  select(Species, d13C, d15N)

# ------------------------- Configurações ----------------------------------
species.names <- c("Ebip", "Talb", "Asol")  # ordem fixa

# Cores pastéis correspondentes
clrs <- c(
  "Ebip" = "#F4A86A",    # salmão
  "Talb" = "#7b9acc",    # azul escuro pastel
  "Asol" = "#7FC8A9"     # verde claro
)

# ------------------------- Parâmetros Bayesianos --------------------------
nsamples <- 1000
ray.par <- lapply(species.names, function(sp) {
  ii <- which(ray$Species == sp)
  niw.post(nsamples = nsamples, X = ray[ii, c("d13C","d15N")])
})
names(ray.par) <- species.names

# ----------------------------- Gráfico de parâmetros ----------------------
tiff("grafico_niche_parameters.tiff", width = 2000, height = 800, res = 300)
par(mar = c(4, 4, 2, 1), mfrow = c(1, 3))
niche.par.plot(ray.par, col = clrs, plot.index = 1)
niche.par.plot(ray.par, col = clrs, plot.index = 2)
niche.par.plot(ray.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = species.names, fill = clrs, cex = 0.8, bg = NA, bty = "n")
dev.off()

# ----------------------------- Projeção 2D --------------------------------
ray.data <- lapply(species.names, function(sp) {
  ray[ray$Species == sp, c("d13C","d15N")]
})
names(ray.data) <- species.names

tiff("grafico_projecao_2D.tiff", width = 3000, height = 1800, res = 300)
niche.plot(niche.par = ray.par, niche.data = ray.data, pfrac = 0.05,
           iso.names = expression(delta^{13}*C, delta^{15}*N),
           col = clrs, xlab = expression("Isotope Ratio (‰)"))
legend("topright", legend = species.names, fill = clrs, cex = 0.6, bty = "n", bg = NA)
dev.off()

# ----------------------------- Sobreposição de nicho ----------------------
over.stat <- overlap(ray.par, nreps = nsamples, nprob = 1000, alpha = 0.95)
names(over.stat) <- species.names

tiff("grafico_overlap.tiff", width = 3000, height = 1800, res = 300)
overlap.plot(over.stat, col = clrs, mean.cred.col = "red", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")
dev.off()

# ----------------------------- Tamanho do nicho ---------------------------
ray.size <- sapply(ray.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = 0.95)
})

tiff("grafico_boxplot_niche_size.tiff", width = 2000, height = 2000, res = 300)
boxplot(ray.size, col = clrs, pch = 16, cex = 0.5,
        ylab = "Niche Size", xlab = "Species",
        names = species.names)
dev.off()

# ----------------------------- Resumo da sobreposição ---------------------
summary(over.stat)







# ------------------------- Preparar dados ---------------------------------
ray <- ray %>%
  mutate(
    d13C = as.numeric(d13C),
    d15N = as.numeric(d15N)
  ) %>%
  select(Species, d13C, d15N)

# ------------------------- Configurações ----------------------------------
species.names <- c("Ebip", "Talb", "Asol")  # ordem fixa
clrs <- c(
  "Ebip" = "#ffb3ab",    # salmão
  "Talb" = "#7b9acc",    # azul escuro pastel
  "Asol" = "#a8e6cf"     # verde claro
)

# ------------------------- Parâmetros Bayesianos --------------------------
nsamples <- 1000
ray.par <- lapply(species.names, function(sp) {
  ii <- which(ray$Species == sp)
  niw.post(nsamples = nsamples, X = ray[ii, c("d13C","d15N")])
})
names(ray.par) <- species.names

# ------------------------- Tamanho do nicho -------------------------------
ray.size <- sapply(ray.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = 0.95)
}, simplify = "array")
dimnames(ray.size) <- list(NULL, species.names)

# Estatísticas resumidas do tamanho de nicho
size_summary <- lapply(species.names, function(sp) {
  vec <- ray.size[, sp]
  c(mean = mean(vec), median = median(vec), min = min(vec), max = max(vec))
})
names(size_summary) <- species.names

# ------------------------- Médias e covariâncias ---------------------------
# Média e covariância de exemplo (primeira amostra)
mu <- lapply(ray.par, function(x) x$Mu[,1])
sigma <- lapply(ray.par, function(x) x$Sigma[,,1])

# ------------------------- Sobreposição de nichos -------------------------
over.stat <- overlap(ray.par, nreps = nsamples, nprob = 1000, alpha = 0.95)

# Média da sobreposição entre espécies
over_mean <- apply(over.stat, c(1,2), mean)  # nspecies x nspecies

# Extrair uma sobreposição específica
over_Ebip_Talb <- over.stat["Ebip","Talb", ]
over_Talb_Asol <- over.stat["Talb","Asol", ]

# Resumo rápido da sobreposição
over_summary <- summary(as.vector(over.stat))

# ------------------------- Exemplo de saída --------------------------------
size_summary
mu
sigma
over_mean
summary(over_Ebip_Talb)
summary(over_Talb_Asol)
over_summary
