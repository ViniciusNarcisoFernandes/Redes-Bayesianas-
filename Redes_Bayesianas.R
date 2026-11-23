library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(vcd)
library(bnlearn)
library(ggplot2)
library(patchwork)
set.seed(42)

# ==============================================================================
# CARREGAMENTO E ENGENHARIA DE VARIÁVEIS (TABELA 1)
# ==============================================================================
dados_brutos <- read_delim("Downloads/Dados Mpox/dados_mpox_limpos.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE,
                           show_col_types = FALSE)

# Garante nome limpo da primeira coluna 
names(dados_brutos)[1] <- "ano"

# Padronização de fatores e conversão obrigatória para base
dados_rede <- dados_brutos %>%
  mutate(
    ano = as.factor(as.character(ano)),
    regiao = as.factor(regiao),
    faixa_etaria = as.factor(faixa_etaria),
    sexo = as.factor(sexo),
    raca = as.factor(raca),
    ident_genero = as.factor(ident_genero),
    orienta_sexual = as.factor(orienta_sexual),
    hiv = as.factor(hiv),
    ist_ativa = as.factor(ist_ativa),
    pac_imunossup = factor(
      case_when(
        pac_imunossup %in% c(1, 2, 3) ~ "Yes",
        pac_imunossup == 4 ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("Yes", "No", "Unknown")
    )
  ) %>%
  select(ano, regiao, faixa_etaria, sexo, raca,
         ident_genero, orienta_sexual,
         hiv, ist_ativa, pac_imunossup) %>%
  droplevels() %>%
  as.data.frame() 

# Conferência de segurança dos níveis
cat("\n=== NÍVEIS POR VARIÁVEL ===\n")
print(sapply(dados_rede, nlevels))

# Frequências descritivas para conferência da Tabela 1
cat("\n=== RESUMO DESCRITIVO (TABELA 1) ===\n")
for (col in names(dados_rede)) {
  cat("\nVariável:", col, "\n")
  tab <- table(dados_rede[[col]])
  pct <- round(prop.table(tab) * 100, 1)
  print(cbind(n = tab, Percentual = pct))
}

# ==============================================================================
# TRIAGEM DE ASSOCIAÇÃO CATEGÓRICA: MATRIZ V DE CRAMÉR
# ==============================================================================
vars_cramer <- c("ano", "regiao", "faixa_etaria", "sexo", "raca",
                 "comp_sexual", "orienta_sexual", "ident_genero",
                 "hiv", "ist_ativa", "pac_imunossup")

df_cramer_base <- dados_brutos %>%
  mutate(
    pac_imunossup = case_when(
      pac_imunossup %in% c(1, 2, 3) ~ "Yes",
      pac_imunossup == 4 ~ "No",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(all_of(vars_cramer)) %>%
  as.data.frame()

matriz_v <- matrix(1, nrow = length(vars_cramer), ncol = length(vars_cramer),
                   dimnames = list(vars_cramer, vars_cramer))

for (i in 1:(length(vars_cramer) - 1)) {
  for (j in (i + 1):length(vars_cramer)) {
    tab_ij <- table(df_cramer_base[[i]], df_cramer_base[[j]])
    v_val <- assocstats(tab_ij)$cramer
    matriz_v[i, j] <- matriz_v[j, i] <- round(v_val, 3)
  }
}

cat("\n=== MATRIZ V DE CRAMÉR ===\n")
print(matriz_v)

# ==============================================================================
# APRENDIZAGEM ESTRUTURAL COM RESTRIÇÕES (BLACKLIST) E BOOTSTRAP (B = 500)
# ==============================================================================
demog  <- c("ano", "regiao", "faixa_etaria", "sexo", "raca")
behav  <- c("ident_genero", "orienta_sexual")
clin   <- c("hiv", "ist_ativa", "pac_imunossup")
raizes <- c("ano", "regiao", "faixa_etaria")

bl_1 <- expand.grid(from = c(clin, behav), to = demog, stringsAsFactors = FALSE)
bl_2 <- expand.grid(from = clin, to = behav, stringsAsFactors = FALSE)
bl_3 <- expand.grid(from = demog, to = raizes, stringsAsFactors = FALSE) %>%
  filter(from != to)

blacklist_mpox <- bind_rows(bl_1, bl_2, bl_3) %>%
  distinct() %>%
  as.data.frame()

cat("\nExecutando 500 reamostragens bootstrap (aguarde alguns segundos)...\n")
boot_rede <- boot.strength(
  data = dados_rede,
  R = 500,
  algorithm = "hc",
  algorithm.args = list(score = "bic", blacklist = blacklist_mpox, maxp = 3)
)

# Força e Direção dos arcos candidatos
tabela_bootstrap <- boot_rede %>%
  filter(strength >= 0.45) %>%
  arrange(desc(strength), desc(direction)) %>%
  mutate(
    Forca = round(strength, 3),
    Direcao = round(direction, 3),
    Status = ifelse(strength >= 0.70, "Incluído (>= 0,70)", "Descartado (< 0,70)")
  ) %>%
  select(from, to, Forca, Direcao, Status)

cat("\n=== TABELA DE ESTABILIDADE ESTRUTURAL (BOOTSTRAP B = 500) ===\n")
print(tabela_bootstrap, row.names = FALSE)

# Obter DAG consensual com limiar >= 0.70 e ajustar parâmetros Bayesianos (ISS = 10)
rede_consenso <- averaged.network(boot_rede, threshold = 0.70)
rede_ajustada <- bn.fit(rede_consenso, data = dados_rede, method = "bayes", iss = 10)

cat("\n=== MANTO DE MARKOV DO DESFECHO HIV ===\n")
print(mb(rede_consenso, node = "hiv"))

# ==============================================================================
# TABELAS CONDICIONAIS (CPTs) E INFERÊNCIA POR PERFIL (TABELA 4)
# ==============================================================================
tab_imuno_hiv <- prop.table(table(dados_rede$pac_imunossup, dados_rede$hiv), margin = 1) * 100
cat("\n=== TABELA 4: DISTRIBUIÇÃO CONDICIONAL IMUNOSSUPRESSÃO x HIV (%) ===\n")
print(round(tab_imuno_hiv, 1))

calc_perfil_hiv <- function(df_sub, nome_perfil) {
  d_valid <- df_sub %>% filter(hiv %in% c("Yes", "No"))
  n_tot <- nrow(d_valid)
  n_pos <- sum(d_valid$hiv == "Yes")
  teste <- prop.test(n_pos, n_tot)
  data.frame(
    Perfil = nome_perfil,
    n_validos = n_tot,
    Prob_HIV = round((n_pos / n_tot) * 100, 1),
    IC_inf = round(teste$conf.int[1] * 100, 1),
    IC_sup = round(teste$conf.int[2] * 100, 1)
  )
}

tabela_perfis <- bind_rows(
  calc_perfil_hiv(dados_rede, "Geral (Base Completa)"),
  calc_perfil_hiv(dados_rede %>% filter(sexo == "Female", orienta_sexual == "Heterosexual"),
                  "Feminino + Heterossexual"),
  calc_perfil_hiv(dados_rede %>% filter(sexo == "Male", orienta_sexual == "Heterosexual"),
                  "Masculino + Heterossexual"),
  calc_perfil_hiv(dados_rede %>% filter(sexo == "Male", orienta_sexual == "Bisexual"),
                  "Masculino + Bissexual"),
  calc_perfil_hiv(dados_rede %>% filter(sexo == "Male", orienta_sexual == "Homosexual"),
                  "Masculino + Homossexual"),
  calc_perfil_hiv(dados_rede %>% filter(sexo == "Male", orienta_sexual == "Homosexual", ist_ativa == "Yes"),
                  "Masc. + Homossexual + IST Ativa")
)

cat("\n=== INFERÊNCIA CONDICIONAL DE HIV POR PERFIL (SOROLOGIA INFORMADA) ===\n")
print(tabela_perfis, row.names = FALSE)

# ==============================================================================
# GERAÇÃO DO GRÁFICO MULTIPAINEL EM GGPLOT2 (FIGURA DO ARTIGO)
# ==============================================================================

tabela_perfis_en <- tabela_perfis %>%
  mutate(
    Perfil_EN = case_when(
      Perfil == "Geral (Base Completa)" ~ "Overall (Full Dataset)",
      Perfil == "Feminino + Heterossexual" ~ "Female + Heterosexual",
      Perfil == "Masculino + Heterossexual" ~ "Male + Heterosexual",
      Perfil == "Masculino + Bissexual" ~ "Male + Bisexual",
      Perfil == "Masculino + Homossexual" ~ "Male + Homosexual",
      Perfil == "Masc. + Homossexual + IST Ativa" ~ "Male + Homosexual + Active STI",
      TRUE ~ Perfil
    )
  )

df_plot_boot <- tabela_bootstrap %>%
  mutate(
    Arco = paste(from, "→", to),
    Forca_pct = Forca * 100,
    Status_EN = ifelse(Forca >= 0.70, "Included (≥ 0.70)", "Discarded (< 0.70)")
  )

gA <- ggplot(df_plot_boot, aes(x = reorder(Arco, Forca_pct), y = Forca_pct, fill = Status_EN)) +
  geom_col(width = 0.68, alpha = 0.92) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "#991B1B", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Forca_pct)), hjust = -0.15, size = 3.1, fontface = "bold", color = "#1E293B") +
  scale_fill_manual(
    values = c("Included (≥ 0.70)" = "#0F2942", "Discarded (< 0.70)" = "#94A3B8"),
    name = NULL
  ) +
  coord_flip() +
  ylim(0, 118) +
  labs(
    title = "A) Structural Arc Stability (B = 500 Bootstraps)",
    x = NULL,
    y = "Bootstrap Appearance Frequency (%)"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#0F2942", size = 11),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

# Panel B: Conditional HIV Probability by Epidemiological Profile
gB <- ggplot(tabela_perfis_en, aes(x = reorder(Perfil_EN, Prob_HIV), y = Prob_HIV)) +
  geom_col(fill = "#1D4ED8", width = 0.62, alpha = 0.88) +
  geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup), width = 0.22, linewidth = 0.75, color = "#0F172A") +
  geom_text(aes(y = IC_sup, label = sprintf("%.1f%%", Prob_HIV)), hjust = -0.2, size = 3.1, fontface = "bold", color = "#1E293B") +
  coord_flip() +
  ylim(0, 82) +
  labs(
    title = "B) Conditional HIV Positivity by Profile (95% CI)",
    x = NULL,
    y = "Conditional Probability of HIV Positive (%)"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#0F2942", size = 11),
    panel.grid.major.y = element_blank()
  )

# Panel C: Conditional Probability Table (CPT) — Immunosuppression vs. HIV
df_plot_cpt <- dados_rede %>%
  count(pac_imunossup, hiv) %>%
  group_by(pac_imunossup) %>%
  mutate(
    Total_grupo = sum(n),
    Pct = (n / Total_grupo) * 100,
    Grupo_label = case_when(
      pac_imunossup == "Yes" ~ "Immunosuppressed (Yes)\n[n = 8,089]",
      pac_imunossup == "No"  ~ "Non-Immunosuppressed (No)\n[n = 48,849]",
      TRUE                   ~ "Unknown Immunosuppression\n[n = 17,416]"
    )
  )

df_plot_cpt$Grupo_label <- factor(
  df_plot_cpt$Grupo_label,
  levels = c(
    "Immunosuppressed (Yes)\n[n = 8,089]",
    "Non-Immunosuppressed (No)\n[n = 48,849]",
    "Unknown Immunosuppression\n[n = 17,416]"
  )
)

gC <- ggplot(df_plot_cpt, aes(x = Grupo_label, y = Pct, fill = hiv)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, alpha = 0.92) +
  geom_text(aes(label = sprintf("%.1f%%", Pct)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3.2, fontface = "bold", color = "#1E293B") +
  scale_fill_manual(
    values = c("Yes" = "#991B1B", "No" = "#0F2942", "Unknown" = "#94A3B8"),
    labels = c("No" = "HIV Negative (No)", "Unknown" = "HIV Unknown", "Yes" = "HIV Positive (Yes)"),
    name = NULL
  ) +
  ylim(0, 102) +
  labs(
    title = "C) Conditional Distribution of HIV Serostatus by Immunosuppression Status (Learned CPT)",
    x = NULL,
    y = "Conditional Probability (%)"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#0F2942", size = 11),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# Combine all 3 panels and export at 300 DPI
painel_final <- (gA | gB) / gC +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(
    title = "Discrete Bayesian Network Modeling Results — Mpox Notifications in Brazil (N = 74,354)",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", color = "#0F2942", hjust = 0.5)
    )
  )

print(painel_final)
