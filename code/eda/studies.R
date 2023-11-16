# studies = list(
#   Th0 = "eda_limma_thp_th0",
#   iTreg = "eda_limma_thp_itreg",
#   Th17 = "eda_limma_thp_th17",
#   Th2 = "limma_elo_thp_th2_lfc_cut",
#   Th1 = "limma_aijoe_thp_th1_lfc_cut"
# )

study.l = list(
  Th0 = list(
    dge.res = "limma_thp_th0",
    samples = list(
      HOURS_05h_vs_0h = c(9, 12),
      HOURS_1h_vs_0h = c(9, 12),
      HOURS_2h_vs_0h = c(12, 12),
      HOURS_4h_vs_0h = c(6, 12),
      HOURS_6h_vs_0h = c(12, 12),
      HOURS_12h_vs_0h = c(9, 12),
      HOURS_24h_vs_0h = c(12, 12),
      HOURS_48h_vs_0h = c(12, 12),
      HOURS_72h_vs_0h = c(9, 12)
    )
  ),
  iTreg = list(
    dge.res = "limma_thp_itreg",
    samples = list(
      HOURS_05h_vs_0h = c(6, 9),
      HOURS_1h_vs_0h = c(6, 9),
      HOURS_2h_vs_0h = c(9, 9),
      HOURS_4h_vs_0h = c(3, 9),
      HOURS_6h_vs_0h = c(9, 9),
      HOURS_12h_vs_0h = c(6, 9),
      HOURS_24h_vs_0h = c(9, 9),
      HOURS_48h_vs_0h = c(9, 9),
      HOURS_72h_vs_0h = c(6, 9)
    )
  ),
  Th17 = list(
    dge.res = "limma_thp_th17",
    samples = list(
      HOURS_05h_vs_0h = c(3, 3),
      HOURS_1h_vs_0h = c(3, 3),
      HOURS_2h_vs_0h = c(3, 3),
      HOURS_4h_vs_0h = c(3, 3),
      HOURS_6h_vs_0h = c(3, 3),
      HOURS_12h_vs_0h = c(3, 3),
      HOURS_24h_vs_0h = c(3, 3),
      HOURS_48h_vs_0h = c(3, 3),
      HOURS_72h_vs_0h = c(3, 3)
    )
  ),
  Th2 = list(
    dge.res = "limma_elo_thp_th2",
    samples = list(
      HOURS_05h_vs_0h = c(3, 3),
      HOURS_1h_vs_0h = c(3, 3),
      HOURS_2h_vs_0h = c(3, 3),
      HOURS_4h_vs_0h = c(3, 3),
      HOURS_6h_vs_0h = c(3, 3),
      HOURS_12h_vs_0h = c(3, 3),
      HOURS_24h_vs_0h = c(3, 3),
      HOURS_48h_vs_0h = c(3, 3),
      HOURS_72h_vs_0h = c(3, 3)
    )
  ),
  Th1 = list(
    dge.res = "limma_aijoe_thp_th1",
    samples = list(
      HOURS_12h_vs_0h = c(3, 3),
      HOURS_24h_vs_0h = c(3, 3),
      HOURS_48h_vs_0h = c(3, 3),
      HOURS_72h_vs_0h = c(2, 3)
    )
  ),
  Th0Cd2 = list(
    dge.res = "limma_izi_thp_th0_act",
    samples = list(
      HOURS_6h_vs_0h = c(4, 4),
      HOURS_12h_vs_0h = c(4, 4),
      HOURS_24h_vs_0h = c(4, 4),
      HOURS_48h_vs_0h = c(4, 4),
      HOURS_72h_vs_0h = c(4, 4)
    )
  ),
  Th0Cd4Mem = list(
    dge.res = "limma_arcelus_cd4_th0",
    samples = list(
      HOURS_2h_vs_0h = c(23, 24),
      HOURS_4h_vs_0h = c(24, 24),
      HOURS_8h_vs_0h = c(24, 24),
      HOURS_12h_vs_0h = c(19, 24),
      HOURS_24h_vs_0h = c(24, 24),
      HOURS_48h_vs_0h = c(22, 24),
      HOURS_72h_vs_0h = c(24, 24)
    )
  )
)

