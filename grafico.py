import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# Usar backend no interactivo (no abre ventanas)
matplotlib.use("Agg")

# Cargar CSV de Streptococcus
df = pd.read_csv("resultados streptococcus/streptococcus7001.csv")

# ==========================
# 1. Histograma de rho
# ==========================
plt.figure(figsize=(8, 5))
df["rho"].hist(bins=30, edgecolor="black")
plt.xlabel("rho (novedad genómica)")
plt.ylabel("Frecuencia")
plt.title("Distribución de rho - Streptococcus")
plt.tight_layout()
plt.savefig("histograma_rho_streptococcus7001.png", dpi=300)
plt.close()

# ==========================
# 2. Scatter |S| vs rho
# ==========================
plt.figure(figsize=(8, 5))
plt.scatter(df["|S|"], df["rho"], s=10)
plt.xlabel("|S| (tamaño estimado del sketch)")
plt.ylabel("rho (novedad genómica)")
plt.title("Novedad vs tamaño del sketch - Streptococcus")
plt.tight_layout()
plt.savefig("scatter_rho_streptococcus7001.png", dpi=300)
plt.close()

print("Gráficos guardados como histograma_rho_streptococcus6001.png y scatter_rho_streptococcus5001.png")
