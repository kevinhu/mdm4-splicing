# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import many
from scipy.stats import mannwhitneyu

merged_ccle_info = pd.read_csv(
    "../data/supplementary/S1_merged-ccle-info.txt", sep="\t", index_col=0
)

# %%
plt.figure(figsize=(4, 3))
ax = plt.subplot()

mut_hue = "#e23e57"
wt_hue = "#eaeaea"

sns.boxplot(
    x=merged_ccle_info["RPL22_MS_protein"],
    y=(merged_ccle_info["RPL22_mutation_classification_collapsed"] == "damaging").map(
        lambda x: {False: "WT/silent/non-conserving", True: "Damaging"}[x]
    ),
    notch=True,
    bootstrap=1000,
    palette={"WT/silent/non-conserving": wt_hue, "Damaging": mut_hue},
)
plt.ylabel("")
plt.xlabel("RPL22 protein")

damaging = merged_ccle_info[
    merged_ccle_info["RPL22_mutation_classification_collapsed"] == "damaging"
]["RPL22_MS_protein"].dropna()

wt = merged_ccle_info[
    merged_ccle_info["RPL22_mutation_classification_collapsed"] != "damaging"
]["RPL22_MS_protein"].dropna()

u, pval = mannwhitneyu(damaging, wt)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.spines["left"].set_position(("axes", -0.025))
ax.spines["bottom"].set_position(("axes", -0.025))

plt.text(0.05, 1.025, "P = " + many.visuals.as_si(pval, 2), transform=ax.transAxes)

# Truncate x range
plt.xlim(-3, 2)

plt.savefig("../figures/rpl22_protein_and_mut.pdf", bbox_inches="tight")
