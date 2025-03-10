#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1])
print(df)
tot_map = {0: 60604, 1: 60862, 3: 61362, 5: 61995, 10: 63153, 20: 65265}
df["map_value"] = df["mutation perc"].apply(lambda x: tot_map.get(int(x), None))
print(df)
print(df["map_value"])

df["phased_ratio"] = df["phased count"] / df["map_value"] * 100

x_values = df["mutation perc"].unique()
x_ticks = [0, 1, 3, 5, 10, 20]
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

beagle = df[df["tool"] == "beagle"]
mupbwt = df[df["tool"] == "mupbwt"]
axes[0].plot(
    mupbwt["mutation perc"],
    mupbwt["switch error rate"],
    label="μ-PBWT",
    marker="o",
    color="#5e81ac",
)
axes[0].plot(
    beagle["mutation perc"],
    beagle["switch error rate"],
    label="Beagle",
    marker="s",
    color="#bf616a",
)
axes[0].set_title("a) Switch Error Rate")
axes[0].set_ylabel("Error Rate")
axes[0].set_xlabel("Mutation %")
# axes[0].legend()
axes[0].set_xticks(x_ticks)
axes[0].set_xticklabels(x_ticks)

axes[1].plot(
    mupbwt["mutation perc"],
    mupbwt["mismatch rate"],
    label="μ-PBWT",
    marker="o",
    color="#5e81ac",
)
axes[1].plot(
    beagle["mutation perc"],
    beagle["mismatch rate"],
    label="Beagle",
    marker="s",
    color="#bf616a",
)
axes[1].set_title("b) Mismatch Rate")
axes[1].set_ylabel("Error Rate")
axes[1].set_xlabel("Mutation %")
# axes[1].legend()
axes[1].set_xticks(x_ticks)
axes[1].set_xticklabels(x_ticks)

axes[2].plot(
    mupbwt["mutation perc"],
    mupbwt["phased_ratio"],
    label="μ-PBWT",
    marker="o",
    color="#5e81ac",
)
axes[2].plot(
    beagle["mutation perc"],
    beagle["phased_ratio"],
    label="Beagle",
    marker="s",
    color="#bf616a",
)
axes[2].set_title("c) % of Heterozygous Phased Sites")
axes[2].set_ylabel("Ratio")
axes[2].set_xlabel("Mutation %")
# axes[2].legend()
axes[2].set_xticks(x_ticks)
axes[2].set_xticklabels(x_ticks)
axes[2].set_ylim(0, 105)

handles, labels = axes[0].get_legend_handles_labels()
print(handles, labels)
plt.tight_layout()

fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.05),
    ncol=2,
)

# Aggiusta i margini per lasciare spazio alla legenda

plt.savefig(sys.argv[2], dpi=500, bbox_inches="tight")
# plt.show()
plt.close()
