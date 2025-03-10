#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python script.py <file_csv> <output_plot>")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_path = sys.argv[2]

    df = pd.read_csv(csv_file)
    x_values = df["mutation perc"].unique()
    x_ticks = [0, 1, 3, 5, 10, 20]

    beagle = df[df["tool"] == "beagle_1"]
    mupbwt = df[df["tool"] == "mupbwt"]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Wall clock time scatter plot
    # axes[0].scatter(
    #     mupbwt["mutation perc"],
    #     mupbwt["wall_clock"],
    #     label="μ-PBWT (wall clock)",
    #     marker="o",
    #     color="#5e81ac",
    # )
    axes[0].plot(
        mupbwt["mutation perc"],
        mupbwt["wall_clock"],
        label="μ-PBWT",
        marker="o",
        color="#5e81ac",
    )

    # axes[0].scatter(
    #     beagle["mutation perc"],
    #     beagle["wall_clock"],
    #     label="Beagle (wall clock)",
    #     marker="s",
    #     color="#bf616a",
    # )

    axes[0].plot(
        beagle["mutation perc"],
        beagle["wall_clock"],
        label="Beagle",
        marker="s",
        color="#bf616a",
    )
    if "beagle_312" in df["tool"].values:
        beagle = df[df["tool"] == "beagle_32"]
        axes[0].plot(
            beagle["mutation perc"],
            beagle["wall_clock"],
            label="Beagle (32 threads)",
            marker="s",
            color="#bf616a",
            linestyle="dashed",
        )

    if "beagle_11" in df["tool"].values:
        beagle = df[df["tool"] == "beagle_1"]
        axes[0].plot(
            beagle["mutation perc"],
            beagle["wall_clock"],
            label="Beagle (1 thread)",
            marker="s",
            color="#bf616a",
            linestyle="dotted",
        )

    beagle = df[df["tool"] == "beagle_1"]
    # axes[0].scatter(
    #     mupbwt["mutation perc"],
    #     mupbwt["user_time"] + mupbwt["sys_time"],
    #     label="μ-PBWT (user + sys)",
    #     marker="o",
    #     color="#5e81ac",
    # )
    # axes[0].plot(
    #     mupbwt["mutation perc"],
    #     mupbwt["user_time"] + mupbwt["sys_time"],
    #     label="μ-PBWT (user + sys)",
    #     marker="o",
    #     linestyle="dashed",
    #     color="#5e81ac",
    # )
    # axes[0].scatter(
    #     beagle["mutation perc"],
    #     beagle["user_time"] + beagle["sys_time"],
    #     color="#bf616a",
    # )

    # axes[0].plot(
    #     beagle["mutation perc"],
    #     beagle["user_time"] + beagle["sys_time"],
    #     label="Beagle (user + sys)",
    #     linestyle="dashed",
    #     marker="s",
    #     color="#bf616a",
    # )

    axes[0].set_xlabel("Mutation Percentage")
    axes[0].set_ylabel("Wall Clock Time (s)")
    axes[0].set_title("a) Execution Time")
    # axes[0].legend()
    # axes[0].grid()
    axes[0].set_ylim(0, 40)
    axes[0].set_xticks(x_ticks)
    axes[0].set_xticklabels(x_ticks)

    # Max memory scatter plot
    axes[1].plot(
        mupbwt["mutation perc"],
        mupbwt["max_mem"] / (1024 * 1024),
        label="μ-PBWT",
        marker="o",
        color="#5e81ac",
    )
    axes[1].plot(
        beagle["mutation perc"],
        beagle["max_mem"] / (1024 * 1024),
        label="Beagle",
        color="#bf616a",
        marker="s",
    )

    if "beagle_312" in df["tool"].values:
        beagle = df[df["tool"] == "beagle_32"]
        axes[1].plot(
            beagle["mutation perc"],
            beagle["max_mem"] / (1024 * 1024),
            label="Beagle (32 threads)",
            color="#bf616a",
            marker="s",
            linestyle="dashed",
        )
    if "beagle_11" in df["tool"].values:
        beagle = df[df["tool"] == "beagle_1"]
        axes[1].plot(
            beagle["mutation perc"],
            beagle["max_mem"] / (1024 * 1024),
            label="Beagle",
            color="#bf616a",
            marker="s",
            linestyle="dotted",
        )
    axes[1].set_xlabel("Mutation Percentage")
    axes[1].set_ylabel("Max Memory (GB)")
    axes[1].set_title("b) Max Memory Usage")
    # axes[1].legend()
    # axes[1].grid()

    # axes[1].semilogy()
    # axes[1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    # axes[1].set_ylim(0, 10)
    axes[1].set_xticks(x_ticks)
    axes[1].set_xticklabels(x_ticks)

    handles, labels = axes[0].get_legend_handles_labels()
    print(handles, labels)
    plt.tight_layout()

    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=4,
    )

    plt.savefig(sys.argv[2], dpi=500, bbox_inches="tight")
    # plt.show()
    plt.close()

    plt.show()
