#%%
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from matplotlib.cm import ScalarMappable
import pandas as pd
import numpy as np

plt.rc("font", size=11)
plt.rc("axes", labelsize="large")

df = pd.read_table(
    "../../build/example_programs/results/2Higgses_pdf2.dat",
    sep=r"\s+",
)

df["deltaChisq"] = df.chisq - np.min(df.chisq)
df["deltaChisqMu"] = df.chisq_mu - np.min(df.chisq_mu)
df["chisqMbar"] = df.chisq_mh - df.chisq_sep
norm = mcol.Normalize(0, 6)
cmap = "YlOrRd"
fig, axes = plt.subplots(figsize=(6, 5.1), ncols=2, nrows=2, sharex=True, sharey=True)
ax = axes[0][0]
ax.set_title(r"$\Delta\chi^2_{\overline{m}}$")
ax.tricontourf(
    df.Mh1,
    df.Mh2,
    df.chisqMbar,
    norm=norm,
    levels=np.linspace(0, 6),
    extend="max",
    cmap=cmap,
)

ax = axes[0][1]
ax.set_title(r"$\Delta\chi^2_\mathrm{sep}$")
ax.tricontourf(
    df.Mh1,
    df.Mh2,
    df.chisq_sep,
    norm=norm,
    levels=np.linspace(0, 6),
    extend="max",
    cmap=cmap,
)


ax = axes[1][0]
ax.set_title(r"$\Delta\chi^2_\mu$")
ax.tricontourf(
    df.Mh1,
    df.Mh2,
    df.deltaChisqMu,
    norm=norm,
    levels=np.linspace(0, 6),
    extend="max",
    cmap=cmap,
)

ax = axes[0][1]
ax.set_title(r"$\chi^2_\mathrm{sep}$")
ax.tricontourf(
    df.Mh1,
    df.Mh2,
    df.chisq_sep,
    norm=norm,
    levels=np.linspace(0, 6),
    extend="max",
    cmap=cmap,
)

ax = axes[1][1]
ax.set_title(r"$\Delta\chi^2$")
ax.tricontourf(
    df.Mh1,
    df.Mh2,
    df.deltaChisq,
    norm=norm,
    levels=np.linspace(0, 6),
    extend="max",
    cmap=cmap,
)
ax.tricontour(
    df.Mh1,
    df.Mh2,
    df.deltaChisq,
    levels=[2.3, 5.99],
    colors=["k", "k"],
    linestyles=["-", "--"],
)

for ax in axes.flatten():
    ax.set_xlim(123.5, 126.5)
    ax.set_ylim(123.5, 126.5)
    ax.set_aspect("equal", "box")

for ax in axes[1]:
    ax.set_xlabel(r"$M_1$ [GeV]")
for ax in axes[:,0]:
    ax.set_ylabel(r"$M_2$ [GeV]")

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
fig.colorbar(sm, ax=axes, label=r"$\chi^2$", extend="max",shrink=0.7)
fig.savefig(
    "/home/jonasw/Dropbox/HiggsSignals/HS-2-documentation/draft/figs/2Higgses.pdf",
    bbox_inches="tight",
)
# %%
