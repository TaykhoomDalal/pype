import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

from .statistics import multiple_testing_correction


REQUIRED_COLUMNS = {"Category", "Independent_Var", "p-val", "beta"}


def _prepare_results(results):
	data = results.copy()
	missing = REQUIRED_COLUMNS - set(data.columns)
	if missing:
		raise ValueError("Missing result columns: " + ", ".join(sorted(missing)))

	data["p-val"] = pd.to_numeric(data["p-val"], errors="coerce")
	data["beta"] = pd.to_numeric(data["beta"], errors="coerce")
	if "-log(p)" in data:
		data["-log(p)"] = pd.to_numeric(data["-log(p)"], errors="coerce")
	else:
		with np.errstate(divide="ignore"):
			data["-log(p)"] = -np.log10(data["p-val"])

	data = data.dropna(subset=["Category", "Independent_Var", "p-val", "beta", "-log(p)"])
	if ((data["p-val"] < 0) | (data["p-val"] > 1)).any():
		raise ValueError("P-values must be in the interval [0, 1].")

	finite = data.loc[np.isfinite(data["-log(p)"]), "-log(p)"]
	cap = finite.max() + 1 if not finite.empty else 1
	data["_plot_logp"] = data["-log(p)"].replace([np.inf, -np.inf], [cap, 0])
	return data


def _significance_level(data, alpha, correction):
	threshold = multiple_testing_correction(data["p-val"], alpha, correction)
	return np.inf if threshold == 0 else -np.log10(threshold)


def _safe_name(value):
	name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
	return name or "plot"


def _label(row):
	value = row.get("Description", row.get("Data_Field", row["Independent_Var"]))
	return textwrap.fill(str(value).replace("_", " "), 30)


def manhattan(results, output, alpha=0.05, correction="bonferroni", title=None, annotate=False, top=2, color_map="tab20", seed=None, width=15, height=8, dpi=300):
	"""Plot PheWAS significance grouped by category."""
	data = _prepare_results(results)
	significance = _significance_level(data, alpha, correction)
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)

	significant = data[data["-log(p)"] >= significance].sort_values(
		["Category", "-log(p)"],
		ascending=[True, False],
	)
	significant.drop(columns="_plot_logp").to_csv(
		output.with_name(output.stem + "_significant_results.tsv"),
		sep="\t",
		index=False,
	)

	order = (
		data.groupby("Category")["_plot_logp"]
		.mean()
		.sort_values()
		.index.tolist()
	)
	positions = {category: index for index, category in enumerate(order)}
	cmap = plt.get_cmap(color_map)
	colors = {
		category: cmap(index / max(1, len(order) - 1))
		for index, category in enumerate(order)
	}

	rng = np.random.default_rng(seed)
	data["_x"] = data["Category"].map(positions).astype(float)
	data["_x"] += rng.uniform(-0.3, 0.3, len(data))
	is_significant = data["-log(p)"] >= significance

	fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
	for mask, marker, label in [
		(~is_significant, "o", "Not significant"),
		(is_significant & (data["beta"] >= 0), "^", "Positive beta"),
		(is_significant & (data["beta"] < 0), "v", "Negative beta"),
	]:
		subset = data[mask]
		if subset.empty:
			continue
		ax.scatter(
			subset["_x"],
			subset["_plot_logp"],
			c=[colors[category] for category in subset["Category"]],
			marker=marker,
			s=18,
			alpha=0.75,
			linewidths=0.2,
			label=label,
		)

	if annotate:
		texts = []
		top_rows = (
			data[is_significant]
			.sort_values(["Category", "-log(p)"], ascending=[True, False])
			.groupby("Category", sort=False)
			.head(top)
		)
		for _, row in top_rows.iterrows():
			texts.append(ax.text(row["_x"], row["_plot_logp"], _label(row), fontsize=7))
		if texts:
			adjust_text(texts, ax=ax, arrowprops={"arrowstyle": "->", "lw": 0.3})

	if np.isfinite(significance):
		ax.axhline(significance, color="red", linestyle="--", linewidth=0.7)
	ax.set_xlabel("Phenotype category")
	ax.set_ylabel("-log10(p-value)")
	ax.set_xticks(range(len(order)))
	ax.set_xticklabels(
		[textwrap.fill(str(category), 15) for category in order],
		rotation=90 if len(order) > 5 else 0,
	)
	if title:
		ax.set_title(title)
	if ax.get_legend_handles_labels()[0]:
		ax.legend(loc="best")
	fig.tight_layout()
	fig.savefig(output, dpi=dpi, bbox_inches="tight")
	plt.close(fig)
	return output


def category_enrichment(results, output, alpha=0.05, correction="bonferroni", title=None, width=15, height=8, dpi=300):
	"""Plot the percentage of significant associations in each category."""
	data = _prepare_results(results)
	significance = _significance_level(data, alpha, correction)
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)

	summary = (
		data.assign(Significant=data["-log(p)"] >= significance)
		.groupby("Category")
		.agg(total=("Category", "size"), significant=("Significant", "sum"))
	)
	summary["percentage"] = 100 * summary["significant"] / summary["total"]
	summary = summary.sort_values("percentage")

	fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
	bars = ax.bar(summary.index.astype(str), summary["percentage"], color="tab:blue")
	ax.bar_label(bars, labels=summary["total"].astype(str), padding=3, fontsize=8)
	ax.set_xlabel("Phenotype category")
	ax.set_ylabel("Significant associations (%)")
	ax.set_ylim(0, max(1, summary["percentage"].max() * 1.15 if not summary.empty else 1))
	ax.set_xticks(range(len(summary)))
	ax.set_xticklabels(
		[textwrap.fill(str(category), 15) for category in summary.index],
		rotation=90 if len(summary) > 5 else 0,
	)
	if title:
		ax.set_title(title)
	fig.tight_layout()
	fig.savefig(output, dpi=dpi, bbox_inches="tight")
	plt.close(fig)
	return output


def volcano(results, output, alpha=0.05, correction="bonferroni", title=None, annotate=False, top=5, width=15, height=8, dpi=300):
	"""Plot effect size against significance for each predictor."""
	data = _prepare_results(results)
	significance = _significance_level(data, alpha, correction)
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)
	paths = []

	for predictor, predictor_data in data.groupby("Independent_Var", sort=False):
		name = _safe_name(predictor)
		plot_path = output.with_name(output.stem + "_" + name + output.suffix)
		significant = predictor_data[predictor_data["-log(p)"] >= significance].sort_values(
			"-log(p)",
			ascending=False,
		)
		significant.drop(columns="_plot_logp").to_csv(
			output.with_name(output.stem + "_" + name + ".tsv"),
			sep="\t",
			index=False,
		)

		is_significant = predictor_data["-log(p)"] >= significance
		fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
		ax.scatter(
			predictor_data.loc[~is_significant, "beta"],
			predictor_data.loc[~is_significant, "_plot_logp"],
			color="black",
			s=14,
			alpha=0.65,
			linewidths=0.2,
		)
		for mask, marker, color, label in [
			(is_significant & (predictor_data["beta"] >= 0), "^", "tab:red", "Positive beta"),
			(is_significant & (predictor_data["beta"] < 0), "v", "tab:blue", "Negative beta"),
		]:
			subset = predictor_data[mask]
			if not subset.empty:
				ax.scatter(
					subset["beta"],
					subset["_plot_logp"],
					color=color,
					marker=marker,
					s=24,
					alpha=0.75,
					linewidths=0.2,
					label=label,
				)

		if annotate:
			texts = [
				ax.text(row["beta"], row["_plot_logp"], _label(row), fontsize=7)
				for _, row in predictor_data[is_significant].nlargest(top, "-log(p)").iterrows()
			]
			if texts:
				adjust_text(texts, ax=ax, arrowprops={"arrowstyle": "->", "lw": 0.3})

		if np.isfinite(significance):
			ax.axhline(significance, color="red", linestyle="--", linewidth=0.7)
		ax.axvline(0, color="grey", linewidth=0.5)
		ax.set_xlabel("Effect size")
		ax.set_ylabel("-log10(p-value)")
		if title:
			ax.set_title(title + " : " + str(predictor))
		if ax.get_legend_handles_labels()[0]:
			ax.legend(loc="best")
		fig.tight_layout()
		fig.savefig(plot_path, dpi=dpi, bbox_inches="tight")
		plt.close(fig)
		paths.append(plot_path)

	return paths
