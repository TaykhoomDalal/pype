import contextlib
import io
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import textalloc

from .statistics import multiple_testing_correction


REQUIRED_COLUMNS = {"category", "predictor", "pvalue", "beta"}


def _prepare_results(results):
	data = results.copy()
	missing = REQUIRED_COLUMNS - set(data.columns)
	if missing:
		raise ValueError("Missing result columns: " + ", ".join(sorted(missing)))

	data["pvalue"] = pd.to_numeric(data["pvalue"], errors="coerce")
	data["beta"] = pd.to_numeric(data["beta"], errors="coerce")
	if "negative_log10_pvalue" in data:
		data["negative_log10_pvalue"] = pd.to_numeric(
			data["negative_log10_pvalue"],
			errors="coerce",
		)
	else:
		with np.errstate(divide="ignore"):
			data["negative_log10_pvalue"] = -np.log10(data["pvalue"])

	data = data.dropna(
		subset=[
			"category",
			"predictor",
			"pvalue",
			"beta",
			"negative_log10_pvalue",
		]
	)
	if ((data["pvalue"] < 0) | (data["pvalue"] > 1)).any():
		raise ValueError("P-values must be in the interval [0, 1].")

	finite_values = data.loc[
		np.isfinite(data["negative_log10_pvalue"]),
		"negative_log10_pvalue",
	]
	maximum_plot_value = finite_values.max() + 1 if not finite_values.empty else 1
	data["_plot_pvalue"] = data["negative_log10_pvalue"].replace(
		[np.inf, -np.inf],
		[maximum_plot_value, 0],
	)
	return data


def _significance_level(data, alpha, correction):
	threshold = multiple_testing_correction(data["pvalue"], alpha, correction)
	return np.inf if threshold == 0 else -np.log10(threshold)


def _safe_name(value):
	name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
	return name or "plot"


def _label(row, annotation_width):
	value = row.get(
		"description",
		row.get("outcome", row["predictor"]),
	)
	return textwrap.fill(
		str(value).replace("_", " "),
		annotation_width,
		break_long_words=False,
	)


def _allocate_annotations(
	ax,
	annotation_rows,
	all_x,
	all_y,
	line_segments,
	directions,
	seed,
	annotation_width,
):
	if annotation_rows.empty:
		return

	x_lines = [
		np.asarray(x_coordinates)
		for x_coordinates, _ in line_segments
	]
	y_lines = [
		np.asarray(y_coordinates)
		for _, y_coordinates in line_segments
	]
	# textalloc prints when a label has no valid placement, even with verbose=False.
	with contextlib.redirect_stdout(io.StringIO()):
		textalloc.allocate(
			ax,
			annotation_rows["_annotation_x"].to_numpy(),
			annotation_rows["_annotation_y"].to_numpy(),
			[
				_label(row, annotation_width)
				for _, row in annotation_rows.iterrows()
			],
			x_scatter=np.asarray(all_x),
			y_scatter=np.asarray(all_y),
			x_lines=x_lines or None,
			y_lines=y_lines or None,
			textsize=7,
			margin=0.004,
			min_distance=0.01,
			max_distance=0.18,
			draw_lines=True,
			linecolor="0.35",
			draw_all=False,
			nbr_candidates=300,
			linewidth=0.4,
			direction=directions,
			priority_strategy=0 if seed is None else seed,
			avoid_label_lines_overlap=True,
			avoid_crossing_label_lines=True,
		)


def manhattan(results, output, alpha=0.05, correction="bonferroni", title=None, annotate=False, annotation_count=2, annotation_width=30, color_map="tab20", seed=None, width=15, height=8, dpi=300):
	"""Plot PheWAS significance grouped by category."""
	if annotation_count < 0 or annotation_width < 1:
		raise ValueError(
			"annotation_count must be non-negative and "
			"annotation_width must be positive."
		)
	data = _prepare_results(results)
	significance = _significance_level(data, alpha, correction)
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)

	significant = data[
		data["negative_log10_pvalue"] >= significance
	].sort_values(
		["category", "negative_log10_pvalue"],
		ascending=[True, False],
	)
	significant.drop(columns="_plot_pvalue").to_csv(
		output.with_name(output.stem + "_significant_results.tsv"),
		sep="\t",
		index=False,
	)

	order = (
		data.groupby("category")["_plot_pvalue"]
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
	data["_x"] = data["category"].map(positions).astype(float)
	data["_x"] += rng.uniform(-0.3, 0.3, len(data))
	is_significant = data["negative_log10_pvalue"] >= significance

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
			subset["_plot_pvalue"],
			c=[colors[category] for category in subset["category"]],
			marker=marker,
			s=18,
			alpha=0.75,
			linewidths=0.2,
			label=label,
		)

	annotation_rows = pd.DataFrame()
	if annotate:
		annotation_rows = (
			data[is_significant]
			.sort_values(
				["category", "negative_log10_pvalue"],
				ascending=[True, False],
			)
			.groupby("category", sort=False)
			.head(annotation_count)
			.assign(
				_annotation_x=lambda frame: frame["_x"],
				_annotation_y=lambda frame: frame["_plot_pvalue"],
			)
		)

	line_segments = []
	if np.isfinite(significance):
		ax.axhline(significance, color="red", linestyle="--", linewidth=0.7)
		line_segments.append(
			(ax.get_xlim(), (significance, significance))
		)
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
	_allocate_annotations(
		ax,
		annotation_rows,
		data["_x"],
		data["_plot_pvalue"],
		line_segments,
		["north"] * len(annotation_rows),
		seed,
		annotation_width,
	)
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
		data.assign(significant=data["negative_log10_pvalue"] >= significance)
		.groupby("category")
		.agg(
			total_count=("category", "size"),
			significant_count=("significant", "sum"),
		)
	)
	summary["percentage"] = (
		100 * summary["significant_count"] / summary["total_count"]
	)
	summary = summary.sort_values("percentage")

	fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
	bars = ax.bar(summary.index.astype(str), summary["percentage"], color="tab:blue")
	ax.bar_label(
		bars,
		labels=summary["total_count"].astype(str),
		padding=3,
		fontsize=8,
	)
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


def volcano(results, output, alpha=0.05, correction="bonferroni", title=None, annotate=False, annotation_count=5, annotation_width=30, width=15, height=8, dpi=300):
	"""Plot effect size against significance for each predictor."""
	if annotation_count < 0 or annotation_width < 1:
		raise ValueError(
			"annotation_count must be non-negative and "
			"annotation_width must be positive."
		)
	data = _prepare_results(results)
	significance = _significance_level(data, alpha, correction)
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)
	paths = []

	for predictor, predictor_data in data.groupby("predictor", sort=False):
		name = _safe_name(predictor)
		plot_path = output.with_name(output.stem + "_" + name + output.suffix)
		significant = predictor_data[
			predictor_data["negative_log10_pvalue"] >= significance
		].sort_values(
			"negative_log10_pvalue",
			ascending=False,
		)
		significant.drop(columns="_plot_pvalue").to_csv(
			output.with_name(output.stem + "_" + name + ".tsv"),
			sep="\t",
			index=False,
		)

		is_significant = (
			predictor_data["negative_log10_pvalue"] >= significance
		)
		fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
		ax.scatter(
			predictor_data.loc[~is_significant, "beta"],
			predictor_data.loc[~is_significant, "_plot_pvalue"],
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
					subset["_plot_pvalue"],
					color=color,
					marker=marker,
					s=24,
					alpha=0.75,
					linewidths=0.2,
					label=label,
				)

		annotation_rows = pd.DataFrame()
		if annotate:
			annotation_rows = (
				predictor_data[is_significant]
				.nlargest(annotation_count, "negative_log10_pvalue")
				.assign(
					_annotation_x=lambda frame: frame["beta"],
					_annotation_y=lambda frame: frame["_plot_pvalue"],
				)
			)

		line_segments = []
		if np.isfinite(significance):
			ax.axhline(significance, color="red", linestyle="--", linewidth=0.7)
			line_segments.append(
				(ax.get_xlim(), (significance, significance))
			)
		ax.axvline(
			0,
			color="grey",
			linestyle="--",
			linewidth=0.8,
			zorder=0,
		)
		line_segments.append(((0, 0), ax.get_ylim()))
		ax.set_xlabel("Effect size")
		ax.set_ylabel("-log10(p-value)")
		if title:
			ax.set_title(title + " : " + str(predictor))
		if ax.get_legend_handles_labels()[0]:
			ax.legend(loc="best")
		_allocate_annotations(
			ax,
			annotation_rows,
			predictor_data["beta"],
			predictor_data["_plot_pvalue"],
			line_segments,
			[
				"northwest" if beta < 0 else "northeast"
				for beta in annotation_rows.get("beta", [])
			],
			0,
			annotation_width,
		)
		fig.tight_layout()
		fig.savefig(plot_path, dpi=dpi, bbox_inches="tight")
		plt.close(fig)
		paths.append(plot_path)

	return paths
